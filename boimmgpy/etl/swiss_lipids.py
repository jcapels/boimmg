import gzip
import io
import pandas as pd
import requests
from joblib import Parallel, delayed
from neo4j import GraphDatabase
from tqdm import tqdm

from boimmg.boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor,set_database_information

log,user,password = CompoundsDBAccessor.read_config_file()

class SwissLipidsExtractor:
    """
    Class to extract information from swiss lipids
    """

    def extract(self) -> pd.DataFrame:
        """
        This class calls the scrape_data method and creates a pandas dataframe with the scraped data returned from
            scrape_data method.
        :return: Dataframe of the whole Swiss Lipids data
        :rtype: pd.DataFrame
        """
        return self._extract(self.scrape_data())

    @staticmethod
    def scrape_data():
        """
        This class downloads the ZIP file and extracts the CSV files of SwissLipids.

        :return: csv file of the whole Swiss Lipids data
        :rtype: csv file
        """
        get_zip = requests.get("https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv")
        file_unzip = gzip.open(io.BytesIO(get_zip.content), 'rb')
        return file_unzip

    @staticmethod
    def _extract(raw_file) -> pd.DataFrame:
        """
        Method to create a pandas dataframe with the scraped data.
        raw_file: csv file of the whole Swiss Lipids database
        :type raw_file: csv file
        :return: Data frame of the whole Swiss Lipids database
        :rtype: pd.DataFrame
        """
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe


class SwissLipidsTransformer:
    """
    Class to transform the lipid maps dataframe
    """

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        This method allows to transform the dataframe previously extracted into a desired data structure
        :param df:  Whole pandas dataframe, previously extracted in the extract class
        :type df: pd.Dataframe
        :return: treated dataframe of Swiss Lipids, only with two columns, synonym and ID
        :rtype: pd.DataFrame
        """

        iteration = len(df)
        parallel_callback = Parallel(8)
        data_treated = parallel_callback(delayed(self.treat)(df.iloc[[i]]) for i in tqdm(range(iteration)))
        data_treated = pd.concat(data_treated)
        return data_treated

    @staticmethod
    def treat(data: pd.DataFrame) -> pd.DataFrame:
        """
        This method uses the whole Swiss Lipids dataframe and creates another with only two desirable columns,
        ID and Synonym, the last one encompasses the abbreviation too :param data:  Whole pandas dataframe,
        previously extracted in the extract class
        :type data: pd.DataFrame
        :return: Dataframe with two columns of the initial data, ID and Synonym
        :rtype: pd.DataFrame
        """
        new_df = pd.DataFrame(columns=['Lipid ID', 'Synonym'])
        counter = 0
        for i, row in data.iterrows():
            lipid_id = row["Lipid ID"]
            abreviation = row["Abbreviation*"]
            synonyms = row["Synonyms*"]
            if abreviation is not None and not pd.isnull(abreviation):
                abbreviation_splits = abreviation.split('|')
                for split in abbreviation_splits:
                    new_df.at[counter, "Lipid ID"] = lipid_id
                    split = split.replace(" ", "")
                    new_df.at[counter, "Synonym"] = split.lower()
                    counter += 1

            if synonyms is not None and not pd.isnull(synonyms):
                synonyms_splits = synonyms.split('|')
                for split in synonyms_splits:
                    new_df.at[counter, "Lipid ID"] = lipid_id
                    split = split.replace(" ", "")
                    new_df.at[counter, "Synonym"] = split.lower()
                    counter += 1
        return new_df


class SwissLipidsLoader:
    """
    Class that loads the treated data into the database
    """

    def load(self, treated_df: pd.DataFrame):
        """
        This method connects and loads the treated data into the database
        :param treated_df: Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
        :type treated_df: pd.DataFrame
        """
        # self.set_synonym(self.get_connection_list(treated_df))
        self.load_multiprocessing(treated_df)

    @staticmethod
    def load_multiprocessing(df: pd.DataFrame):
        itera = len(df)
        parallel_callback = Parallel(8)
        parallel_callback(delayed(get_connection_list)(df.iloc[[i]]) for i in tqdm(range(itera)))


data_base_connection = GraphDatabase.driver(uri=log,
                                            auth=(user, password))


def get_connection_list(df: pd.DataFrame):
    """
    This method creates the querys necessary to upload the treated data into the database
    :param df:  Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
    :type df: pd.DataFrame
    :return: List of querys necessary to the upload of the whole dataframe
    :rtype: list
    """
    with data_base_connection.session() as session1:
        for i, row in df.iterrows():
            swiss_lipids_id = row["Lipid ID"]
            sl_synonym = row["Synonym"]
            session1.run('MERGE (s: Synonym {synonym:"%s"})' % str(sl_synonym))
            session1.run("match (l:SwissLipidsCompound),(s:Synonym) where l.swiss_lipids_id=$sl_id and "
                         "s.synonym=$synonym merge (s)-[:is_synonym_of]->(l)", synonym=sl_synonym,
                         sl_id=swiss_lipids_id)


"""
dag=DAG(dag_id="dag_etl_lm",schedule_interval="@once",
        start_date=datetime(2022, 1, 1),
        catchup=False,
        )

def extract(**kwargs):
    
    Function where the extract method will be called.
    :param kwargs:
    
    ti = kwargs['ti']
    extractor=SwissLipidsExtractor()
    raw_df=extractor.extract()
    ti.xcom_push('order_data', raw_df)

def transform(**kwargs):

    Function where the transform method will be called.
    :param kwargs:
    
    ti = kwargs['ti']
    extract_df = ti.xcom_pull(task_ids='extract', key='order_data')
    transformer=SwissLipidsTransformer()
    treated_data=transformer.transform(extract_df)
    

def load(**kwargs):
    
    Function where the load method will be called.
    :param kwargs:
    
    ti = kwargs['ti']
    transformed_df = ti.xcom_pull(task_ids='transform', key='order_data')
    loader = SwissLipidsLoader()
    loader.load(transformed_df)          
    
with dag:
    extract_task = PythonOperator(
        task_id='extract',
        dag=dag,
        python_callable=extract,
    )
    transform_task = PythonOperator(
        task_id='transform',
        dag=dag,
        python_callable=transform,
    )
    load_task = PythonOperator(
        task_id='load',
        dag=dag,
        python_callable=load,
    )
    extract_task >> transform_task >> load_task
    """
