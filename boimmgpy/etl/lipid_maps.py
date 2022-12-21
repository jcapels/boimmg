import io
import pandas as pd
import requests
import zipfile
from joblib import Parallel, delayed
from rdkit.Chem import PandasTools
from tqdm import tqdm

from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager
from boimmgpy.etl._utils import insert_in_database_lipid_maps


class LipidMapsExtractor:
    """
    Class to extract information from lipid maps
    """

    def extract(self) -> pd.DataFrame:
        """
        This class calls the scrape_data method and creates a pandas dataframe with the scraped data returned from
        scrape_data method.
        :return: Dataframe of the whole Lipid Maps data
        :rtype: pd.DataFrame
        """
        return self._extract(self.scrape_data())

    @staticmethod
    def _extract(raw) -> pd.DataFrame:
        """
        Method to create a pandas dataframe with the scraped data.

        raw: sdf file of the whole lipid maps database
        :type raw: Sdf file
        :return: Data frame of the whole Lipid Maps database
        :rtype: pd.DataFrame
        """

        raw_lipid_maps_data = PandasTools.LoadSDF(raw)
        lm_dataframe = pd.DataFrame(raw_lipid_maps_data)
        return lm_dataframe

    @staticmethod
    def scrape_data():
        """
        This method downloads the ZIP file and extracts the sdf file of lipid maps.
        :return: sdf file of the whole lipid maps database
        :rtype: Sdf file
        """
        raw_file = requests.get("https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip")
        file_unziped = zipfile.ZipFile(io.BytesIO(raw_file.content))
        return file_unziped.open('structures.sdf')


class LipidMapsTransformer:
    """Class to transform the lipid maps dataframe
    """

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """

        This method allows to transform the dataframe previously extracted into a desired data structure
        :param df:  Whole pandas dataframe, previosly extracted in the extract class
        :type df: pd.Dataframe
        :return: treated dataframe of lipid maps, only with two columns, synonym and ID
        :rtype: pd.DataFrame
        """
        iteration = len(df)
        parallel_callback = Parallel(8)
        data_treated = parallel_callback(delayed(self.treat_lm_dataframe)(df.iloc[[i]]) for i in tqdm(range(iteration)))
        data_treated = pd.concat(data_treated)
        return data_treated

    @staticmethod
    def treat_lm_dataframe(lm_dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        This method uses the whole lipid maps dataframe and creates another with only two desirable columns,
        ID and Synonym, the last one encompasses the abbreviation too :param lm_dataframe:  Whole pandas dataframe,
        previously extracted in the extract class :type lm_dataframe: pd.DataFrame :return: Dataframe with two
        columns of the initial data, ID and Synonym
        :rtype: pd.DataFrame
        """
        new_df = pd.DataFrame(columns=['LM_ID', 'SYNONYMS'])
        counter = 0
        for i, row in lm_dataframe.iterrows():
            lipid_id = row["LM_ID"]
            abbreviation = row["ABBREVIATION"]
            synonyms = row["SYNONYMS"]
            if abbreviation is not None and not pd.isnull(abbreviation):
                abbreviation_splits = abbreviation.split(';')
                for split in abbreviation_splits:
                    new_df.at[counter, "LM_ID"] = lipid_id
                    split = split.replace(" ", "")
                    new_df.at[counter, "SYNONYMS"] = split.lower()
                    counter += 1

            if synonyms is not None and not pd.isnull(synonyms):
                synonyms_splits = synonyms.split(';')
                for split in synonyms_splits:
                    new_df.at[counter, "LM_ID"] = lipid_id
                    split = split.replace(" ", "")
                    new_df.at[counter, "SYNONYMS"] = split.lower()
                    counter += 1
        return new_df


class LipidMapsLoader:
    """
    Class that loads the treated data into the database
    """

    def __init__(self):
        self.driver = DatabaseAccessManager().connect()

    def load(self, df: pd.DataFrame):
        """
        This method connects and loads the treated data into the database
        :param df: Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be
            load
        :type df: pd.DataFrame
        """
        self.load_multiprocessing(df)

    @staticmethod
    def load_multiprocessing(df: pd.DataFrame, n_jobs: int = 8):
        n_iterations = len(df)
        parallel_callback = Parallel(n_jobs)
        parallel_callback(delayed(insert_in_database_lipid_maps)(df.iloc[[i]]) for i in tqdm(range(n_iterations)))

    def insert_in_database_lipid_maps(self, df: pd.Series):
        """
        This method creates the queries necessary to upload the treated data into the database
        :param df:  Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
        :type df: pd.DataFrame
        :return: List of queries necessary to the upload of the whole dataframe
        :rtype: list
        """

        with self.driver.session() as session:
            for i, row in df.iterrows():
                lipid_maps_id = row["LM_ID"]
                lm_synonym = row["SYNONYMS"]
                session.run('MERGE (s: Synonym {synonym:"%s"})' % str(lm_synonym))
                session.run(
                    "match (l:LipidMapsCompound),(s:Synonym) where l.lipidmaps_id=$lipid_maps_id and s.synonym=$synonym "
                    "merge (s)-[:is_synonym_of]->(l)",
                    synonym=lm_synonym, lipid_maps_id=lipid_maps_id)


"""
dag=DAG(dag_id="dag_etl_lm",schedule_interval="@once",
        start_date=datetime(2022, 1, 1),
        catchup=False,
        )


def extract(**kwargs):

    Function where the extract method will be called.
    :param kwargs:
   
    ti = kwargs['ti']
    extractor=LipidMapsExtractor()
    raw_df=extractor.extract()
    ti.xcom_push('order_data', raw_df)

def transform(**kwargs):
   
    Function where the load transform will be called.
    :param kwargs:
    
    ti = kwargs['ti']
    extract_df = ti.xcom_pull(task_ids='extract', key='order_data')
    transformer=LipidMapsTransformer()
    treated_data=transformer.transform(extract_df)
    ti.xcom_push('order_data', treated_data)
    

def load(**kwargs):
    
    Function where the load method will be called.
    :param kwargs:
    
    ti = kwargs['ti']
    transformed_df = ti.xcom_pull(task_ids='transform', key='order_data')
    loader = LipidMapsLoader()
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
