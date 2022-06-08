import pandas as pd
import pendulum
from airflow.decorators import task
from neo4j import GraphDatabase
from rdkit.Chem import PandasTools
from airflow.models.dag import dag
from airflow.operators.python import PythonOperator
import requests, zipfile, io, gzip
from boimmgpy.etl.airflow_interfaces import AirflowExtractor, AirflowTransformer, AirflowLoader, AirflowPipeline


class SwissLipidsExtractor(AirflowExtractor):
    """
    Class to extract information from lipid maps
    """

    def extract(self)->pd.DataFrame:
        """
        This class calls the scrape_data method and creates a pandas dataframe with the scraped data.
        :return: pandas data frame of lipid maps data 
        """
        return self._extract(self.scrape_data())


    
    def scrape_data(self):
        """This class downloads the ZIP file and extracts the CSV files of SwissLipids.

        Returns:
            CSV file: Csv file with swisss lipids data
        """
        get_zip = requests.get("https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv")
        raw_file = gzip.open(io.BytesIO(get_zip.content),'rb')
        return raw_file
        
        


    def _extract(self,raw_file) -> pd.DataFrame:
        """Method to create a pandas dataframe with the scraped data.

        Args:
            raw_file (csv file):  Csv file of the swiss lipids data

        Returns:
            pd.DataFrame: Swiss lipids data tranformed into a dataframe structure
        """
        table= pd.read_table(raw_file,engine='python',encoding='ISO-8859-1')
        df=pd.DataFrame(table)
        return df

    

class SwissLipidsTransformer(AirflowTransformer):
    """
    Class to transform the lipid maps dataframe
    """
    
    def transform(self,data: pd.DataFrame)->pd.DataFrame:
        """This method allows to transform the dataframe previously extracted into a desired datastructure

        Args:
            data (pd.DataFrame): Whole dataframe of swiss lipids
        Returns:
            pd.DataFrame: Treated dataframe
        """
        data_treated=self.treat(data)
        return data_treated
        
    def treat(self,data: pd.DataFrame)->pd.DataFrame:
        """Method thad treats the dataframe of swiss lipids and creates another with only two columns, one for lipid id and another 
        for synonyms and abbrebiations 

        Args:
            data (pd.DataFrame): whole dataframe of swiss lipids data

        Returns:
            pd.DataFrame: treated dataframe, 2 colums (id and synonym), of the swiss lipids data
        """
        new_df=pd.DataFrame(columns=['Lipid ID','Synonym'])
        counter=0
        for i, row in data.iloc[:100].iterrows():
            lipid_id = row["Lipid ID"]
            abreviation = row["Abbreviation*"]
            synonyms = row["Synonyms*"]
            if abreviation is not None and not pd.isnull(abreviation):
                abbreviation_splits = abreviation.split('|')
                for split in abbreviation_splits:
                    new_df.at[counter, "Lipid ID"] = lipid_id
                    new_df.at[counter, "Synonym"] = split
                    counter+=1

            if synonyms is not None and not pd.isnull(synonyms):
                synonyms_splits = synonyms.split('|')
                for split in synonyms_splits:
                    new_df.at[counter, "Lipid ID"] = lipid_id
                    new_df.at[counter, "Synonym"] = split
                    counter+=1
        return new_df


class SwissLipidsLoader(AirflowLoader):
    """
    Class that loads the treated data into the database
    """

    def load(self, treated_df: pd.DataFrame):
        """Method that loads the treated dataframe into our database

        Args:
            df (pd.DataFrame): treated dataframe of swiss lipids
        """
        self.set_synonym(self.get_connection_list(treated_df))

    def set_synonym(self, connections:list):
        data_base_connection = GraphDatabase.driver(uri="bolt://localhost:7687",auth=("neo4j","potassio19"))
        session = data_base_connection.session()
        for i in connections:
            session.run(i)

    def get_connection_list(self,df : pd.DataFrame)->list:
        """This method creates a list of querys necessary to do the upload of the dataframe in our database

        Args:
            df (pd.DataFrame): treated dataframe

        Returns:
            list: Querys list that will be use in the upload into our database
        """
        creation_nodes_list=[]
        for i,row in df.iterrows():
            swiss_lipids_id=row["Lipid ID"]
            sl_synonym=row["Synonym"]
            creat_node_connection='MATCH (u:SwissLipidsCompound) WHERE u.swiss_lipids_id="' + str(swiss_lipids_id) + '" merge (s: Synonym {swiss_lipids_id : "' + str(swiss_lipids_id) + '", synonym:"' + str(sl_synonym) + '"} )-[:is_synonym_off]->(u);'
            creation_nodes_list.append(creat_node_connection)
        return creation_nodes_list

def extract(**kwargs):
        ti = kwargs['ti']
        extractor=SwissLipidsExtractor()
        raw_df=extractor.extract()
        ti.xcom_push('order_data', raw_df)

#@task  # add parameters to decorator eventually
def transform(**kwargs):
    ti = kwargs['ti']
    extract_df = ti.xcom_pull(task_ids='extract', key='order_data')
    transformer=SwissLipidsTransformer()
    treated_data=transformer.transform(extract_df)
    

#  # add parameters to decorator eventually
def load(**kwargs):
    """
    Method where the load method will be added.
    """
    loader = SwissLipidsLoader()
    loader.load()        
    
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