from numpy import diag_indices
import pandas as pd
from datetime import datetime
from airflow.decorators import task
from rdkit.Chem import PandasTools
from airflow.models.dag import dag
from airflow.operators.python import PythonOperator
import requests, zipfile, io
from neo4j import GraphDatabase
import sys
sys.path.insert(1,'.')
from airflow import DAG

from boimmgpy.etl.airflow_interfaces import AirflowExtractor, AirflowTransformer, AirflowLoader, AirflowPipeline


class LipidMapsExtractor(AirflowExtractor):
    """
    Class to extract information from lipid maps
    """

    def extract(self) -> pd.DataFrame:
        """This class calls the scrape_data method and creates a pandas dataframe with the scraped data returned from scrape_data method.

        Returns:
            pd.DataFrame: Pandas dataframe of Lipid Maps data
        """
        return self._extract(self.scrape_data())


    def _extract(self, raw ) -> pd.DataFrame:
        """
        Method to create a pandas dataframe with the scraped data.

        Args:
            raw (sdf file): sdf file of the whole lipid maps database 

        Returns:
            pd.DataFrame:  pandas data frame of lipid maps data
        """
        raw_lipid_maps_data=PandasTools.LoadSDF(raw)
        lm_dataframe=pd.DataFrame(raw_lipid_maps_data)
        return lm_dataframe

    def scrape_data(self):
        """This method downloads the ZIP file and extracts the sdf file of lipid maps.

        Returns:
            sdf file: Sdf file of lipid maps data
        """
        raw_file = requests.get("https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip")
        file_unziped = zipfile.ZipFile(io.BytesIO(raw_file.content))
        return file_unziped.open('structures.sdf')  



class LipidMapsTransformer(AirflowTransformer):
    """Class to transform the lipid maps dataframe
    """

    def transform(self, df : pd.DataFrame)->pd.DataFrame:
        """
        This method allows to transform the dataframe previously extracted into a desired data structure
        Args:
            df (pd.dataframe): raw_dataframe of lipid maps data
        :return: treated dataframe of lipid maps, only with two columns, synonym and ID
        """
        data_treated=self.treat_lm_dataframe(df)
        return data_treated
    
    def three_col_sl(self, data: pd.DataFrame) -> pd.DataFrame:
        """ This method uses the whole lipid maps dataframe to create one with only three columns desirable

        Args:
            data (pd.DataFrame): 

        Returns:
            pd.DataFrame: _description_
        """
        df_three_col=data[['LM_ID','ABBREVIATION','SYNONYMS']]
        return df_three_col
    
    def treat_lm_dataframe(self,lm_dataframe:pd.DataFrame)->pd.DataFrame:
        """ This method uses the whole lipid maps dataframe and creates another with only two desirable columns, ID and Synonym, 
        the last one englobes the abbreviation two

        Args:
            lm_dataframe (pd.DataFrame): Whole pandas dataframe, previosly extracted in the extract class

        Returns:
            pd.DataFrame: pandas dataframe with two columns of the initial data, ID and Synonym
        """
        new_df=pd.DataFrame(columns=['LM_ID','SYNONYMS'])
        counter=0
        for i, row in lm_dataframe.iterrows():
            lipid_id = row["LM_ID"]
            abreviation = row["ABBREVIATION"]
            synonyms = row["SYNONYMS"]
            if abreviation is not None and not pd.isnull(abreviation):
                abbreviation_splits = abreviation.split(';')
                for split in abbreviation_splits:
                    new_df.at[counter, "LM_ID"] = lipid_id
                    new_df.at[counter, "SYNONYMS"] = split
                    counter+=1

            if synonyms is not None and not pd.isnull(synonyms):
                synonyms_splits = synonyms.split(';')
                for split in synonyms_splits:
                    new_df.at[counter, "LM_ID"] = lipid_id
                    new_df.at[counter, "SYNONYMS"] = split
                    counter+=1
        return new_df


class LipidMapsLoader(AirflowLoader):
    """
    Class that loads the treated data into the database
    """

    def load(self, df: pd.DataFrame):
        """This method connects and loads the treated data into the database

        Args:
            df (pd.DataFrame): Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
        """
        self.set_synonym(self.get_connection_list(df))

    def set_synonym(self, connections:list):
        """ Method that link to the graph database and uploads all the data previously treated 

        Args:
            connections (list): List of the querys to upload the whole treated dataframe to the database
        """
        data_base_connection = GraphDatabase.driver(uri="bolt://localhost:7687",auth=("neo4j","potassio19"))
        session = data_base_connection.session()
        for i in connections:
            session.run(i)

    def get_connection_list(self,df : pd.DataFrame)->list:
        """ This method creates the querys necessary to upload the treated data into the database 

        Args:
            df (pd.DataFrame): Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load

        Returns:
            _type_: List of querys necessary to the upload of the whole dataframe
        """
        creation_nodes_list=[]
        for i,row in df.iterrows():
            lipid_maps_id=row["LM_ID"]
            lm_synonym=row["SYNONYMS"]
            creat_node_connection='MATCH (u:LipidMapsCompound)WHERE u.lipidmaps_id="' + str(lipid_maps_id) + '" merge (s: Synonym {lm_id : "' + str(lipid_maps_id) + '", synonym:"' + str(lm_synonym) + '"} )-[:is_synonym_off]->(u);'
            creation_nodes_list.append(creat_node_connection)
        return creation_nodes_list

dag=DAG(dag_id="dag_etl_lm",schedule_interval="@once",
        start_date=datetime(2022, 1, 1),
        catchup=False,
        )

#@task  # add parameters to decorator eventually
def extract(**kwargs):
        ti = kwargs['ti']
        extractor=LipidMapsExtractor()
        raw_df=extractor.extract()
        ti.xcom_push('order_data', raw_df)

#@task  # add parameters to decorator eventually
def transform(**kwargs):
    ti = kwargs['ti']
    extract_df = ti.xcom_pull(task_ids='extract', key='order_data')
    transformer=LipidMapsTransformer()
    treated_data=transformer.transform(extract_df)
    

#  # add parameters to decorator eventually
def load(**kwargs):
    """
    Method where the load method will be added.
    """
    loader = LipidMapsLoader()
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
