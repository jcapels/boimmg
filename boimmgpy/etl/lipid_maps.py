import pandas as pd
from rdkit.Chem import PandasTools
import pendulum
from airflow.decorators import task
from airflow.models.dag import dag
from airflow.operators.python import PythonOperator
import requests, zipfile, io

from airflow_interfaces import AirflowExtractor, AirflowTransformer, AirflowLoader, AirflowPipeline


class LipidMapsExtractor(AirflowExtractor):
    """
    Class to extract information from lipid maps
    """

    def extract(self, **kwargs):
        """
        This class calls the scrape_data method and creates a pandas dataframe with the scraped data.
        :return: pandas data frame of lipid maps data 
        """
        return self._extract(self.scrape_data())


    def _extract(self, raw, **kwargs) -> pd.DataFrame:
        """
        Method to create a pandas dataframe with the scraped data.
        :return: pandas data frame of lipid maps data 
        """
        raw=PandasTools.LoadSDF(raw)
        df=pd.DataFrame(raw)
        
        return df

    def scrape_data(self):
        """
        This class downloads the ZIP file and extracts the sdf file of lipid maps.
        :return: raw sdf file of lipid maps data
        """
        r = requests.get("https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip")
        z = zipfile.ZipFile(io.BytesIO(r.content))
    
        return z.open('structures.sdf')  



class LipidMapsTransformer(AirflowTransformer):

    def transform(self, **kwargs):
        """
        This method allows to transform the dataframe previously extracted into a desired datastructure
        :return:
        """
        pass


class LipidMapsLoader(AirflowLoader):

    def load(self, **kwargs):
        """
        This method loads the data into the database
        :return:
        """
        pass


class LipidMapsETLPipeline(AirflowPipeline):

    @task  # add parameters to decorator eventually
    def extract(self, **kwargs):
        """
        Method where the extraction method will be added.
        """
        extractor = LipidMapsExtractor()
        extractor.extract(**kwargs)

    @task  # add parameters to decorator eventually
    def transform(self, **kwargs):
        """
        Method where the transform method will be added.
        """
        transformer = LipidMapsTransformer()
        transformer.transform(**kwargs)

    @task  # add parameters to decorator eventually
    def load(self, **kwargs):
        """
        Method where the load method will be added.
        """
        loader = LipidMapsLoader()
        loader.load(**kwargs)

    @dag(schedule_interval=None,
         start_date=pendulum.datetime(2022, 1, 1, tz="UTC"),
         catchup=False,
         tags=['example'], )
    def run(self):
        extract_task = PythonOperator(
            task_id='extract',
            python_callable=self.extract,
        )
        transform_task = PythonOperator(
            task_id='transform',
            python_callable=self.transform,
        )
        load_task = PythonOperator(
            task_id='load',
            python_callable=self.load,
        )
        extract_task >> transform_task >> load_task
