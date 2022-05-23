import pandas as pd
import pendulum
from airflow.decorators import task
from rdkit.Chem import PandasTools
from airflow.models.dag import dag
from airflow.operators.python import PythonOperator
import requests, zipfile, io


from boimmgpy.etl.airflow_interfaces import AirflowExtractor, AirflowTransformer, AirflowLoader, AirflowPipeline


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


    def _extract(self, raw) -> pd.DataFrame:
        """
        Method to create a pandas dataframe with the scraped data.
        :return: pandas data frame of lipid maps data 
        """
        raw=PandasTools.LoadSDF(raw)
        df=pd.DataFrame(raw)
        #print(df.head())
        return df

    def scrape_data(self):
        """
        This class downloads the ZIP file and extracts the sdf file of lipid maps.
        :return: raw sdf file of lipid maps data
        """
        raw_file = requests.get("https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip")
        file_unziped = zipfile.ZipFile(io.BytesIO(raw_file.content))
        return file_unziped.open('structures.sdf')  



class LipidMapsTransformer(AirflowTransformer):

    def __init__(self) -> None:
        self.treated_dataframe=pd.DataFrame(columns=['LM_ID','ABBREVIATION','SYNONYMS'])

    def transform(self):
        """
        This method allows to transform the dataframe previously extracted into a desired datastructure
        :return:
        """
        data=LipidMapsExtractor()
        data_treated=self.treat(data.extract())
        return data_treated
    
    def three_col_lm(self,data_three):
        df_three_col=data_three[['LM_ID','ABBREVIATION','SYNONYMS']]
        return df_three_col

    def treat(self,data):
        df=self.three_col_lm(data)
        first_col=df[['LM_ID']].values
        second_col=df[['ABBREVIATION']].values
        third_col=df[['SYNONYMS']].values

        for n in range (len(first_col)): #len(first_col)
            abrevs=second_col[n]
            synonyms=third_col[n]
            lista_ab=[]
            lista_sy=[]
            lista_id=[]
            if abrevs != None:
                for value in abrevs:
                    value=str(value)
                    value=value.split(';')
                    for i in value:
                        lista_ab.append(i)

            if synonyms != None:    
                for value in synonyms:
                    value=str(value)
                    value=value.split(';')
                    for i in value:
                        lista_sy.append(i)

            for a in range(max([len(lista_ab),len(lista_sy)])):
                for value in first_col[n]:
                    value=str(value)
                    lista_id.append(value)

            # pegar em cada linha e adicionar ao df
            treated_data=pd.DataFrame.from_dict({'LM_ID': lista_id, 'ABBREVIATION': lista_ab, 'SYNONYMS':lista_sy},orient='index').T
            #df.append(treated_data,ignore_index=True,)
            self.treated_dataframe=pd.concat([self.treated_dataframe,treated_data],ignore_index=True,sort=False,axis=0)
        print(self.treated_dataframe.head())
        return self.treated_dataframe


class LipidMapsLoader(AirflowLoader):

    def load(self, **kwargs):
        """
        This method loads the data into the database
        :return:
        """
        pass


class LipidMapsETLPipeline(AirflowPipeline):

    #@task  # add parameters to decorator eventually
    def extract(self, **kwargs):
        """
        Method where the extraction method will be added.
        """
        extractor = LipidMapsExtractor()
        extractor.extract(**kwargs)

    #@task  # add parameters to decorator eventually
    def transform(self, **kwargs):
        """
        Method where the transform method will be added.
        """
        transformer = LipidMapsTransformer()
        transformer.transform(**kwargs)

    #@task  # add parameters to decorator eventually
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
