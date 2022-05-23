import pandas as pd
import pendulum
from airflow.decorators import task
from rdkit.Chem import PandasTools
from airflow.models.dag import dag
from airflow.operators.python import PythonOperator
import requests, zipfile, io, gzip
from boimmgpy.etl.airflow_interfaces import AirflowExtractor, AirflowTransformer, AirflowLoader, AirflowPipeline


class SwissLipidsExtractor(AirflowExtractor):
    """
    Class to extract information from lipid maps
    """

    def extract(self, **kwargs):
        """
        This class calls the scrape_data method and creates a pandas dataframe with the scraped data.
        :return: pandas data frame of lipid maps data 
        """
        return self._extract(self.scrape_data())


    
    def scrape_data(self):
        """
        This class downloads the ZIP file and extracts the CSV files of SwissLipids.
        :return:
        """
        r = requests.get("https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv")
        z = gzip.open(io.BytesIO(r.content),'rb')
        return z
        
        


    def _extract(self,raw) -> pd.DataFrame:
        """
        Method to create a pandas dataframe with the scraped data.
        :return:
        """
        table= pd.read_table(raw,engine='python',encoding='ISO-8859-1')
        df=pd.DataFrame(table)
        return df

    

class SwissLipidsTransformer(AirflowTransformer):
    
    def __init__(self) -> None:
        self.treated_dataframe=pd.DataFrame(columns=['Lipid ID','Abbreviation*','Synonyms*'])

    def transform(self):
        """
        This method allows to transform the dataframe previously extracted into a desired datastructure
        :return:
        """
        data=SwissLipidsExtractor()
        data_treated=self.treat(data.extract())
        return data_treated

    def three_col_lm(self):
        df=pd.read_csv(r'C:\Users\ampsi\OneDrive\Ambiente de Trabalho\projeto\Projeto\Pratica\sl_treated.csv')
        df=df[['Lipid ID','Abbreviation*','Synonyms*']]
        return df

    def treat(self,df):
        first_col=df[['Lipid ID']].values
        second_col=df[['Abbreviation*']].values
        third_col=df[['Synonyms*']].values

        for n in range (len(first_col)):
            data=pd.DataFrame(columns=['Lipid ID','Abbreviation*','Synonyms*'])
            abrevs=second_col[n]
            synonyms=third_col[n]
            id=[]
            lista_ab=[]
            lista_sy=[]
            lista_id=[]
            if abrevs != None:
                for value in abrevs:
                    value=str(value)
                    value=value.split('|')
                    for i in value:
                        lista_ab.append(''.join(i))

            if synonyms != None:    
                for value in synonyms:
                    value=str(value)
                    value=value.split('|')
                    for i in value:
                        lista_sy.append(''.join(i))

            for a in range(max([len(lista_ab),len(lista_sy)])):
                for value in first_col[n]:
                    value=str(value)
                    lista_id.append(''.join(value))
        
            treated_data=pd.DataFrame.from_dict({'Lipid ID': lista_id, 'Abbreviation*': lista_ab, 'Synonym*':lista_sy},orient='index').T
            #df.append(treated_data,ignore_index=True,)
            
            self.treated_dataframe=pd.concat([self.treated_dataframe,treated_data],ignore_index=True,sort=False)
            
        return self.treated_dataframe


class SwissLipidsLoader(AirflowLoader):

    def load(self, **kwargs):
        """
        This method loads the data into the database
        :return:
        """
        pass


class SwissLipidsETLPipeline(AirflowPipeline):

    #@task  # add parameters to decorator eventually
    def extract(self, **kwargs):
        """
        Method where the extraction method will be added.
        """
        extractor = SwissLipidsExtractor()
        extractor.extract(**kwargs)

    #@task  # add parameters to decorator eventually
    def transform(self, **kwargs):
        """
        Method where the transform method will be added.
        """
        transformer = SwissLipidsTransformer()
        transformer.transform(**kwargs)

    #  # add parameters to decorator eventually
    def load(self, **kwargs):
        """
        Method where the load method will be added.
        """
        loader = SwissLipidsLoader()
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