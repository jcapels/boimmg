import urllib
import pandas as pd
from io import StringIO
import requests

class ModelSeedCompoundsDBScraper:
    """Class that extract Model SEED compounds database and creates a readable dataframe with it
    """

    def extract(self)->pd.DataFrame:
        """Method that extracts and transforms the Model SEED compounds database into a Pandas Dataframe by calling all static methods in the class

        Returns:
            pd.DataFrame: Pandas Dataframe with Model SEED compounds database
        """
        csv_file = self.scrape_data()
        data_frame = self.extract_data(csv_file)
        return data_frame
    
    @staticmethod
    def scrape_data():
        """Method that acesses and extracts Model SEED compounds database and turns it into a python readable format.
        Returns:
            _type_: Readable Model SEED compounds database
        """
        get_file = urllib.request.urlopen("https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.tsv")
        bytes_data = get_file.read()
        data = str(bytes_data,'utf-8')
        readable_data = StringIO(data) 
        return readable_data

    @staticmethod
    def extract_data (raw_file)->pd.DataFrame:
        """Method that receives raw Model SEED compouds database and transforms it into a Pandas Dataframe

        Args:
            raw_file (_type_): raw database

        Returns:
            pd.DataFrame: Pandas Dataframe with Model Seed compounds Database
        """
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe

class ModelSeedStructuresDBScraper:
    """Class that extract Model SEED structures database and creates a readable dataframe with it
    """
    def extract(self)->pd.DataFrame:
        """Method that extracts and transforms the Model SEED structures database into a Pandas Dataframe by calling all static methods in the class

        Returns:
            pd.DataFrame: Pandas Dataframe with Model SEED structures database
        """
        csv_file = self.scrape_data()
        data_frame = self.extract_data(csv_file)
        return data_frame
    
    @staticmethod
    def scrape_data():
        """Method that acesses and extracts Model SEED structures database and turns it into a python readable format.
        Returns:
            _type_: Readable Model SEED structures database
        """
        get_file = urllib.request.urlopen("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Structures/Unique_ModelSEED_Structures.txt")
        bytes_data = get_file.read()
        data = str(bytes_data,'utf-8')
        readable_data = StringIO(data) 
        return readable_data

    @staticmethod
    def extract_data (raw_file):
        """Method that receives raw Model SEED structures database and transforms it into a Pandas Dataframe

        Args:
            raw_file (_type_): raw database

        Returns:
            pd.DataFrame: Pandas Dataframe with Model Seed structures Database
        """
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe
    