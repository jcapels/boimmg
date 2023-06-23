import urllib
import pandas as pd
from io import StringIO
import requests

class ModelSeedCompoundsDBScraper:
    """Extracts Model SEED compounds database and creates a readable dataframe with it.
    It acesses online Model SEED databas and performs a scrape to extract database files.
    """

    def extract(self)->pd.DataFrame:
        """Method that extracts and transforms the Model SEED compounds database into a Pandas Dataframe by calling Scraper and Extracter methods.

        :return: Original Model SEED dataframe arranged in a usefull Pandas Dataframe
        :rtype: pd.DataFrame
        """
        csv_file = self.scrape_data()
        data_frame = self.extract_data(csv_file)
        return data_frame
    
    @staticmethod
    def scrape_data()->StringIO:
        """
        Method that accesses and extracts the Model SEED compounds database and turns it into a Python-readable format.

        :return: A StringIO object containing the Model SEED compounds data.
        :rtype: StringIO
        """
        get_file = urllib.request.urlopen("https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.tsv")
        bytes_data = get_file.read()
        data = str(bytes_data,'utf-8')
        readable_data = StringIO(data) 
        return readable_data

    @staticmethod
    def extract_data (raw_file:StringIO)->pd.DataFrame:
        """
        Method that receives the raw Model SEED compounds database and transforms it into a Pandas DataFrame.

        :param raw_file: The raw database.
        :type raw_file: StringIO

        :return: Pandas DataFrame with the Model Seed compounds database.
        :rtype: pd.DataFrame
        """
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe

class ModelSeedStructuresDBScraper:
    """Class that extract Model SEED structures database and creates a readable dataframe with it
    """
    def extract(self)->pd.DataFrame:
        """Method that extracts and transforms the Model SEED structures database into a Pandas Dataframe by calling Scraper and Extracter methods.

        :return: Original Model SEED structures dataframe arranged in a usefull Pandas Dataframe
        :rtype: pd.DataFrame
        """
        csv_file = self.scrape_data()
        data_frame = self.extract_data(csv_file)
        return data_frame
    
    @staticmethod
    def scrape_data()->StringIO:
        """
        Method that accesses and extracts the Model SEED structures database and turns it into a Python-readable format.

        :return: A StringIO object containing the Model SEED structures data.
        :rtype: StringIO
        """
        get_file = urllib.request.urlopen("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Structures/Unique_ModelSEED_Structures.txt")
        bytes_data = get_file.read()
        data = str(bytes_data,'utf-8')
        readable_data = StringIO(data) 
        return readable_data

    @staticmethod
    def extract_data (raw_file:StringIO)->pd.DataFrame:
        """
        Method that receives the raw Model SEED structures database and transforms it into a Pandas DataFrame.

        :param raw_file: The raw database.
        :type raw_file: StringIO

        :return: Pandas DataFrame with the Model Seed structures database.
        :rtype: pd.DataFrame
        """
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe
    