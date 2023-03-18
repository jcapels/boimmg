import urllib
import pandas as pd
from io import StringIO
import requests

class ModelSeedCompoundsDBScraper:

    def extract(self):
        csv_file = self.scrape_data()
        data_frame = self.extract_data(csv_file)
        return data_frame
    
    @staticmethod
    def scrape_data():
        get_file = urllib.request.urlopen("https://github.com/ModelSEED/ModelSEEDDatabase/raw/master/Biochemistry/compounds.tsv")
        bytes_data = get_file.read()
        data = str(bytes_data,'utf-8')
        readable_data = StringIO(data) 
        return readable_data

    @staticmethod
    def extract_data (raw_file):
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe

class ModelSeedStructuresDBScraper:
    def extract(self):
        csv_file = self.scrape_data()
        data_frame = self.extract_data(csv_file)
        print(data_frame.head())
        return data_frame
    
    @staticmethod
    def scrape_data():
        get_file = urllib.request.urlopen("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Structures/Unique_ModelSEED_Structures.txt")
        bytes_data = get_file.read()
        data = str(bytes_data,'utf-8')
        readable_data = StringIO(data) 
        return readable_data

    @staticmethod
    def extract_data (raw_file):
        read_table = pd.read_table(raw_file, engine='python', encoding='ISO-8859-1')
        sl_dataframe = pd.DataFrame(read_table)
        return sl_dataframe
    