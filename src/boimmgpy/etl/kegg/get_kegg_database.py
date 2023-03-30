from io import StringIO
import re
import pandas as pd
import requests


class KeggDatabase:
    def get_kegg_database(self):
        self.scrape_kegg()
        return self.extract_kegg()

    @staticmethod
    def extract_kegg(self):
        kegg_database = pd.DataFrame(columns=["type","id","name"])
        txt_file = open("response.txt")
        lines = txt_file.readlines()
        counter = 0
        for line in lines:
            _type,line = self.catch_columns(line)
            if type != False:
                kegg_database.at[counter,"type"] = _type
                kegg_id,line = self.catch_columns(line)
                kegg_database.at[counter,"id"] = kegg_id
                name = line.strip()
                name = re.sub(r'\n',"", name)
                kegg_database.at[counter,"name"] = name
                counter+=1
            else:
                pass
        return kegg_database

    @staticmethod
    def catch_columns(line):
        match = re.search(r'^\s*\S+', line)
        line = re.sub(r'^\s*\S+',"",line)
        if match:
            word = match.group()
            word = word.strip()
            return word,line
        else:
            return False,line

    @staticmethod
    def scrape_kegg():
        raw_file = requests.get("https://www.genome.jp/kegg-bin/download_htext?htext=br08002.keg&format=htext&filedir=")
        with open("response.txt", "w") as f:
            f.write(raw_file.text)
        with open("response.txt", "r+") as f:
            lines = f.readlines()
            del lines[:4]
            lines = lines[:-5]  
            f.seek(0)
            f.truncate()
            f.writelines(lines)
        