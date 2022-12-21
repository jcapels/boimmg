import shutil
import urllib
from datetime import datetime

from src.boimmgpy.utilities import file_utilities
from src.boimmgpy.definitions import ROOT_DIR, DATABASE_CONFIGS


class ReactionsIDConverterScraper:

    def __init__(self):
        self.url_db = "https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/Aliases" \
                      "/Unique_ModelSEED_Reaction_Aliases.txt "

        self.__database_config = file_utilities.read_conf_file(DATABASE_CONFIGS)
        up_to_date = self.__check_time()
        if not up_to_date:
            self.__download_modelseed_database()

    def __check_time(self):
        last_update_time = self.__database_config["time_model_seed_reactions_aliases"]
        last_update_time = datetime.strptime(last_update_time, "%Y-%m-%d")
        time_now = datetime.now()
        delta = time_now - last_update_time
        if delta.days >= 7:
            return False
        else:
            return True

    def __download_modelseed_database(self):

        with urllib.request.urlopen(self.url_db) as response, open(
                ROOT_DIR + self.__database_config["path_model_seed_reactions_aliases"],
                'wb') as out_file:
            shutil.copyfileobj(response, out_file)

        self.__database_config["time_model_seed_reactions_aliases"] = str(datetime.now()).split(" ")[0]
        file_utilities.change_conf_file(DATABASE_CONFIGS, self.__database_config)

    def __call__(self):
        return ROOT_DIR + self.__database_config["path_model_seed_reactions_aliases"]
