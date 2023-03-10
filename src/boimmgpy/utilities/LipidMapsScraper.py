import shutil
import urllib
from datetime import datetime

from boimmgpy.utilities import file_utilities
from boimmgpy.definitions import ROOT_DIR, DATABASE_CONFIGS


class LipidMapsScraper:

    def __init__(self):
        self.url = "https://www.lipidmaps.org/rest/compound/lm_id/LM/all/download"
        self.__database_config = file_utilities.read_conf_file(DATABASE_CONFIGS)
        up_to_date = self.__check_time()
        if not up_to_date:
            self.__download_lipids_database()

    def __check_time(self):
        last_update_time = self.__database_config["time_lipid_maps"]
        last_update_time = datetime.strptime(last_update_time, "%Y-%m-%d")
        time_now = datetime.now()
        delta = time_now - last_update_time
        if delta.days >= 7:
            return False
        else:
            return True

    def __download_lipids_database(self):
        with urllib.request.urlopen(self.url) as response, open(ROOT_DIR + self.__database_config["path_lipid_maps"],
                                                                'wb') as out_file:
            shutil.copyfileobj(response, out_file)

        self.__database_config["time_lipid_maps"] = str(datetime.now()).split(" ")[0]
        file_utilities.change_conf_file(DATABASE_CONFIGS, self.__database_config)

        new_lines = []
        i_line = 0
        with open(ROOT_DIR + self.__database_config["path_lipid_maps"], "r", encoding="ISO-8859-1") as file:

            line = file.readline()
            new_lines.append(line)
            lines = file.readlines()
            for line in lines:
                i_line += 1
                line_list = line.split("\t")
                line_list.pop(1)
                line = "\t".join(line_list)
                new_lines.append(line)

                while len(line_list) >= 23:
                    line_list.pop(20)

                line = "\t".join(line_list)
                new_lines.append(line)

        with open(ROOT_DIR + self.__database_config["path_lipid_maps"], "w", encoding="utf-8") as file:
            file.writelines(new_lines)

    def __call__(self):
        return ROOT_DIR + self.__database_config["path_lipid_maps"]
