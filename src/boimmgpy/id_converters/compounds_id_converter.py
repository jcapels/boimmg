import re
import urllib

from rdkit import Chem

from src.boimmgpy.utilities import file_utilities
from src.boimmgpy.id_converters.id_converter import IDConverter
from src.boimmgpy.definitions import TOOL_CONFIG_PATH, ROOT_DIR, DATABASE_CONFIGS


class CompoundsIDConverter(IDConverter):

    def __init__(self):
        super().__init__()

        self.__database_config = file_utilities.read_conf_file(DATABASE_CONFIGS)
        path = ROOT_DIR + self.__database_config["path_model_seed_compounds_aliases"]
        self.__configs = file_utilities.read_conf_file(TOOL_CONFIG_PATH)
        self.__home_path__ = ROOT_DIR
        self.construct_compounds_converter(path)

    def construct_compounds_converter(self, path):

        info = open(path, "r")
        extra_info = open(self.__home_path__ + self.__configs.get("compound_converter_extra"), "r")
        url = False

        previousModelSeed = ""
        previousOtherDb = ""
        new_dic = {}
        for line in info:

            if url:
                line = line.decode("utf-8")

            l = line.split("\t")

            modelSeedId = l[0]
            externalId = l[1]
            source = l[2].strip()

            if re.search("\AModelSEED", line) is None:

                if previousModelSeed != modelSeedId:
                    self.modelSeedToDb[previousModelSeed] = new_dic
                    new_dic = {source: [externalId]}
                    previousModelSeed = modelSeedId
                    previousOtherDb = source

                elif previousOtherDb != source:
                    new_dic[source] = [externalId]
                    previousOtherDb = source

                else:
                    new_dic[source].append(externalId)

                if source in self.dbToModelSeed.keys():
                    if externalId in self.dbToModelSeed[source].keys():
                        self.dbToModelSeed[source][externalId].append(modelSeedId)
                    else:
                        self.dbToModelSeed[source][externalId] = [modelSeedId]

                else:
                    dic = {externalId: [modelSeedId]}
                    self.dbToModelSeed[source] = dic
        info.close()

        for line in extra_info:

            if url:
                line = line.decode("utf-8")

            l = line.split("\t")

            modelSeedId = l[0]
            externalId = l[1]
            source = l[2].strip()

            if re.search("\AModelSEED", line) is None:

                if modelSeedId in self.modelSeedToDb.keys():
                    if source in self.modelSeedToDb.get(modelSeedId) and \
                            externalId not in self.modelSeedToDb[modelSeedId][source]:
                        self.modelSeedToDb[modelSeedId][source].append(externalId)
                    else:
                        self.modelSeedToDb[modelSeedId][source] = [externalId]

                if source in self.dbToModelSeed.keys():
                    if externalId in self.dbToModelSeed[source].keys():
                        self.dbToModelSeed[source][externalId].append(modelSeedId)
                    else:
                        self.dbToModelSeed[source][externalId] = [modelSeedId]

                else:
                    dic = {externalId: [modelSeedId]}
                    self.dbToModelSeed[source] = dic

        extra_info.close()

    @staticmethod
    def convert_molfile_from_kegg_into_inchikey(compound_id):

        data = None
        try:
            data = urllib.request.urlopen("https://www.genome.jp/dbget-bin/www_bget?-f+m+compound+" + compound_id)
        except:
            print("error connecting")
        molformat = ""

        for line in data:
            molformat += line.decode('utf-8')

        try:
            mol = Chem.MolFromMolBlock(molformat)
            inchikey = Chem.MolToInchiKey(mol)
        except:
            return None
        return inchikey

    @staticmethod
    def find_in_model_seed_by_inchikey(inchikey, modelseed_database):
        compound_id = ""
        for key in modelseed_database:
            if key == inchikey:
                compound_id = modelseed_database[key].get("model_seed_id")
                break
        return compound_id
