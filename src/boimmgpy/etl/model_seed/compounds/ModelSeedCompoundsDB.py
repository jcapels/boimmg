from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles

from boimmgpy.model_seed.model_seed_compound import ModelSeedCompound
from boimmgpy.model_seed.model_seed_compounds_db_scraper import ModelSeedCompoundsDBScraper
from boimmgpy.utilities import chemo_utilities


class ModelSeedCompoundsDB:

    def __init__(self):
        self.__inchikey_database = {}
        self.__smiles_database = {}
        scraper = ModelSeedCompoundsDBScraper()
        (path_db, path_structures) = scraper()
        self.__compounds = []
        self.__read_model_seed_compounds_database(path_db)
        self.__read_model_seed_structures(path_structures)
        self.__set_inchi_key_and_smiles_database()

    def get_compound_by_id(self, id):
        return self.__compounds_database[id]

    def get_compound_by_inchi_key(self, inchikey):
        if inchikey[:-1] in self.__inchikey_database.keys():
            return self.__inchikey_database[inchikey[:-1]]
        return None

    def get_compound_by_canonical_smiles(self, smiles):
        if smiles in self.__smiles_database.keys():
            return self.__smiles_database[smiles]
        return None

    def get_compounds(self):
        return self.__compounds

    def find_compound_by_inchikey(self, inchikey):
        for compound in self.__compounds_database:
            compound_container = self.__compounds_database[compound]
            db_compound_inchikey = compound_container.getInchikey()
            if db_compound_inchikey:
                if db_compound_inchikey[:-2] == inchikey[:-2]:
                    return compound_container
        return None

    def getSmilesDatabase(self):
        return self.__smiles_database

    def __set_inchi_key_and_smiles_database(self):

        self.isomer_smiles_database = {}
        for compound in self.__compounds:
            inchikey = compound.getInchikey()
            smiles = compound.getSmiles()

            if inchikey or smiles:

                self.__inchikey_database[inchikey[:-1]] = compound
                mol = MolFromSmiles(smiles)
                if mol:

                    canonical_smiles_not_isomer = MolToSmiles(mol, isomericSmiles=False)

                    canonical_smiles = MolToSmiles(mol)

                    try:
                        neutralized, _ = chemo_utilities.neutralise_charges(canonical_smiles)

                        neutralized_not_isomer, _ = chemo_utilities.neutralise_charges(canonical_smiles_not_isomer)

                        if neutralized in self.__smiles_database:

                            self.__smiles_database[neutralized].append(compound)

                        else:

                            self.__smiles_database[neutralized] = [compound]

                        if neutralized_not_isomer in self.isomer_smiles_database:

                            self.isomer_smiles_database[neutralized_not_isomer].append(compound)

                        else:

                            self.isomer_smiles_database[neutralized_not_isomer] = [compound]


                    except:
                        pass

    def __read_model_seed_structures(self, path):
        with open(path, "r") as f:
            f.readline()
            compounds = f.readlines()
            previous = compounds[0].split("\t")[0]
            mapped_features = {}
            for compound in compounds:
                features = compound.split("\t")
                id = features[0]
                if previous != id:
                    smiles = mapped_features["smile"]
                    if "inchikey" in mapped_features.keys():
                        inchikey = mapped_features["inchikey"]

                    if "inchi" in mapped_features.keys():
                        inchi = mapped_features["inchi"]

                    self.__compounds_database[previous].setStructures(smiles, inchikey, inchi)
                    mapped_features = {}
                mapped_features[features[1].lower()] = features[5].strip()
                previous = id

    def __read_model_seed_compounds_database(self, path):
        self.__compounds_database = {}
        with open(path, "r") as f:
            f.readline()
            compounds = f.readlines()
            for compound in compounds:
                features = compound.split("\t")

                name = features[2]
                formula = features[3]
                charge = int(features[7])
                inchikey = features[6]
                smiles = features[19]
                modelseedCompound = ModelSeedCompound(features[0], name, charge, formula)
                modelseedCompound.setStructures(smiles, inchikey)
                self.__compounds_database[features[0]] = modelseedCompound
                self.__compounds.append(modelseedCompound)
