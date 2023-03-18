
from boimmgpy.etl.model_seed.compounds.ModelSeedCompound import ModelSeedCompound
from boimmgpy.etl.model_seed.compounds.model_seed_db_scraper import ModelSeedCompoundsDBScraper,ModelSeedStructuresDBScraper


class ModelSeedCompoundsDB:

    def __init__(self) -> None:
        self.__inchikey_database = {}
        self.__smiles_database = {}
        self.__compounds = []
        self.read_model_seed_compounds()


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


    def read_model_seed_structures():
        extractor = ModelSeedStructuresDBScraper()
        structure_dataframe = extractor.extract()
    

    def read_model_seed_compounds(self):
        self.__compounds_database = {}
        extractor = ModelSeedCompoundsDBScraper()
        compound_dataframe = extractor.extract()
        for i,row in compound_dataframe.iterrows():
            name = row["name"]
            formula = row["formula"]
            charge = int(row["charge"])
            inchikey = row["inchikey"]
            smiles = row["smiles"]
            _id = row["id"]
            modelseed_compound = ModelSeedCompound(_id, name, charge, formula)
            modelseed_compound.setStructures(smiles, inchikey)
            self.__compounds_database[_id] = modelseed_compound
            self.__compounds.append(modelseed_compound)




ModelSeedCompoundsDB()