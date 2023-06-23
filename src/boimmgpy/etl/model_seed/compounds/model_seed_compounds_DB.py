import re
from typing import List
import pandas as pd
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles
from boimmgpy.etl.model_seed.compounds.ModelSeedCompound import ModelSeedCompound
from boimmgpy.etl.model_seed.compounds.model_seed_db_scraper import ModelSeedCompoundsDBScraper,ModelSeedStructuresDBScraper
from boimmgpy.utilities import chemo_utilities
from tqdm import tqdm

class ModelSeedCompoundsDB:

    def __init__(self) -> None:
        self.__inchikey_database = {}
        self.__smiles_database = {}
        self.__compounds = []
        self._set_compounds_db()
        

    def _set_compounds_db(self):
        """Implements all necessary methods to set the Model SEED database
        """
        self.read_model_seed_compounds()
        self.read_model_seed_structures()
        self.set_inchi_key_and_smiles_database()
    
    def get_compound_by_id(self, id:str)->ModelSeedCompound:
        """
        Searches for a compound by the given ID.

        :param id: The ID to be searched.
        :type id: str

        :return: The Model Seed compound object with the given ID.
        :rtype: ModelSeedCompound
        """
        
        return self.__compounds_database[id]

    def get_compound_by_inchi_key(self, inchikey:str)->ModelSeedCompound:
        """
        Searches for a Model Seed compound by the given InChIKey.

        :param inchikey: The InChIKey to be searched.
        :type inchikey: str

        :return: The Model Seed compound object with the given InChIKey, or None if not found.
        :rtype: ModelSeedCompound or None
        """
        if inchikey[:-1] in self.__inchikey_database.keys():
            return self.__inchikey_database[inchikey[:-1]]
        return None

    def get_compound_by_canonical_smiles(self, smiles:str)->ModelSeedCompound:
        """
        Searches for a Model Seed compound by the given canonical SMILES.

        :param smiles: The canonical SMILES to be searched.
        :type smiles: str

        :return: The Model Seed compound object with the given canonical SMILES, or None if not found.
        :rtype: ModelSeedCompound or None
        """
        if smiles in self.__smiles_database.keys():
            return self.__smiles_database[smiles]
        return None

    def get_compounds(self) -> List[ModelSeedCompound]:
        """
        Method to access all Model Seed compounds.

        :return: List containing all Model Seed compound objects.
        :rtype: List[ModelSeedCompound]
        """
        return self.__compounds


    def get_transformed_ID(self) -> pd.DataFrame:
        """
        Method to access Model Seed transformations to KEGG IDs.

        :return: DataFrame with KEGG IDs and their corresponding Model Seed IDs.
        :rtype: pd.DataFrame
        """
        return self.__kegg_to_modelseed

    def find_compound_by_inchikey(self, inchikey:str)->ModelSeedCompound:
        """
        Method to access a compound with a given InChIKey.

        :param inchikey: The InChIKey to be searched.
        :type inchikey: str

        :return: The Model Seed compound object with the given InChIKey, or None if not found.
        :rtype: ModelSeedCompound or None
        """

        for compound in self.__compounds_database:
            compound_container = self.__compounds_database[compound]
            db_compound_inchikey = compound_container.getInchikey()
            if db_compound_inchikey:
                if db_compound_inchikey[:-2] == inchikey[:-2]:
                    return compound_container
        return None


    def read_model_seed_structures(self):
        """Method to read and set model seed structures
        """
        extractor = ModelSeedStructuresDBScraper()
        structure_dataframe = extractor.extract()
        mapped_features = {}
        previous_id = structure_dataframe['ID'].iloc[0]
        for i,row in structure_dataframe.iterrows():
            _id = row["ID"]
            if previous_id != _id:
                smiles = mapped_features["smile"]
                if "inchikey" in mapped_features.keys():
                        inchikey = mapped_features["inchikey"]

                if "inchi" in mapped_features.keys():
                    inchi = mapped_features["inchi"]

                self.__compounds_database[previous_id].setStructures(smiles, inchikey, inchi)
                mapped_features = {}
            mapped_features[row["Type"].lower()] = row["Structure"].strip()
            previous_id = _id

    def read_model_seed_compounds(self):
        """Method to read and set Model Seed compounds database with usefull information such as name, formula, charge, inchikey, smiles and id.
        """
        kegg_dataframe = []
        self.__compounds_database = {}
        extractor = ModelSeedCompoundsDBScraper()
        compound_dataframe = extractor.extract()
        for i,row in tqdm(compound_dataframe.iterrows()):
            name = row["name"]
            formula = row["formula"]
            charge = int(row["charge"])
            inchikey = row["inchikey"]
            smiles = row["smiles"]
            _id = row["id"]
            aliases = row["aliases"]
            if not pd.isna(aliases):
                match = re.search(r"KEGG:\s*([^|]+)",aliases)
                if match:
                    kegg_id = match.group()
                    kegg_id=re.sub("KEGG: ","",kegg_id)
                    kegg_dataframe.append(self.set_kegg_modelseed_id(kegg_id,_id))
            modelseed_compound = ModelSeedCompound(_id, name, charge, formula)
            modelseed_compound.setStructures(smiles, inchikey)
            self.__compounds_database[_id] = modelseed_compound
            self.__compounds.append(modelseed_compound)
        self.__kegg_to_modelseed = pd.concat(kegg_dataframe)

    @staticmethod
    def set_kegg_modelseed_id (kegg_ids:str,modelseed_id:str)->pd.DataFrame:
        """
        Set KEGG ID to ModelSeed ID relationships.

        :param kegg_ids: The KEGG IDs separated by semicolons.
        :type kegg_ids: str
        :param modelseed_id: The ModelSeed ID.
        :type modelseed_id: str

        :return: A DataFrame containing the KEGG ID to ModelSeed ID mappings.
        :rtype: pd.DataFrame
        """
        kegg_ids = kegg_ids.split(";")
        kegg_to_modelseed = pd.DataFrame(columns=["kegg_id","modelseed_id"])
        counter = 0
        for kegg_id in kegg_ids:
            kegg_to_modelseed.at[counter,"kegg_id"] = kegg_id.strip()
            kegg_to_modelseed.at[counter,"modelseed_id"] = modelseed_id
            counter += 1
        return kegg_to_modelseed


    def set_inchi_key_and_smiles_database(self):
        """Method to set Model Seed compounds inchi key and smiles databases
        """
        self.isomer_smiles_database = {}
        for compound in self.__compounds:
            inchikey = compound.getInchikey()
            smiles = compound.getSmiles()
            if inchikey or smiles:
                    if not pd.isna(inchikey):
                        self.__inchikey_database[inchikey[:-1]] = compound
                    
                    mol = None
                    
                    if not pd.isna(smiles):
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


if __name__ == '__main__':
    teste = ModelSeedCompoundsDB()