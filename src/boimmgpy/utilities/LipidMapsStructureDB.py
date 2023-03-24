from boimmgpy.utilities.Lipid import Lipid
import pandas as pd
from boimmgpy.etl.lipid_maps.lipid_maps_synonyms import LipidMapsExtractor

class LipidMapsStructureDB:
    def __init__(self) -> None:
        self.read_database()

    def getDatabase(self):
        return self.__database

    def getInchiKeyDatabase(self):
        return self.__inchikey_database

    def getSmilesDatabase(self):
        return self.__smiles_database

    def get_compound_by_smiles(self, smiles):
        if smiles in self.__inchikey_database.keys():
            return self.__inchikey_database[smiles]

        else:
            return None

    def get_compound_by_inchikey(self, inchikey):
        if inchikey[:-1] in self.__inchikey_database.keys():
            return self.__inchikey_database[inchikey[:-1]]

        else:
            return None

    def get_compound_by_id(self, id):
        if id in self.__database:
            return self.__database[id]

        else:
            return None
    
    def read_database(self):
        scraper = LipidMapsExtractor()
        lipid_maps_db = scraper.extract()
        self.__database = {}
        self.__inchikey_database = {}
        self.__smiles_database = {}
        for i,row in lipid_maps_db.iterrows():
            inchikey = row["INCHI_KEY"]
            smiles = row["SMILES"]
            lipid_id = row["LM_ID"]
            sys_name = row["SYSTEMATIC_NAME"]
            name = row["NAME"]
            formula = row["FORMULA"]
            inchi = row["INCHI"]
            mass = row["EXACT_MASS"]
            kegg_id = row["KEGG_ID"]
            hmdb_id = row["HMDB_ID"]
            chebi_id = row["CHEBI_ID"]
            lipidbank_id = row["LIPIDBANK_ID"]
            pubchem_cid = row["PUBCHEM_CID"]

            aliases = {}

            if not pd.isna(kegg_id):
                aliases["KEGG"] = [kegg_id]

            if not pd.isna(hmdb_id):
                aliases["HMDB"] = [hmdb_id]

            if not pd.isna(chebi_id):
                aliases["ChEBI"] = [str(chebi_id)]

            if not pd.isna(lipidbank_id):
                aliases["LipidBank"] = [lipidbank_id]

            if not pd.isna(pubchem_cid):
                aliases["PubChem"] = [str(pubchem_cid)]

            new_lipid = Lipid(id, name, sys_name, mass, formula, inchi, inchikey, smiles, aliases)

            self.__database[id] = new_lipid
            if not pd.isna(inchikey):
                self.__inchikey_database[inchikey[:-1]] = new_lipid

            if not pd.isna(smiles):
                self.__smiles_database[smiles] = new_lipid

