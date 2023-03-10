from .Lipid import Lipid
from .LipidMapsScraper import LipidMapsScraper
import pandas as pd


class LipidMapsStructureDB:

    def __init__(self):
        scraper = LipidMapsScraper()
        path = scraper()
        self.__read_database(path)

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

    def __read_database(self, path):
        self.__database = {}
        self.__inchikey_database = {}
        self.__smiles_database = {}
        data = pd.read_csv(path, header=0, delimiter="\t", encoding='ISO-8859-1')

        for i in range(data.shape[0]):

            inchikey = data.loc[:, "inchi_key"].iloc[i]
            smiles = data.loc[:, "smiles"].iloc[i]

            id = data.loc[:, "regno"].iloc[i]
            sys_name = data.loc[:, "sys_name"].iloc[i]
            name = data.loc[:, "name"].iloc[i]
            formula = data.loc[:, "formula"].iloc[i]
            inchi = data.loc[:, "inchi"].iloc[i]
            mass = data.loc[:, "exactmass"].iloc[i]
            kegg_id = data.loc[:, "kegg_id"].iloc[i]
            hmdb_id = data.loc[:, "hmdb_id"].iloc[i]
            chebi_id = data.loc[:, "chebi_id"].iloc[i]
            lipidbank_id = data.loc[:, "lipidbank_id"].iloc[i]
            pubchem_cid = data.loc[:, "pubchem_cid"].iloc[i]

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
