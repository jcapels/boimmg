import abc
import time
from typing import Union

import requests
from neo4j import GraphDatabase

from boimmgpy import definitions
from boimmgpy.model_seed.model_seed_compound import ModelSeedCompound
from boimmgpy.utilities import file_utilities


class ModelSeedCompoundsDB:

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'login') and
                callable(subclass.login) and
                hasattr(subclass, 'get_predecessors_by_ont_id') and
                callable(subclass.get_predecessors_by_ont_id) and
                hasattr(subclass, 'get_predecessors_by_ont_id_rel_type') and
                callable(subclass.get_predecessors_by_ont_id_rel_type) and
                hasattr(subclass, 'get_node_by_ont_id') and
                callable(subclass.get_node_by_ont_id)
                or
                NotImplemented)

    @abc.abstractmethod
    def get_compound_by_id(self, ont_id: Union[int, str]) -> ModelSeedCompound:
        """Get predecessors using as parameter the database identifier"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_compound_by_inchi_key(self, inchikey: str) -> list:
        """Get predecessors using as parameter the database identifier"""
        raise NotImplementedError


class ModelSeedCompoundsDBRest(ModelSeedCompoundsDB):

    def __init__(self):
        conf = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)
        self.rest_uri = conf["rest_uri_model_seed"]

        self.sleep_time = 0.5

    def deserialize_compound(self, compound_dict):
        return ModelSeedCompound(compound_dict)

    def get_compound_by_id(self, id):
        time.sleep(self.sleep_time)
        res = requests.get(self.rest_uri + "db_id/" + id)

        return self.deserialize_compound(res)

    def get_compound_by_inchi_key(self, inchikey):
        time.sleep(self.sleep_time)
        res = requests.get(self.rest_uri + "inchi_key/" + inchikey)

        return self.deserialize_compound(res)


class ModelSeedCompoundsDBRaw(ModelSeedCompoundsDB):

    def __init__(self, uri="", user="", password=""):

        if not uri or not user or not password:
            uri, user, password = self.read_config_file()

        self.__uri = uri
        self.__user = user
        self.__password = password

    def read_config_file(self):
        configs = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)

        uri = configs["uri"]
        user = configs["user"]
        password = configs["password"]

        return uri, user, password

    def login(self):

        driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__password), encrypted=False)
        self.tx = driver

    @property
    def tx(self):
        return self._tx

    @tx.setter
    def tx(self, tx):
        self._tx = tx

    def get_compound_by_id(self, id):

        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:ModelSeedCompound) "
                                 "WHERE c.model_seed_id = $db_id "
                                 "RETURN p",
                                 db_id=id)

            node = result.single()

            if node:
                node_properties = node.get("p")

                return ModelSeedCompound(node_properties)

        return None

    def get_compound_by_inchi_key(self, inchikey):

        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:ModelSeedCompound) "
                                 "WHERE c.inchikey contains $inchikey "
                                 "RETURN c "
                                 ,
                                 inchikey=inchikey[:-1])
            data = result.single()
            if data:
                node_properties = data.get("p")

                return ModelSeedCompound(node_properties)

        return None

    # def get_compound_by_canonical_smiles(self,smiles):
    #     if smiles in self.__smiles_database.keys():
    #         return self.__smiles_database[smiles]
    #     return None

    # def get_compounds(self):
    #     return self.__compounds

    # def find_compound_by_inchikey(self,inchikey):
    #     for compound in self.__compounds_database:
    #         compound_container = self.__compounds_database[compound]
    #         db_compound_inchikey = compound_container.getInchikey()
    #         if db_compound_inchikey:
    #             if db_compound_inchikey[:-2] == inchikey[:-2]:
    #                 return compound_container
    #     return None

    # def __set_inchi_key_and_smiles_database(self,compound):
    #
    #     inchikey = compound.getInchikey()
    #     smiles = compound.getSmiles()
    #
    #     if inchikey or smiles:
    #
    #         self.__inchikey_database[inchikey[:-1]] = compound
    #         mol = MolFromSmiles(smiles)
    #         if mol:
    #
    #             canonical_smiles_not_isomer = MolToSmiles(mol,isomericSmiles=True)
    #
    #             canonical_smiles = MolToSmiles(mol)
    #
    #             try:
    #                 neutralized, _ = chemo_utilities.NeutraliseCharges(canonical_smiles)
    #
    #                 neutralized_not_isomer, _ = chemo_utilities.NeutraliseCharges(canonical_smiles_not_isomer)
    #
    #                 if neutralized in self.__smiles_database:
    #
    #                     self.__smiles_database[neutralized].append(compound)
    #
    #                 else:
    #
    #                     self.__smiles_database[neutralized] = [compound]
    #
    #                 if neutralized_not_isomer in self.isomer_smiles_database:
    #
    #                     self.isomer_smiles_database[neutralized_not_isomer].append(compound)
    #
    #                 else:
    #
    #                     self.isomer_smiles_database[neutralized_not_isomer] = [compound]
    #
    #
    #             except:
    #                 pass
    #
    #
    #
    # def __read_model_seed_structures(self, path):
    #     with open(path,"r") as f:
    #         f.readline()
    #         compounds = f.readlines()
    #         previous = compounds[0].split("\t")[0]
    #         mapped_features = {}
    #         for compound in compounds:
    #             features = compound.split("\t")
    #             id = features[0]
    #             if previous != id:
    #                 smiles = mapped_features["smile"]
    #                 if "inchikey" in mapped_features.keys():
    #                     inchikey = mapped_features["inchikey"]
    #
    #                 if "inchi" in mapped_features.keys():
    #                     inchi = mapped_features["inchi"]
    #
    #                 self.__compounds_database[previous].setStructures(smiles,inchikey,inchi)
    #
    #                 self.__set_inchi_key_and_smiles_database(self.__compounds_database[previous])
    #
    #                 mapped_features = {}
    #             mapped_features[features[1].lower()] = features[5].strip()
    #             previous = id
    #
    #
    #
    # def __read_model_seed_compounds_database(self, path):
    #     self.__compounds_database = {}
    #     with open(path,"r") as f:
    #         f.readline()
    #         compounds = f.readlines()
    #         for compound in compounds:
    #
    #             features = compound.split("\t")
    #
    #             name = features[2]
    #             formula = features[3]
    #             charge = int(features[7])
    #             inchikey = features[6]
    #             smiles = features[19]
    #             modelseedCompound = ModelSeedCompound(features[0],name,charge,formula)
    #             modelseedCompound.setStructures(smiles,inchikey)
    #
    #
    #
    #             self.__compounds_database[features[0]] = modelseedCompound
    #             self.__compounds.append(modelseedCompound)
