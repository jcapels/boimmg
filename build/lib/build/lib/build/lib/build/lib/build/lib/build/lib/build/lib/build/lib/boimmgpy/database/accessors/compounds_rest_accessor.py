from boimmgpy import definitions
from boimmgpy.database.containers.compound_node import CompoundNode
from boimmgpy.database.interfaces.boimmg_database_accessor import BOIMMGDatabaseAccessor
from boimmgpy.utilities import file_utilities
import requests

class CompoundsRestAccessor(BOIMMGDatabaseAccessor):

    def __init__(self):

        conf = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)
        self.rest_uri = conf["rest_uri"]

    def deserialize_compound(self,compound_dict):

        return CompoundNode(compound_dict.get("id"),compound_dict,{})


    def get_node_from_model_seed_id(self, modelseed_id):

        res = requests.get(self.rest_uri + "rest/modelseed/node/" + modelseed_id)

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            node = self.deserialize_compound(results)
            return node

        return None

    def get_node_id_from_model_seed_id(self, modelseed_id):

        res = requests.get(self.rest_uri + "rest/modelseed/id/" + modelseed_id)

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return None

    def get_conjugates(self, lipid_id):

        res = requests.get(self.rest_uri + "rest/conjugates" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_predecessors_by_ont_id(self, lipid_id):

        res = requests.get(self.rest_uri + "rest/predecessors/" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_compounds_with_only_one_component(self, lipid_id,component):

        res = requests.post(self.rest_uri + "rest/compounds/same_components/" + str(lipid_id),
                            data={"components" : component})

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_compounds_with_specific_parent_within_set_of_components(self, lipid_id,components):

        res = requests.post(self.rest_uri + "rest/compounds/diff_components/" + str(lipid_id),
                            data={"components":str(components)})

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_predecessors_by_ont_id_rel_type(self, lipid_id, relationship_type):

        res = requests.get(self.rest_uri + "rest/predecessors/" + relationship_type + "/" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_all_predecessors_by_ont_id_rel_type(self, lipid_id, relationship_type):

        res = requests.get(self.rest_uri + "rest/predecessors/all/" + relationship_type + "/" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []


    def get_node_by_ont_id(self, lipid_id):

        res = requests.get(self.rest_uri + "rest/compounds/" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            node = self.deserialize_compound(results)
            return node

        return None

    def get_compounds_with_specific_parent_set_of_components(self, parent,components):

        res = requests.post(self.rest_uri + "rest/compounds/set_components/" + str(parent),
                            data={"components": components})

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_leaves_from_ont_id(self, lipid_id):

        res = requests.get(self.rest_uri + "rest/leaves/" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_successors_by_ont_id_rel_type(self, lipid_id, relationship_type):

        res = requests.get(self.rest_uri + "rest/successors/" + relationship_type + "/" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []