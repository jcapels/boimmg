from boimmgpy import definitions
from boimmgpy.utilities import file_utilities
import requests

class CompoundsRestAccessor:

    def __init__(self):

        conf = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)
        self.rest_uri = conf["rest_uri"]

    def get_node_from_model_seed_id(self, modelseed_id):

        res = requests.get(self.rest_uri + "rest/modelseed/node/" + modelseed_id)

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            #### TRASFORM json into node
            return results

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

        res = requests.get(self.rest_uri + "rest/conjugates" + str(lipid_id))

        if res.status_code == 200:
            res_dict = res.json()
            results = res_dict.get("result")
            return results

        return []

    def get_compounds_with_only_one_component(self, lipid_id,component):

        pass

    def get_compounds_with_specific_parent_within_set_of_components(self, lipid_id):

        pass

    def get_predecessors_by_ont_id_rel_type(self, lipid_id, relationship_type):

        pass

    def get_all_predecessors_by_ont_id_rel_type(self, lipid_id, relationship_type):
        pass

    def get_node_by_ont_id(self, lipid_id):

        pass

    def get_compounds_with_specific_parent_set_of_components(self, parent):

        pass

    def get_leaves_from_ont_id(self, lipid_id):

        pass

    def get_successors_by_ont_id_rel_type(self, lipid_id, relationship_type):

        pass