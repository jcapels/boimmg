import requests
from boimmgpy import definitions
from boimmgpy.utilities import file_utilities


class RestUtils:

    @staticmethod
    def map_model(model,database):

        conf = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)
        rest_uri = conf["rest_uri"]

        with open(model, 'rb') as file:

            file = {'model': file}
            response = requests.post(rest_uri + "rest/map_model/" + database, files = file)

        return response