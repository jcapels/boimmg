import pathlib
from boimmgpy.utilities import file_utilities

def read_config_file():

    if pathlib.Path("/code/boimmgpy/configs/my_database_settings.conf").exists():
        configs = file_utilities.read_conf_file("/code/boimmgpy/configs/my_database_settings.conf")

        if "uri" in configs.keys() and "user" in configs.keys() and "password" in configs.keys():
            uri = configs["uri"]
            user = configs["user"]
            password = configs["password"]

            return uri, user, password

    else:

        raise Exception("Please insert the required database information using set_database_information function")

