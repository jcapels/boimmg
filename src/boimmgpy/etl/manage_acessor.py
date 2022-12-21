from src.boimmgpy.utilities import file_utilities
import pathlib
import os

print(os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs/my_database_settings.conf"))
class ManageAcess:

    def __init__(self):

        uri, user, password = self.read_config_file()

        self.__uri = uri
        self.__user = user
        self.__password = password

    @staticmethod
    def read_config_file():

        if pathlib.Path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs/my_database_settings.conf")).exists():
            configs = file_utilities.read_conf_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs/my_database_settings.conf"))

            if "uri" in configs.keys() and "user" in configs.keys() and "password" in configs.keys():
                uri = configs["uri"]
                user = configs["user"]
                password = configs["password"]

                return uri, user, password

        else:

            raise Exception("Please insert the required database information using set_database_information function")


acess=ManageAcess
print(acess.read_config_file())
