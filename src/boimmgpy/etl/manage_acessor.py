from boimmgpy.utilities import file_utilities
import pathlib
import os

class ManageAcess:
    """Class for managing access to a database.
    """
    def __init__(self):
        """
        Initializes the class.

        :param: None
        :return: None
        """
        uri, user, password = self.read_config_file()

        self.__uri = uri
        self.__user = user
        self.__password = password

    @staticmethod
    def read_config_file()->tuple(str):
        """
        Reads the configuration file and returns the database settings.

        :raises Exception: If the configuration file is missing or does not contain the required information.
        :return: A tuple containing the URI, user, and password for the database.
        :rtype: tuple[str]
        """
        if pathlib.Path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs/my_database_settings.conf")).exists():
            configs = file_utilities.read_conf_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs/my_database_settings.conf"))

            if "uri" in configs.keys() and "user" in configs.keys() and "password" in configs.keys():
                uri = configs["uri"]
                user = configs["user"]
                password = configs["password"]

                return uri, user, password

        else:

            raise Exception("Please insert the required database information using set_database_information function")


