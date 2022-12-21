import pathlib
import re

from boimmgpy import definitions
from boimmgpy.utilities import file_utilities

from neo4j import GraphDatabase, basic_auth


class DatabaseAccessManager:
    _driver = None

    @property
    def driver(self):
        if self._driver is None:
            raise Exception("Driver not initialized, call establish_connection first")
        else:
            return self._driver

    def establish_connection(self, uri: str = None, user: str = None, password: str = None):
        if uri is None or user is None or password is None:
            uri, user, password = self._get_credentials_from_file()
        self._driver = GraphDatabase.driver(uri, auth=basic_auth(user, password))

    def connect(self):
        self.establish_connection()
        return self.driver

    @staticmethod
    def write_config_file(uri: str, user: str, password: str):
        uri_bool = False
        user_bool = False
        pass_bool = False
        lines = []
        if pathlib.Path(definitions.BOIMMG_DATABASE).exists():
            with open(definitions.BOIMMG_DATABASE, "r") as f:
                lines = f.readlines()

        with open(definitions.BOIMMG_DATABASE, "w") as f:

            if lines:
                for line in lines:
                    if re.match("^uri=", line):
                        uri_bool = True
                        f.write("uri=" + uri + "\n")
                    elif re.match("^user=", line):
                        user_bool = True
                        f.write("user=" + user + "\n")
                    elif re.match("^password=", line):
                        pass_bool = True
                        f.write("password=" + password + "\n")
                    else:
                        f.write(line)

            if not uri_bool:
                f.write("uri=" + uri + "\n")

            if not user_bool:
                f.write("user=" + user + "\n")

            if not pass_bool:
                f.write("password=" + password + "\n")

    def close(self):
        self._driver.close()

    def run_query(self, query, params=None):
        with self._driver.session() as session:
            result = session.run(query, params)
        return result

    def _get_credentials_from_file(self):
        if not pathlib.Path(definitions.BOIMMG_DATABASE).exists():
            raise Exception("Database credentials file not found")
        else:
            uri, user, password = self.read_database_credentials()
            return uri, user, password

    @staticmethod
    def read_database_credentials():

        if pathlib.Path(definitions.BOIMMG_DATABASE).exists():
            configs = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)

            if "uri" in configs.keys() and "user" in configs.keys() and "password" in configs.keys():
                uri = configs["uri"]
                user = configs["user"]
                password = configs["password"]

                return uri, user, password

        else:

            raise Exception("Please insert the required database information using set_database_information function")


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
