from unittest import TestCase

from boimmgpy.database.accessors.compounds_database_accessor import set_database_information


class TestRestAccess(TestCase):

    def test_access(self):
        set_database_information(uri="bolt://palsson.di.uminho.pt:6094", user="neo4j", password="Tuning999")