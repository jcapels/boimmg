from unittest import TestCase

from src.boimmgpy.database.accessors.compounds_database_accessor import set_database_information, CompoundsDBAccessor


class TestRestAccess(TestCase):

    def test_access(self):
        set_database_information(uri="bolt://localhost:7687", user="neo4j", password="Tuning999")
        print(CompoundsDBAccessor().get_predecessors_by_ont_id(0))