from unittest import TestCase

from boimmgpy.database.accessors.compounds_rest_accessor import CompoundsRestAccessor


class TestRestAccess(TestCase):

    def test_get_node_from_model_seed_id(self):

        accessor = CompoundsRestAccessor()

        node = accessor.get_node_from_model_seed_id("cpd27342")

        assert node.model_seed_id == "cpd27342"
        assert node.id == 98
        assert node.formula == "C8H11NO10PR2"
        assert node.inchikey == None

    def test_get_node_id_from_model_seed_id(self):

        accessor = CompoundsRestAccessor()

        id = accessor.get_node_id_from_model_seed_id("cpd27342")

        assert id == 98

    def test_get_conjugates(self):

        accessor = CompoundsRestAccessor()

        id = accessor.get_conjugates(749395)

        assert id[0] == 749396

    def test_get_predecessors_by_ont_id(self):

        accessor = CompoundsRestAccessor()

        ids = accessor.get_predecessors_by_ont_id(305731)

        assert len(ids) == 2

    def test_get_compounds_with_only_one_component(self):


        accessor = CompoundsRestAccessor()

        ids = accessor.get_compounds_with_only_one_component(11,[293815])

        assert 312156 in ids

    def test_get_compounds_with_specific_parent_within_set_of_components(self):
        accessor = CompoundsRestAccessor()

        ids = accessor.get_compounds_with_specific_parent_within_set_of_components(11, [293816,5497])
        assert 296497 in ids

