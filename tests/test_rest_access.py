from unittest import TestCase

import pandas

from src.boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor, set_database_information
from src.boimmgpy.database.accessors.compounds_rest_accessor import CompoundsRestAccessor
from src.boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from src.boimmgpy.kegg.kegg_compound import KeggCompound


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

        ids = accessor.get_compounds_with_only_one_component(11, [293815])

        assert 312156 in ids

    def test_get_compounds_with_specific_parent_within_set_of_components(self):
        accessor = CompoundsRestAccessor()

        ids = accessor.get_compounds_with_specific_parent_within_set_of_components(11, [293816, 5497])
        assert 296497 in ids

    def test_files_change(self):

        set_database_information(uri="COISO",user="hello",password="password")

    def test_get_all_cofactors(self):

        accessor = CompoundsDBAccessor()
        id_converter = CompoundsIDConverter()

        ids = accessor.get_all_predecessors_by_ont_id_rel_type(749388, "is_a")
        names = []
        KEGG = []
        modelseed = []
        bigg = []

        for boimmg_id in ids:

            node = accessor.get_node_by_ont_id(boimmg_id)

            if node.model_seed_id:
                kegg_ids = id_converter.convert_modelSeedId_into_other_dbID(node.model_seed_id, "KEGG")
                if kegg_ids:
                    print(kegg_ids)
                    KEGG.append(kegg_ids[0])
                    names.append(node.name)
                    modelseed.append(node.model_seed_id)
                    bigg_ids = id_converter.convert_modelSeedId_into_other_dbID(node.model_seed_id, "BiGG")
                    bigg.append("\t".join(bigg_ids))

        ids = accessor.get_all_predecessors_by_ont_id_rel_type(749376, "is_a")

        for boimmg_id in ids:

            node = accessor.get_node_by_ont_id(boimmg_id)

            if node.model_seed_id:
                kegg_ids = id_converter.convert_modelSeedId_into_other_dbID(node.model_seed_id, "KEGG")
                if kegg_ids:
                    print(kegg_ids)
                    KEGG.append(kegg_ids[0])
                    names.append(node.name)
                    modelseed.append(node.model_seed_id)
                    bigg_ids = id_converter.convert_modelSeedId_into_other_dbID(node.model_seed_id, "BiGG")
                    bigg.append("\t".join(bigg_ids))

        df = pandas.read_csv("cofactors_merlin", sep="\t", header=0)
        for row, i in df.iterrows():
            kegg_id = row[2]
            modelseed_id = row[3]
            name = row[0]
            if kegg_id not in KEGG:
                bigg_ids = id_converter.convert_modelSeedId_into_other_dbID(modelseed_id, "BiGG")
                KEGG.append(kegg_id)
                names.append(name)
                modelseed.append(modelseed_id)
                bigg.append("\t".join(bigg_ids))

        base_kegg_id = "C000"
        for i in range(1, 95):
            kegg_id = base_kegg_id + str(i).zfill(2)
            try:
                compound = KeggCompound(kegg_id)

            except:
                compound = None

            if kegg_id not in KEGG and compound \
                    and "L-" not in compound.name and \
                    "D-" not in compound.name:
                modelseed_id = id_converter.convert_db_id_to_model_seed_by_db_id(kegg_id)

                bigg_ids = id_converter.convert_modelSeedId_into_other_dbID(modelseed_id[0], "BiGG")
                KEGG.append(kegg_id)
                bigg.append("\t".join(bigg_ids))
                modelseed.append(modelseed_id[0])
                names.append(compound.name)

        df = pandas.DataFrame({"name": names,
                               "kegg": KEGG,
                               "modelseed": modelseed,
                               "bigg": bigg})

        df.to_csv("cofactors.csv", index=False)
