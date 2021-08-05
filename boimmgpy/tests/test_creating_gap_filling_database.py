import cobra
from biocyc import biocyc, Reaction
from boimmgpy.utilities import model_utilities

from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB

compounds_converter = CompoundsIDConverter()
compounds_model_seed = ModelSeedCompoundsDB()
compounds_database = CompoundsDBAccessor()


def get_degradation_pathways():
    biocyc.set_organism("meta")

    pathway = biocyc.get("Fatty-Acid-and-Lipid-Degradation")

    instances = pathway.instances

    return instances


def transform_reaction_into_its_instance(reaction: Reaction):
    left = reaction.compounds_left
    right = reaction.compounds_right

    generic_reactants, generic_products = get_generic_nodes(left, right)
    new_reactions = []

    for generic_reactant in generic_reactants:

        generic_reactant_children = compounds_database.get_all_predecessors_by_ont_id_rel_type(generic_reactant, "is_a")
        for child in generic_reactant_children:

            precursors = compounds_database.get_all_predecessors_by_ont_id_rel_type(child, "precursor_of")
            successors = compounds_database.get_successors_by_ont_id_rel_type(child, "precursor_of")
            precursors_successors = precursors + successors

            complete_reactants = [compounds_database.get_node_by_ont_id(child)]
            complete_products = []

            for generic_product in generic_products:

                for precursor in precursors_successors:
                    precursor_parent = compounds_database.get_successors_by_ont_id_rel_type(precursor, "is_a")

                    if precursor_parent:

                        if precursor_parent[0] == generic_product:
                            precursor_node = compounds_database.get_node_by_ont_id(precursor)
                            complete_products.append(precursor_node)
                            # TODO: add metabolite to a model reaction (write function)

            generate_reaction_for_all_databases(reaction, complete_reactants, complete_products)


# def generate_reaction_for_all_databases(reaction, complete_reactants, complete_products):
#
#     if container.model_seed_id:
#         ms_container = self.__modelseedCompoundsDb.get_compound_by_id(container.model_seed_id)
#         self.generate_new_metabolite(ms_container)
#
#     else:
#         self.generate_new_boimmg_metabolite(container)


def get_generic_nodes(left, right):
    generic_reactants = []
    generic_products = []

    for compound in left:
        metacyc_id = compound.id

        modelseed_id = compounds_converter.convert_db_id_to_model_seed_by_db_id(metacyc_id)

        if modelseed_id:
            modelseed_id = modelseed_id[0]

            node_id = compounds_database.get_node_id_from_model_seed_id(modelseed_id)
            generic_reactants.append(node_id)

    for compound in right:
        metacyc_id = compound.id

        modelseed_id = compounds_converter.convert_db_id_to_model_seed_by_db_id(metacyc_id)

        if modelseed_id:
            modelseed_id = modelseed_id[0]

            node_id = compounds_database.get_node_id_from_model_seed_id(modelseed_id)
            generic_products.append(node_id)

    return generic_reactants, generic_products


get_degradation_pathways()
