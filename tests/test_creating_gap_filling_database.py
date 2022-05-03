import cobra
from biocyc import biocyc, Reaction
from cobra import Metabolite
from cobrababel import bigg

from boimmgpy.utilities import model_utilities

from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from boimmgpy.utilities.annotation_utils import AnnotationUtils

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


def generate_reaction_for_all_databases(reaction, complete_reactants, complete_products):
    pass

def generate_bigg_metabolite(model_seed_id,compoundsIdConverter,compartment):
    if model_seed_id and "BiGG" in compoundsIdConverter.get_modelSeedIdToDb().get(model_seed_id) \
            or "BiGG1" in compoundsIdConverter.get_modelSeedIdToDb().get(model_seed_id):
        if "BiGG" in compoundsIdConverter.get_modelSeedIdToDb().get(model_seed_id):
            bigg_ids = compoundsIdConverter.convert_modelSeedId_into_other_dbID(model_seed_id, "BiGG")
        else:
            bigg_ids = compoundsIdConverter.convert_modelSeedId_into_other_dbID(model_seed_id, "BiGG1")

        bigg_metabolite = None
        bigg_id = None
        found = False
        i = 0
        while not found and i < len(bigg_ids):
            try:
                bigg_id = bigg_ids[i]
                bigg_metabolite = bigg.get_bigg_metabolite(bigg_id)
                found = True
            except:
                i += 1

        if bigg_metabolite:
            compounds = []
            model_metabolite = Metabolite(bigg_id + "_" + compartment)
            aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(model_seed_id)
            annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)
            model_metabolite.annotation = annotation
            if bigg_metabolite.get("inchikey"):
                model_metabolite.annotation["inchikey"] = bigg_metabolite.get("inchikey")

            model_metabolite.formula = bigg_metabolite.get("formulae")[0]
            model_metabolite.charge = bigg_metabolite.get("charges")[0]
            model_metabolite.name = bigg_metabolite.get("name")
            model_metabolite.compartment = compartment
            compounds.append(model_metabolite)

            return compounds

        else:
            return generate_modelseed_metabolite(model, model_seed_id, compoundsIdConverter, modelseedCompoundsDb)

    else:
        return generate_modelseed_metabolite(model, model_seed_id, compoundsIdConverter, modelseedCompoundsDb)

def generate_kegg_metabolite():
    pass

def generate_modelseed_metabolite(model_seed_id,compoundsIdConverter):
    modelseed_compound = modelseedCompoundsDb.get_compound_by_id(model_seed_id)

    aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(model_seed_id)
    annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)
    compounds = []

    if not compartment:

        for compartment in model.compartments:
            db_id = modelseed_compound.getDbId()
            model_metabolite = Metabolite(db_id + "_" + compartment)
            model_metabolite.annotation = annotation
            model_metabolite.annotation["inchikey"] = modelseed_compound.getInchikey()
            model_metabolite.formula = modelseed_compound.getFormula()
            model_metabolite.charge = modelseed_compound.getCharge()
            model_metabolite.name = modelseed_compound.getName()
            model_metabolite.compartment = compartment
            compounds.append(model_metabolite)

    else:
        db_id = modelseed_compound.getDbId()
        model_metabolite = Metabolite(db_id + "_" + compartment)
        model_metabolite.annotation = annotation
        model_metabolite.annotation["inchikey"] = modelseed_compound.getInchikey()
        model_metabolite.formula = modelseed_compound.getFormula()
        model_metabolite.charge = modelseed_compound.getCharge()
        model_metabolite.name = modelseed_compound.getName()
        model_metabolite.compartment = compartment
        compounds.append(model_metabolite)

    return compounds


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
