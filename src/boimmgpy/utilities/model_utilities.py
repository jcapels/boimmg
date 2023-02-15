from cobra import Metabolite, Reaction, Model
from cobrababel import bigg

from boimmgpy.database.containers.compound_node import CompoundNode
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from boimmgpy.utilities.annotation_utils import AnnotationUtils
from boimmgpy.kegg.kegg_compound import KeggCompound
from boimmgpy.kegg.kegg_reaction import KeggReaction

import numpy as np


def check_if_elem_exist_in_list_of_cobra_containers(lst, elem):
    """
    Function to check if a given element is in a list with cobrapy containers (Metabolite, Reaction, etc)
    :param list lst: list of cobrapy containers
    :param elem: element
    :return boolean: if it exists in the list returns True, otherwise returns False
    """

    for el in lst:
        if elem.id == el.id:
            return True

    return False


def set_objective_function(model, reaction_id):
    reaction = model.reactions.get_by_id(reaction_id)
    model.objective = reaction


def evalSlimSol(solution, tol):
    if np.isnan(solution):
        return False

    if abs(solution) < tol:
        return False

    return True


def get_unique_list_of_cobra_containers(lst):
    """
    This function checks whether there are repetitive cobrapy containers in a given list and returns a list with unique
    containers.

    :param list lst: list of cobrapy containers
    :return list: list with unique cobrapy containers
    """

    res = []
    for elem in lst:
        if not check_if_elem_exist_in_list_of_cobra_containers(res, elem):
            res.append(elem)
    return res


#################################################### Generate model compounds ###############################################


def generate_model_compounds_by_database_format(model: Model,
                                                modelseedid: str,
                                                compoundsIdConverter: CompoundsIDConverter,
                                                modelseedCompoundsDb: ModelSeedCompoundsDB,
                                                database_format: str,
                                                compartment=""):
    """
    This function allows the generation of a given compound in the chosen database format.

    :param cobrapy.Model model: a cobrapy model
    :param string modelseedid: the model seed identifier
    :param CompoundsIDConverter compoundsIdConverter: a converter
    :param ModelSeedCompoundsDB modelseedCompoundsDb: a compound database
    :param string database_format: database name (only "ModelSEED", "BiGG" and "KEGG" supported)
    :param (Optional) string compartment: specific compartment to introduce compound
    :return list: a list of the generated compounds (in all the model compartments)
    """

    translation = compoundsIdConverter.get_modelSeedIdToDb().get(modelseedid)
    new_compounds = []
    if translation:
        if database_format == "BiGG" and \
                ("BiGG" in translation.keys() or "BiGG1" in translation.keys()):
            new_compounds = generate_bigg_compound(model,
                                                   modelseedid,
                                                   compoundsIdConverter,
                                                   modelseedCompoundsDb,
                                                   compartment)
            if new_compounds:
                return new_compounds

        elif database_format == "KEGG" and \
                "KEGG" in translation.keys():
            new_compounds = generate_kegg_compound(model,
                                                   modelseedid,
                                                   compoundsIdConverter,
                                                   modelseedCompoundsDb,
                                                   compartment)

            if new_compounds:
                return new_compounds

        if not new_compounds:
            new_compounds = generate_modelseed_compound(model,
                                                        modelseedid,
                                                        compoundsIdConverter,
                                                        modelseedCompoundsDb,
                                                        compartment)
            if new_compounds:
                return new_compounds

    else:
        new_compounds = generate_modelseed_compound(model,
                                                    modelseedid,
                                                    compoundsIdConverter,
                                                    modelseedCompoundsDb,
                                                    compartment)

        if new_compounds:
            return new_compounds


def generate_bigg_compound(model, model_seed_id, compoundsIdConverter, modelseedCompoundsDb, compartment=""):
    """
    Method to generate metabolites for each compartment of the model in the BiGG database format.
    This method tries to find the BiGG format of the new compounds. If it does not find, it will change it into the
    Model SEED format.

    :param cobrapy.Model model: cobrapy Model
    :param string model_seed_id: model seed identifier of the compound
    :param CompoundsIDConverter compoundsIdConverter: a converter
    :param ModelSeedCompoundsDB modelseedCompoundsDb: the compounds Model SEED database
    :param (Optional) string compartment: compartment to introduce the new compound
    :return list: new compounds list
    """

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
            if not compartment:
                for compartment in model.compartments:
                    model_metabolite = Metabolite(bigg_id + "_" + compartment)
                    aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(model_seed_id)
                    annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)
                    model_metabolite.annotation = annotation
                    if bigg_metabolite.get("inchikey"):
                        bigg_inchikey = bigg_metabolite.get("inchikey")
                        if bigg_inchikey:
                            model_metabolite.annotation["inchikey"] = bigg_inchikey

                    model_metabolite.formula = bigg_metabolite.get("formulae")[0]
                    model_metabolite.charge = bigg_metabolite.get("charges")[0]
                    model_metabolite.name = bigg_metabolite.get("name")
                    model_metabolite.compartment = compartment
                    compounds.append(model_metabolite)

            else:
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
            return generate_modelseed_compound(model, model_seed_id, compoundsIdConverter, modelseedCompoundsDb)

    else:
        return generate_modelseed_compound(model, model_seed_id, compoundsIdConverter, modelseedCompoundsDb)


def generate_modelseed_compound(model, model_seed_id, compoundsIdConverter, modelseedCompoundsDb, compartment=""):
    """
    Method to generate metabolites for each compartment of the model in the Model SEED database format.
    This method tries to find the Model SEED format of the new compounds.

    :param cobrapy.Model model: cobrapy Model
    :param string model_seed_id: model seed identifier of the compound
    :param CompoundsIDConverter compoundsIdConverter: a converter
    :param ModelSeedCompoundsDB modelseedCompoundsDb: the compounds Model SEED database
    :param (Optional) string compartment: compartment to introduce the new compound
    :return list: new compounds list
    """

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


def generate_boimmg_metabolites(model: Model, compound: CompoundNode,
                                database_format: str, compoundsIdConverter: CompoundsIDConverter,
                                compoundsAnnotationConfigs) -> list:
    """
    Generate an arbitrary format for the metabolites

    :param model:
    :param compound:
    :param database_format:
    :param compoundsIdConverter:
    :param compoundsAnnotationConfigs:
    :return:
    """

    compounds = []
    modelseedid = compound.model_seed_id

    id = compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + str(compound.id)

    formula = compound.formula
    name = compound.name
    charge = compound.charge
    if charge:
        charge = int(charge)

    model_compartments = model.compartments

    for compartment in model_compartments:
        new_metabolite = Metabolite(id=id + "_" + compartment, name=name, charge=charge, formula=formula)
        annotation = {}
        if modelseedid:
            aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(modelseedid)

            annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)

        annotation["boimmg.compound"] = compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + \
                                        str(compound.id)

        new_metabolite.annotation = annotation
        new_metabolite.annotation["smiles"] = compound.smiles
        new_metabolite.annotation["inchikey"] = compound.inchikey
        new_metabolite.compartment = compartment
        compounds.append(new_metabolite)

    return compounds


def generate_boimmg_metabolites_in_compartment(compound: CompoundNode,
                                               compartment: str, compoundsIdConverter: CompoundsIDConverter,
                                               compoundsAnnotationConfigs) -> Metabolite:
    """
    Generate an arbitrary format for the metabolites

    :param model:
    :param compound:
    :param database_format:
    :param compoundsIdConverter:
    :return:
    """

    modelseedid = compound.model_seed_id

    id = compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + str(compound.id)

    formula = compound.formula
    name = compound.name
    charge = compound.charge

    if charge:
        charge = int(charge)
    else:
        charge = 0

    new_metabolite = Metabolite(id=id + "_" + compartment, name=name, charge=charge, formula=formula)

    aliases = compound.aliases
    annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)
    if modelseedid:
        aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(modelseedid)

        annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)

    annotation["boimmg.compound"] = compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + \
                                    str(compound.id)

    new_metabolite.annotation = annotation
    new_metabolite.annotation["smiles"] = compound.smiles
    new_metabolite.annotation["inchikey"] = compound.inchikey
    new_metabolite.compartment = compartment

    return new_metabolite


def generate_kegg_compound(model, modelseedid, compoundsIdConverter, modelseedCompoundsDb, compartment=""):
    """
    Method to generate metabolites for each compartment of the model in the KEGG database format.
    This method tries to find the KEGG format of the new compounds. If it does not find, it will change it into the
    Model SEED format.

    :param cobrapy.Model model: cobrapy Model
    :param string modelseedid: model seed identifier of the compound
    :param CompoundsIDConverter compoundsIdConverter: a converter
    :param ModelSeedCompoundsDB modelseedCompoundsDb: the compounds Model SEED database
    :param (Optional) string compartment: compartment to introduce the new compound
    :return list: new compounds list
    """

    if modelseedid and "KEGG" in compoundsIdConverter.get_modelSeedIdToDb().get(modelseedid):
        kegg_id = compoundsIdConverter.convert_modelSeedId_into_other_dbID(modelseedid, "KEGG")[0]
        kegg_metabolite = KeggCompound(kegg_id)

        aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(modelseedid)
        annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)

        compounds = []
        if not compartment:
            for compartment in model.compartments:
                model_metabolite = Metabolite(kegg_id + "_" + compartment)
                model_metabolite.annotation = annotation
                model_metabolite.annotation["inchikey"] = kegg_metabolite.get_inchikey()
                model_metabolite.annotation["smiles"] = kegg_metabolite.get_smiles()
                model_metabolite.formula = kegg_metabolite.get_formula()
                model_metabolite.name = kegg_metabolite.get_name()
                model_metabolite.charge = 0
                model_metabolite.compartment = compartment
                compounds.append(model_metabolite)

        else:
            model_metabolite = Metabolite(kegg_id + "_" + compartment)
            model_metabolite.annotation = annotation
            model_metabolite.annotation["inchikey"] = kegg_metabolite.get_inchikey()
            model_metabolite.annotation["smiles"] = kegg_metabolite.get_smiles()
            model_metabolite.formula = kegg_metabolite.get_formula()
            model_metabolite.name = kegg_metabolite.get_name()
            model_metabolite.charge = 0
            model_metabolite.compartment = compartment
            compounds.append(model_metabolite)

        return compounds



    else:
        return generate_modelseed_compound(model, modelseedid, compoundsIdConverter, modelseedCompoundsDb)


#################################################### Generate model reactions ###############################################


def generate_reactions_to_model(model, reactions, babel, compartment,
                                compoundsIdConverter, modelseedCompoundsDb,
                                reactionsIdConverter, database_format, compoundsAnnotationConfigs):
    """
    This function will generate model reactions taking into account the :param reactions and the :param babel.
    The :param babel is a dictionary carrying the information about the compounds swapped in the model.
    It will generate the reactions that have the available reactants and products in a specific model compartment

    :param cobrapy.Model model: a cobrapy model
    :param list reactions: a list of reactions to possibly add to the model
    :param dict babel: a dictionary with the model seed id as a key and a list of model compounds.
    :param Model.Compartment compartment:
    :param CompoundsIDConverter compoundsIdConverter: a compounds ID converter
    :param ModelSeedCompoundsDB modelseedCompoundsDb: the Model SEED compounds database
    :param ReactionsIDConverter reactionsIdConverter: a reactions ID converter
    :param string database_format: database name (only "ModelSEED", "BiGG" and "KEGG" supported)
    :return:
    """

    reactions_to_add = []

    for modelseed_reaction in reactions:

        compounds = modelseed_reaction.getCompounds()
        reaction_found_compounds = babel.copy()
        for compound in compounds:
            aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(compound)
            modelseed_compound = modelseedCompoundsDb.get_compound_by_id(compound)
            inchikey = modelseed_compound.getInchikey()

            compounds_in_model = check_if_metabolite_exists_in_model(model, inchikey, aliases,
                                                                     compoundsAnnotationConfigs)
            if compounds_in_model:
                i = 0
                found = False
                compartment_compound = None
                while not found and i < len(compounds_in_model):
                    if compounds_in_model[i].compartment == compartment:
                        found = True
                        compartment_compound = compounds_in_model[i]
                    i += 1

                if compartment_compound:
                    reaction_found_compounds[compound] = compounds_in_model

        if len(reaction_found_compounds.keys()) == len(compounds):
            print(modelseed_reaction.getName())
            translation = reactionsIdConverter.get_modelSeedIdToDb().get(modelseed_reaction.getDbId())
            new_reaction = None
            if translation:
                if database_format == "BiGG" and \
                        ("BiGG" in translation.keys() or "BiGG1" in translation.keys()):
                    new_reaction = generate_bigg_reaction(modelseed_reaction, reaction_found_compounds,
                                                          compartment, reactionsIdConverter, compoundsIdConverter)
                    if new_reaction:
                        reactions_to_add.append(new_reaction)

                elif database_format == "KEGG" and \
                        "KEGG" in translation.keys():
                    new_reaction = generate_kegg_reaction(modelseed_reaction,
                                                          reaction_found_compounds,
                                                          compartment,
                                                          reactionsIdConverter, compoundsIdConverter)
                    if new_reaction:
                        reactions_to_add.append(new_reaction)

                if not new_reaction:
                    new_reaction = generate_modelseed_reaction(modelseed_reaction,
                                                               reaction_found_compounds,
                                                               compartment,
                                                               reactionsIdConverter)
                    if new_reaction:
                        reactions_to_add.append(new_reaction)

            else:
                new_reaction = generate_modelseed_reaction(modelseed_reaction,
                                                           reaction_found_compounds,
                                                           compartment,
                                                           reactionsIdConverter)
                if new_reaction:
                    reactions_to_add.append(new_reaction)

    model.add_reactions(reactions_to_add)


def generate_modelseed_reaction(modelseed_reaction, reaction_found_compounds, compartment, reactionsIdConverter):
    """
    This function generates a new Model SEED format reaction in a given compartment.

    :param ModelSeedReaction modelseed_reaction: a Model SEED reaction container
    :param dict reaction_found_compounds: dictionary with model seed compounds identifiers as keys and the
    respective Model.Metabolite as value.
    :param Model.Compartment compartment:
    :param ReactionsIDConverter reactionsIdConverter: a reaction identifier converter
    :return cobrapy.Reaction: the new model reaction
    """

    stoichiometry = modelseed_reaction.getStoichiometry()
    new_stoichiometry = {}

    for compound_id in stoichiometry:
        if compound_id in reaction_found_compounds.keys():
            model_compound = reaction_found_compounds[compound_id]
            new_stoichiometry[model_compound] = stoichiometry[compound_id]

    id = modelseed_reaction.getDbId() + "_" + compartment
    name = modelseed_reaction.getName()
    reversibility = modelseed_reaction.getReversibility()
    direction = modelseed_reaction.getDirection()
    if direction == "<":
        for compound in new_stoichiometry:
            new_stoichiometry[compound] = new_stoichiometry[compound] * -1

    if reversibility:
        new_model_reaction = Reaction(id, name, lower_bound=None)

    else:
        new_model_reaction = Reaction(id, name)

    aliases = reactionsIdConverter.get_all_aliases_by_modelSeedID(modelseed_reaction.getDbId())

    if aliases:
        new_model_reaction.annotation = AnnotationUtils.get_reaction_annotation_format_by_aliases(aliases)
    new_model_reaction.add_metabolites(new_stoichiometry)

    return new_model_reaction


def generate_bigg_reaction(modelseed_reaction, reaction_found_compounds,
                           compartment, reactionsIdConverter,
                           compoundsIdConverter):
    """
    This function generates a new BiGG format reaction in a given compartment.

    :param ModelSeedReaction modelseed_reaction: a Model SEED reaction container
    :param dict reaction_found_compounds: dictionary with model seed compounds identifiers as keys and the
    respective Model.Metabolite as value.
    :param Model.Compartment compartment:
    :param ReactionsIDConverter reactionsIdConverter: a reaction identifier converter
    :param CompoundsIDConverter compoundsIdConverter: a compounds identifier converter
    :return cobrapy.Reaction: the new model reaction
    """

    bigg_reaction_ids = \
        reactionsIdConverter.convert_modelSeedId_into_other_dbID(modelseed_reaction.getDbId(), "BiGG")

    bigg_reaction = None
    bigg_id = None
    found = False
    i = 0
    while not found and i < len(bigg_reaction_ids):
        try:
            bigg_reaction = bigg.get_bigg_reaction(bigg_reaction_ids[i])
            bigg_id = bigg_reaction_ids[i]
            found = True
        except:
            i += 1

    if bigg_reaction:
        bigg_metabolites = bigg_reaction.get("metabolites")

        new_stoichiometry = {}

        for compound in bigg_metabolites:
            compound_id = compound.get("bigg_id")
            modelseedid = compoundsIdConverter.convert_dbID_into_modelSeedId("BiGG", compound_id)[0]
            model_compound = reaction_found_compounds[modelseedid]
            new_stoichiometry[model_compound] = compound.get("stoichiometry")

        id = bigg_id + "_" + compartment
        name = bigg_reaction.get("name")
        reversibility = modelseed_reaction.getReversibility()

        if reversibility:
            new_reaction_model = Reaction(id, name, lower_bound=None)

        else:
            new_reaction_model = Reaction(id, name)

        new_reaction_model.add_metabolites(new_stoichiometry)
        aliases = reactionsIdConverter.get_all_aliases_by_modelSeedID(modelseed_reaction.getDbId())
        new_reaction_model.annotation = AnnotationUtils.get_reaction_annotation_format_by_aliases(aliases)
        return new_reaction_model


def generate_kegg_reaction(modelseed_reaction, reaction_found_compounds,
                           compartment, reactionsIdConverter, compoundsIdConverter):
    """
    This function generates a new KEGG format reaction in a given compartment.

    :param ModelSeedReaction modelseed_reaction: a Model SEED reaction container
    :param dict reaction_found_compounds: dictionary with model seed compounds identifiers as keys and the
    respective Model.Metabolite as value.
    :param Model.Compartment compartment:
    :param ReactionsIDConverter reactionsIdConverter: a reaction identifier converter
    :param CompoundsIDConverter compoundsIdConverter: a compounds identifier converter
    :return cobrapy.Reaction: the new model reaction
    """

    kegg_reaction_id = reactionsIdConverter.convert_modelSeedId_into_other_dbID(modelseed_reaction.getDbId(), "KEGG")[0]
    try:
        kegg_reaction = KeggReaction(kegg_reaction_id)

        stoichiometry = kegg_reaction.get_stoichiometry()
        new_stoichiometry = {}

        for compound_id in stoichiometry:
            modelseedids = compoundsIdConverter.convert_dbID_into_modelSeedId("KEGG", compound_id)

            found = False
            i = 0
            found_compound = None
            while not found and i < len(modelseedids):
                if modelseedids[i] in reaction_found_compounds.keys():
                    found = True
                    found_compound = modelseedids[i]

                i += 1

            if found:
                model_compound = reaction_found_compounds[found_compound]
                new_stoichiometry[model_compound] = stoichiometry[compound_id]

            else:
                return None

        id = kegg_reaction.get_id() + "_" + compartment

        name = kegg_reaction.get_name()
        if not name:
            name = modelseed_reaction.getName()

        reversibility = kegg_reaction.get_reversibility()
        direction = kegg_reaction.get_direction()
        if direction == "<=":
            for compound in new_stoichiometry:
                new_stoichiometry[compound] = new_stoichiometry[compound] * -1

        if reversibility:
            new_model_reaction = Reaction(id, name, lower_bound=None)

        else:
            new_model_reaction = Reaction(id, name)

        aliases = reactionsIdConverter.get_all_aliases_by_modelSeedID(modelseed_reaction.getDbId())
        new_model_reaction.annotation = AnnotationUtils.get_reaction_annotation_format_by_aliases(aliases)
        new_model_reaction.add_metabolites(new_stoichiometry)
        return new_model_reaction

    except:
        return None


def generate_reaction_in_compartimentalized_model(modelseed_reaction,
                                                  compartment, babel, reactionsIdConverter,
                                                  compoundsIdConverter, database_format):
    """
    This function will generate a model reaction taking into account the :param reactions and the :param babel.
    The :param babel is a dictionary carrying the information about the compounds swapped in the model.
    It will generate one reaction that have the available reactants and products in a specific model compartment.

    :param ModelSeedReaction modelseed_reaction:
    :param Model.Compartment compartment:
    :param dict babel: dictionary with specific information about some metabolites in the model
    :param CompoundsIDConverter compoundsIdConverter: a compounds ID converter
    :param ReactionsIDConverter reactionsIdConverter: a reactions ID converter
    :param string database_format: database name (only "ModelSEED", "BiGG" and "KEGG" supported)
    :return cobrapy.Reaction: new reaction

    """

    translation = reactionsIdConverter.get_modelSeedIdToDb().get(modelseed_reaction.getDbId())
    new_reaction = None
    if translation:
        if database_format == "BiGG" and \
                ("BiGG" in translation.keys() or "BiGG1" in translation.keys()):

            new_reaction = generate_bigg_reaction(modelseed_reaction,
                                                  babel,
                                                  compartment,
                                                  reactionsIdConverter,
                                                  compoundsIdConverter)
            if new_reaction:
                return new_reaction

        elif database_format == "KEGG" and \
                "KEGG" in translation.keys():
            new_reaction = generate_kegg_reaction(modelseed_reaction, babel, compartment,
                                                  reactionsIdConverter,
                                                  compoundsIdConverter)
            if new_reaction:
                return new_reaction

        if not new_reaction:
            new_reaction = generate_modelseed_reaction(modelseed_reaction,
                                                       babel, compartment,
                                                       reactionsIdConverter)
            if new_reaction:
                return new_reaction

    else:
        new_reaction = generate_modelseed_reaction(modelseed_reaction,
                                                   babel,
                                                   compartment,
                                                   reactionsIdConverter)
        if new_reaction:
            return new_reaction


########################################################## Checkers ####################################################

def check_if_reaction_exists_in_model(model, modelseed_reaction, reactionsIdConverter, reactionsAnnotationConfigs):
    """
    This function checks whether a given reaction is present in a given model or not.
    It returns the model reaction or None if not found
    :param model:
    :param modelseed_reaction:
    :param reactionsIdConverter:
    :param reactionsAnnotationConfigs:
    :return:
    """

    aliases = reactionsIdConverter.get_all_aliases_by_modelSeedID(modelseed_reaction.getDbId())

    for reaction in model.reactions:
        annotation = reaction.annotation

        for database in aliases:

            if database in reactionsAnnotationConfigs.keys():
                annotation_db = reactionsAnnotationConfigs[database]
            else:
                annotation_db = database

            if annotation_db in annotation.keys():
                if type(annotation[annotation_db]) == list:
                    for alias in annotation[annotation_db]:
                        if alias in aliases[database]:
                            return reaction
                elif annotation[annotation_db] in aliases[database]:
                    return reaction
    return None


def check_if_metabolite_exists_in_model(model, inchikey, aliases, compoundsAnnotationConfigs, boimmg_container=None):
    """
    This method searches for a given metabolite in the Model

    :param string inchikey: InChIKey of the metabolite to be found in the model
    :param dict aliases: dictionary with the metabolite aliases ({database : list with aliases})

    :return Model.metabolite: if the metabolite is found returns the metabolite, otherwise returns None
    """

    metabolites = []
    for metabolite in model.metabolites:
        found = False

        if inchikey:
            if "inchikey" in metabolite.annotation.keys():
                inchi_key = metabolite.annotation.get("inchikey")
                if inchi_key and inchi_key != "":
                    if inchi_key[:-1] == inchikey[:-1]:
                        # found_inchi_key_metabolite = True
                        if metabolite not in metabolites:
                            metabolites.append(metabolite)

        if not found:
            for alias in aliases:

                if alias in compoundsAnnotationConfigs.keys():
                    new_alias = compoundsAnnotationConfigs[alias]

                    if new_alias in metabolite.annotation.keys():
                        aliases_in_model = metabolite.annotation.get(new_alias)

                        if type(aliases_in_model) == list:
                            for alias_in_model in aliases_in_model:
                                if alias_in_model in aliases[alias]:
                                    # found_alias = True
                                    if metabolite not in metabolites:
                                        metabolites.append(metabolite)


                        elif aliases_in_model in aliases[alias]:
                            if metabolite not in metabolites:
                                metabolites.append(metabolite)

        if not found and boimmg_container:
            new_alias = compoundsAnnotationConfigs["BOIMMG"]

            if new_alias in metabolite.annotation.keys():
                aliases_in_model = metabolite.annotation.get(new_alias)
                boimmg_id = compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + str(boimmg_container.id)
                if aliases_in_model == boimmg_id:
                    metabolites.append(metabolite)

    return metabolites


############################################################################ Other ################################################################################

def get_compartments_babel(babel):
    """
    This function trawl through the babel in order to organize them by compartments.

    :param dict babel: a dictionary with the model seed id as a key and a list of model compounds.
    :return dict: dictionary of a dictionary with the compartments as keys and the babel for the compounds
    """
    res = {}
    modelseedids = list(babel.keys())
    for compound_in_model in babel[modelseedids[0]]:
        compartment = compound_in_model.compartment

        compartment_dict = {}
        compartment_dict[modelseedids[0]] = compound_in_model
        for i in range(1, len(modelseedids)):

            found = False
            j = 0
            compounds_in_model2 = babel[modelseedids[i]]
            while not found and j < len(compounds_in_model2):
                compound_compartment = compounds_in_model2[j].compartment
                if compound_compartment is None:
                    compound_compartment = "c"
                if compound_compartment == compartment:
                    compartment_dict[modelseedids[i]] = compounds_in_model2[j]
                    found = True
                j += 1

        if len(compartment_dict.keys()) == len(modelseedids):
            res[compartment] = compartment_dict
    return res


def check_if_is_exchange_reaction(reaction):
    """
    This function checks whether a given reaction is an exchange reaction.

    :param reaction:
    :return boolean:
    """
    reactants = False
    products = False
    stoich = reaction.metabolites
    for metabolite in stoich:
        if stoich[metabolite] < 0:
            reactants = True
        elif stoich[metabolite] > 0:
            products = True
    if reactants and products:
        return False
    else:
        return True


def swap_only_metabolites_in_reaction(metabolite, other, reaction):
    """
    This method swap a given metabolite into another in a specific reaction.

    :param Model.Metabolite metabolite: metabolite to be changed
    :param Model.Metabolite other: metabolite that will replace the other
    :param Model.Reaction reaction: reaction where the metabolites will be swapped.
    :return:
    """

    met_to_add = {}
    met_to_subtract = {}
    stoichiometry = reaction.metabolites
    coef = stoichiometry.get(metabolite)
    met_to_add[other] = coef
    met_to_subtract[metabolite] = coef
    reaction.subtract_metabolites(met_to_subtract)
    reaction.add_metabolites(met_to_add)


def get_model_metabolites_by_name(model, names):
    res = []
    for name in names:
        res.extend(model.metabolites.query(name, attribute='name'))
    return res


def get_model_seed_id_by_model_compound(compound, configs, converter) -> list:
    annotation = compound.annotation
    model_seed_annotation_key = configs.get("ModelSEED")
    if model_seed_annotation_key in annotation.keys():

        model_seed_id = annotation.get(model_seed_annotation_key)

        if type(model_seed_id) == str:
            return model_seed_id

        else:
            return model_seed_id[0]

    else:
        annotation_keys = list(configs.keys())
        for key in annotation_keys:
            annotation_value = configs.get(key)
            if annotation_value in annotation.keys():
                alias = annotation.get(annotation_value)

                if type(alias) == str:
                    model_seed_id = converter.convert_db_id_to_model_seed_by_db_id(alias)
                    if model_seed_id:
                        return model_seed_id

                else:
                    i = 0
                    while i < len(alias):
                        db_id = alias[i]
                        model_seed_id = converter.convert_db_id_to_model_seed_by_db_id(db_id)
                        if model_seed_id:
                            return model_seed_id
                        i += 1
    return None


def convert_model_seed_reactions_into_model_reactions(model, reactions,
                                                      reactionsIdConverter, reactionsAnnotationConfigs,
                                                      compoundsIdConverter, modelseedCompoundsDb,
                                                      compoundsAnnotationConfigs, database_format,
                                                      model_container=None, modelseedid=None):
    filtered_reactions = []
    for reaction in reactions:
        inModel = check_if_reaction_exists_in_model(model,
                                                    reaction,
                                                    reactionsIdConverter,
                                                    reactionsAnnotationConfigs)

        go = False
        if not inModel:
            stoichiometry = reaction.getStoichiometry()

            if model_container and modelseedid:

                if stoichiometry[modelseedid] > 0:
                    go = True

            else:
                go = True

            if go:
                compounds = reaction.getCompounds()
                reaction_found_compounds = {}

                for compound in compounds:
                    aliases = compoundsIdConverter.get_all_aliases_by_modelSeedID(compound)
                    modelseed_compound = modelseedCompoundsDb.get_compound_by_id(compound)
                    inchikey = modelseed_compound.getInchikey()

                    compound_in_model = check_if_metabolite_exists_in_model(model,
                                                                            inchikey,
                                                                            aliases,
                                                                            compoundsAnnotationConfigs)
                    if compound_in_model:
                        reaction_found_compounds[compound] = compound_in_model

                if len(reaction_found_compounds.keys()) == len(compounds):
                    print(reaction.getName(), reaction.getDbId())
                    compartment_babel = get_compartments_babel(reaction_found_compounds)
                    for compartment in compartment_babel:
                        new_reaction = generate_reaction_in_compartimentalized_model(reaction,
                                                                                     compartment,
                                                                                     compartment_babel[
                                                                                         compartment],
                                                                                     reactionsIdConverter,
                                                                                     compoundsIdConverter,
                                                                                     database_format)

                        if new_reaction and not check_if_is_exchange_reaction(new_reaction):
                            # try:
                            #     model.reactions.get_by_id(new_reaction.id)
                            # except:

                            filtered_reactions.append(new_reaction)

    return filtered_reactions
