import logging
from logging.handlers import TimedRotatingFileHandler
from typing import Union

from cobra import Model, Metabolite

from src.boimmgpy.database.containers.compound_node import CompoundNode
from src.boimmgpy.database.interfaces.boimmg_database_accessor import BOIMMGDatabaseAccessor
from src.boimmgpy.model_seed.model_seed_compound import ModelSeedCompound
from src.boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from src.boimmgpy.utilities import file_utilities, model_utilities
from src.boimmgpy.utilities.annotation_utils import AnnotationUtils
from src.boimmgpy.service.network_modifiers.reactions_swapper import ReactionsChanger
from src.boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from src.boimmgpy.kegg.kegg_compound import KeggCompound
from cobrababel import bigg
from src.boimmgpy.definitions import ROOT_DIR, TOOL_CONFIG_PATH, REACTIONS_ANNOTATION_CONFIGS_PATH, \
    COMPOUNDS_ANNOTATION_CONFIGS_PATH

logPath = ROOT_DIR + "/logs"

# format the log entries, backup limit and others
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
handler = TimedRotatingFileHandler(logPath, when='midnight', backupCount=20)
handler.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


def change_reaction_by_swapping_existing_metabolites(model_metabolite, new_metabolite):
    """
    Method to change a given metabolite and the respective reactions.

    :param Model.Metabolite model_metabolite: metabolite in the model to be replaced
    :param Model.Metabolite new_metabolite:  metabolite to replace.
    :return: list<Reaction> reactions_to_change: list with the reactions to be changed (format, id, annotation, etc...)
    """

    reactions_to_change = []
    if model_metabolite.id != new_metabolite.id:

        reactions = model_metabolite.reactions
        for reaction in reactions:
            met_to_add = {}
            met_to_subtract = {}
            stoichiometry = reaction.metabolites
            coef = stoichiometry.get(model_metabolite)
            met_to_add[new_metabolite] = coef
            met_to_subtract[model_metabolite] = coef
            reaction.subtract_metabolites(met_to_subtract)
            reaction.add_metabolites(met_to_add)
            reactions_to_change.append(reaction)

    return reactions_to_change


class MetaboliteSwapper:

    def __init__(self, model: Model, metabolite_in_model: Union[int, None], new_metabolite: Union[int, None],
                 type_change: int,
                 modelSeedCompoundsDb: ModelSeedCompoundsDB, db_accessor: BOIMMGDatabaseAccessor,
                 model_database: str = "ModelSEED",
                 compoundsIdConverter: Union[CompoundsIDConverter, None] = None,
                 universal_model: Union[Model, None] = None,
                 not_to_change_compounds: Union[list, None] = None):

        """
       Class constructor

       :param Model model: COBRApy model
       :param int metabolite_in_model: ontology id of the metabolite to be replaced
       :param int new_metabolite: ontology id of the metabolite that will replace
       :param int type_change: type of change (0, 1, 2 or 3)
       :param ModelSeedDB modelSeedCompoundsDb: Model SEED compounds database
       :param BOIMMGDatabaseAccessor: accessor (either raw or restAPI)
       :param string model_database: database format of metabolites and reactions
       :param (optional) CompoundsIDConverter compoundsIdConverter: compound id converter
       :param (optional) Model universal_model: an universal model
       :param (optional) list not_to_change_compounds: list with compound that are not supposed to be changed
       """

        if not_to_change_compounds is None:
            not_to_change_compounds = []
        if not universal_model:
            self.__universal_model = universal_model

        else:
            self.__universal_model = Model("universal_reactions")

        if not compoundsIdConverter:
            self.__compoundsIdConverter = CompoundsIDConverter()
        else:
            self.__compoundsIdConverter = compoundsIdConverter

        self.__compounds_ontology = db_accessor

        self.__model_database = model_database
        self.__compounds_db = modelSeedCompoundsDb
        self.__type = type_change
        self.__metabolite_in_model = metabolite_in_model
        self.__new_metabolite = new_metabolite
        self.__model = model
        self.__metabolites_not_to_change = not_to_change_compounds

        self.define_instance_variables()

    @property
    def model(self):
        return self.__model

    @model.setter
    def model(self, value):
        self.__model = value

    @property
    def mapper(self):
        return self.__mapper

    @mapper.setter
    def mapper(self, value):
        self.__mapper = value

    @property
    def changed_reactions(self):
        return self._changed_reactions

    @changed_reactions.setter
    def changed_reactions(self, value):
        self._changed_reactions = value

    def define_instance_variables(self):

        """
        Define instance variables

        :return:
        """

        self.report_material = {}
        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

        self.__home_path__ = ROOT_DIR
        self.__configs = file_utilities.read_conf_file(
            TOOL_CONFIG_PATH)

        self.changed_reactions = []
        self.reactions_swapper = ReactionsChanger(self.model,
                                                  self.__type, self.__model_database, self.__universal_model,
                                                  compounds_converter=self.__compoundsIdConverter)

    def set_type(self, new_type):
        """
        Method to change the type of swap.
        :param int new_type: new type of swap
        :return:
        """
        self.__type = new_type
        self.reactions_swapper.set_type(new_type)

    def set_new_compound(self, new_compound):
        """
        Method to set the new compound that will replace the other in the model.
        :param int new_compound: ontology identifier of the new compound.
        :return:
        """
        self.__new_metabolite = new_compound

    def set_compound_in_model(self, compound_in_model):
        """
        Method to set the compound to be replaced.
        :param int compound_in_model: ontology identifier of the compound in the model.
        :return:
        """
        self.__metabolite_in_model = compound_in_model

    def get_universal_model(self):
        """
        Get the instance variable __universal_model.
        :return: Model __universal_model: return the universal model.
        """
        return self.__universal_model

    def swap_metabolites(self) -> list:
        """
        Algorithm to swap metabolites:
            1º Set the game changer: dictionary with the metabolites' intermediates, conjugated acid or/and base.
                - {Metabolite in model : metabolite replacer}
            2º Swap all the metabolite's intermediates, conjugated acid or/and base and the metabolite itself.
            3º Swap the reactions associated to the replaced metabolite using the ReactionSwapper class.

        :return list : changed reactions
        """

        compounds_with_reactions_to_change = []
        game_changer = self.__get_game_changer()

        game_changer_in_model = {}
        for compound in game_changer.keys():
            compound_container = self.__compounds_ontology.get_node_by_ont_id(compound)
            in_model_compound_inchikey = compound_container.inchikey

            modelseed_id = compound_container.model_seed_id
            in_model_aliases = {}
            if modelseed_id:
                in_model_aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(modelseed_id)

            metabolitesInModel = self.mapper.check_metabolites_in_model(in_model_compound_inchikey,
                                                                        in_model_aliases,
                                                                        boimmg_container=compound_container)

            if metabolitesInModel:
                for metabolite in metabolitesInModel:
                    game_changer_in_model[metabolite] = game_changer[compound]

        for metabolite in game_changer_in_model.keys():

            inModel = game_changer_in_model[metabolite][0]
            if not inModel:
                new_compound_ontology_id = game_changer_in_model[metabolite][1]

                if self.__model_database == "ModelSEED":
                    self.__swap_modelseed_compound(metabolite, new_compound_ontology_id)

                elif self.__model_database == "BiGG":
                    self.__swap_bigg_metabolite(metabolite, new_compound_ontology_id)

                elif self.__model_database == "KEGG":
                    self.__swap_kegg_metabolite(metabolite, new_compound_ontology_id)

                compounds_with_reactions_to_change.append(metabolite)

            else:
                mets = game_changer_in_model[metabolite][1]
                for met in mets:
                    if metabolite.compartment == met.compartment:
                        change_reaction_by_swapping_existing_metabolites(metabolite, met)
                        compounds_with_reactions_to_change.append(met)

        changed_reactions = self.reactions_swapper.swap_reactions_by_compounds(compounds_with_reactions_to_change)
        return changed_reactions

    def __swap_conjugates(self, model_metabolite=None, new_metabolite=None):
        """
        Method to swap the conjugated base and acid of the compounds being swapped. If the model_metabolite and the
        new_metabolite are provided, the method will replace directly the model_metabolite for the new_metabolite.

        :param (optional) Metabolite model_metabolite: metabolite in the model to be replaced
        :param (optional) Metabolite new_metabolite:  metabolite to replace.
        :return:
        """
        conjugates_game_changer = self.__get_game_changer_conjugates()

        if not new_metabolite and not model_metabolite:

            compounds_with_reactions_to_change = []

            game_changer_in_model = {}
            for compound in conjugates_game_changer.keys():
                compound_container = self.__compounds_ontology.get_node_by_ont_id(compound)
                in_model_compound_inchikey = compound_container.inchikey
                in_model_aliases = compound_container.aliases
                metabolitesInModel = self.mapper.check_metabolites_in_model(in_model_compound_inchikey,
                                                                            in_model_aliases)
                if metabolitesInModel:
                    for metabolite_to_be_changed in metabolitesInModel:
                        game_changer_in_model[metabolite_to_be_changed] = conjugates_game_changer[compound]

            for metabolite_to_be_changed in game_changer_in_model.keys():

                inModel = game_changer_in_model[metabolite_to_be_changed][0]
                if not inModel:
                    new_compound_ontology_id = game_changer_in_model[metabolite_to_be_changed][1]

                    if self.__model_database == "ModelSEED":
                        self.__swap_modelseed_compound(metabolite_to_be_changed, new_compound_ontology_id)

                    elif self.__model_database == "BiGG":
                        self.__swap_bigg_metabolite(metabolite_to_be_changed, new_compound_ontology_id)

                    elif self.__model_database == "KEGG":
                        self.__swap_kegg_metabolite(metabolite_to_be_changed, new_compound_ontology_id)

                    compounds_with_reactions_to_change.append(metabolite_to_be_changed)

            self.reactions_swapper.swap_reactions_by_compounds(compounds_with_reactions_to_change)

        else:
            reactions_to_change = []
            for compound in conjugates_game_changer.keys():
                conjugates = conjugates_game_changer[compound][1]
                for conjugate in conjugates:
                    if conjugate.compartment == new_metabolite.compartment:

                        compound_container = self.__compounds_ontology.get_node_by_ont_id(compound)

                        in_model_compound_inchikey = compound_container.inchikey
                        in_model_aliases = compound_container.aliases

                        metabolitesInModel = self.mapper.check_metabolites_in_model(in_model_compound_inchikey,
                                                                                    in_model_aliases,
                                                                                    compound_container)

                        for metabolite_to_be_changed in metabolitesInModel:
                            if metabolite_to_be_changed.compartment == new_metabolite.compartment:
                                reactions_to_change = change_reaction_by_swapping_existing_metabolites(
                                    metabolite_to_be_changed, conjugate)

                    else:
                        db_type = self.__compoundsAnnotationConfigs.get(self.__model_database)
                        db_id = conjugate.annotation.get(db_type)

                        if type(db_id) == str:
                            modelseed_id = self.__compoundsIdConverter.convert_db_id_to_model_seed_by_db_id(db_id)
                        else:
                            modelseed_id = self.__compoundsIdConverter.convert_db_id_to_model_seed_by_db_id(db_id[0])

                        cont = self.__compounds_db.get_compound_by_id(modelseed_id[0])
                        new_compound = self.generate_new_metabolites(cont)

                        compound_ont_container = self.__compounds_ontology.get_node_by_ont_id(compound)

                        compound_to_be_changed = \
                            self.mapper.check_metabolites_in_model(compound_ont_container.inchikey,
                                                                   compound_ont_container.aliases,
                                                                   compound_ont_container)

                        if compound_to_be_changed:
                            reactions_to_change = change_reaction_by_swapping_existing_metabolites(
                                compound_to_be_changed[0], new_compound[0])

            if reactions_to_change:
                self.reactions_swapper.swap_reactions(reactions_to_change)

    def swap_only_metabolites_and_conjugates(self, model_metabolite=None, new_metabolite=None):
        """
        Method to swap only the metabolites and their conjugated acid and/or base. This method does not change the
        metabolite intermediates/precursors.

        :param (optional) Model.Metabolite model_metabolite: metabolite in the model to be replaced
        :param (optional) Model.Metabolite new_metabolite:  metabolite to replace.
        :return:
        """

        compounds_with_reactions_to_change = []

        if not model_metabolite:
            self.__swap_conjugates()

            ontology_container = self.__compounds_ontology.get_node_by_ont_id(self.__metabolite_in_model)
            inchikey = ontology_container.inchikey
            aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(ontology_container.model_seed_id)
            model_metabolites = self.mapper.check_metabolites_in_model(inchikey, aliases)

            replacer_metabolites_in_model = \
                self.__check_if_compound_exists_in_model_by_ontology_id(self.__new_metabolite)

            if not replacer_metabolites_in_model:

                if self.__model_database == "ModelSEED":
                    for model_metabolite in model_metabolites:
                        self.__swap_modelseed_compound(model_metabolite, self.__new_metabolite)

                elif self.__model_database == "BiGG":
                    for model_metabolite in model_metabolites:
                        self.__swap_bigg_metabolite(model_metabolite, self.__new_metabolite)

                elif self.__model_database == "KEGG":
                    for model_metabolite in model_metabolites:
                        self.__swap_kegg_metabolite(model_metabolite, self.__new_metabolite)

                compounds_with_reactions_to_change.append(model_metabolite)

                self.reactions_swapper.swap_reactions_by_compounds(compounds_with_reactions_to_change)

            else:

                for model_metabolite in model_metabolites:
                    for replacer_metabolite in replacer_metabolites_in_model:
                        if model_metabolite.compartment == replacer_metabolite.compartment:
                            reactions_to_change = change_reaction_by_swapping_existing_metabolites(
                                model_metabolite, replacer_metabolite)

                            self.reactions_swapper.swap_reactions(reactions_to_change)

        elif new_metabolite:
            self.__swap_conjugates(model_metabolite, new_metabolite)

            reactions_to_change = change_reaction_by_swapping_existing_metabolites(
                model_metabolite, new_metabolite)

            changed_reactions = self.reactions_swapper.swap_reactions(reactions_to_change)
            self.changed_reactions.extend(changed_reactions)

    def __swap_bigg_metabolite(self, model_metabolite, new_compound_ontology_id):
        """
        Method to swap a metabolite in the BiGG database format.
        This method tries to find the BiGG format of the new compound. If it does not find, it will change into the
        Model SEED format.

        :param Model.Metabolite model_metabolite: metabolite to be replaced.
        :param int new_compound_ontology_id: ontology identifier of the metabolite that will replace the other
        :return:
        """

        model_metabolite.annotation = {}
        compound_container = self.__compounds_ontology.get_node_by_ont_id(new_compound_ontology_id)
        model_seed_id = compound_container.model_seed_id

        if model_seed_id:
            bigg_ids = None
            aliases = self.__compoundsIdConverter.get_modelSeedIdToDb().get(model_seed_id)
            if "BiGG" in aliases:
                bigg_ids = self.__compoundsIdConverter.convert_modelSeedId_into_other_dbID(model_seed_id, "BiGG")
            elif "BiGG1" in aliases:
                bigg_ids = self.__compoundsIdConverter.convert_modelSeedId_into_other_dbID(model_seed_id, "BiGG1")

            if bigg_ids:

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

                    old_inchikey = ""
                    if "inchi_key" in model_metabolite.annotation.keys():
                        old_inchikey = model_metabolite.annotation["inchi_key"]

                    old_aliases = AnnotationUtils.get_annotation_from_cobra_annotation(model_metabolite)
                    old_id = model_metabolite.id

                    new_aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(
                        compound_container.model_seed_id)
                    annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(new_aliases)
                    model_metabolite.annotation = annotation

                    new_inchikey = ""
                    if not compound_container.generic:
                        new_inchikey = compound_container.inchikey
                        model_metabolite.annotation["inchi_key"] = compound_container.inchikey

                    model_metabolite.annotation["smiles"] = compound_container.smiles
                    model_metabolite.formula = bigg_metabolite.get("formulae")[0]

                    charges = bigg_metabolite.get("charges")
                    if charges:
                        model_metabolite.charge = charges[0]

                    compartment = model_metabolite.compartment

                    if compartment is None:
                        model_metabolite.compartment = "c"

                    model_metabolite.name = bigg_metabolite.get("name")

                    self.__check_if_id_is_used_in_model_and_delete_it(bigg_id + "_" + model_metabolite.compartment)

                    self.report_material[model_metabolite.id] = bigg_id + "_" + model_metabolite.compartment
                    model_metabolite.id = bigg_id + "_" + model_metabolite.compartment

                    self.mapper.update_maps(old_inchikey, new_inchikey, old_id,
                                            model_metabolite.id, new_compound_ontology_id,
                                            old_aliases, new_aliases)

                else:
                    self.__change_boimmg_format_metabolite(model_metabolite, compound_container)

            else:
                self.__change_boimmg_format_metabolite(model_metabolite, compound_container)

        else:
            self.__change_boimmg_format_metabolite(model_metabolite, compound_container)

    def __swap_kegg_metabolite(self, model_metabolite, new_compound_ontology_id):
        """
        Method to swap a metabolite in the KEGG database format.
        This method tries to find the KEGG format of the new compound. If it does not find, it will change into the
        Model SEED format.

        :param Model.Metabolite model_metabolite: metabolite to be replaced.
        :param int new_compound_ontology_id: ontology identifier of the metabolite that will replace the other
        :return:
        """

        old_inchikey = ""
        if "inchi_key" in model_metabolite.annotation.keys():
            old_inchikey = model_metabolite.annotation["inchi_key"]

        old_aliases = AnnotationUtils.get_annotation_from_cobra_annotation(model_metabolite)
        old_id = model_metabolite.id

        model_metabolite.annotation = {}
        compound_container = self.__compounds_ontology.get_node_by_ont_id(new_compound_ontology_id)
        model_seed_id = compound_container.model_seed_id
        new_aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(model_seed_id)

        if model_seed_id in self.__compoundsIdConverter.get_modelSeedIdToDb().keys():

            if "KEGG" in self.__compoundsIdConverter.get_modelSeedIdToDb().get(model_seed_id):

                kegg_id = self.__compoundsIdConverter.convert_modelSeedId_into_other_dbID(model_seed_id, "KEGG")[0]
                kegg_metabolite = KeggCompound(kegg_id)

                annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(new_aliases)
                model_metabolite.annotation = annotation

                new_inchikey = ""
                if not compound_container.generic:
                    model_metabolite.annotation["inchi_key"] = compound_container.inchikey
                    new_inchikey = compound_container.inchikey

                model_metabolite.annotation["smiles"] = compound_container.smiles
                model_metabolite.formula = kegg_metabolite.get_formula()
                model_metabolite.name = kegg_metabolite.get_name()
                model_metabolite.charge = 0
                compartment = model_metabolite.compartment

                self.__check_if_id_is_used_in_model_and_delete_it(kegg_id + "_" + compartment)

                self.report_material[model_metabolite.id] = kegg_id + "_" + compartment

                model_metabolite.id = kegg_id + "_" + compartment

                self.mapper.update_maps(old_inchikey, new_inchikey, old_id,
                                        model_metabolite.id, compound_container.id,
                                        old_aliases, new_aliases)

            else:
                self.__change_boimmg_format_metabolite(model_metabolite, compound_container)

        else:
            self.__change_boimmg_format_metabolite(model_metabolite, compound_container)

    def __change_boimmg_format_metabolite(self, model_metabolite: Metabolite, compound_container):
        """
        This method changes the id of the metabolite arbitrarily. The format of the metabolite is set for the
        ontology format.
        It is worth noting that the :param compound_container can be either a ModelSeedCompound or an OntologyCompound

        :param Model.Metabolite model_metabolite: metabolite to be replaced.
        :param CompoundContainer compound_container: compound that will replace the other.
        :return:
        """

        old_inchikey = ""
        if "inchi_key" in model_metabolite.annotation.keys():
            old_inchikey = model_metabolite.annotation["inchi_key"]

        old_aliases = AnnotationUtils.get_annotation_from_cobra_annotation(model_metabolite)
        old_id = model_metabolite.id

        aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(compound_container.model_seed_id)
        annotation = {}
        if aliases:
            annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(aliases)

        annotation["boimmg.compound"] = \
            self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + str(compound_container.id)

        new_aliases = {"BOIMMG": [annotation["boimmg.compound"]]}
        model_metabolite.annotation = annotation

        new_inchikey = ""
        if not compound_container.generic:
            model_metabolite.annotation["inchi_key"] = compound_container.inchikey
            new_inchikey = compound_container.inchikey

        model_metabolite.annotation["smiles"] = compound_container.smiles
        model_metabolite.formula = compound_container.formula
        model_metabolite.charge = compound_container.charge
        model_metabolite.name = compound_container.name
        compartment = model_metabolite.compartment

        metabolite_id = self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"] + str(compound_container.id)

        new_metabolite_id = metabolite_id + "_" + compartment

        self.__check_if_id_is_used_in_model_and_delete_it(new_metabolite_id)

        self.report_material[model_metabolite.id] = new_metabolite_id
        model_metabolite.id = new_metabolite_id

        self.mapper.update_maps(old_inchikey, new_inchikey, old_id,
                                model_metabolite.id, compound_container.id,
                                old_aliases, new_aliases)

    def __swap_modelseed_compound(self, model_metabolite: Metabolite, new_compound_ontology_id: int):

        """
        Method to swap a metabolite in the Model SEED database format.
        This method tries to find the Model SEED format of the new compound. If it does not find, it will change into an
        arbitrary id with the Model SEED or the Ontology format.

        :param Model.Metabolite model_metabolite: metabolite to be replaced.
        :param int new_compound_ontology_id: ontology identifier of the metabolite that will replace the other
        :return:
        """

        old_inchikey = ""
        if "inchi_key" in model_metabolite.annotation.keys():
            old_inchikey = model_metabolite.annotation["inchi_key"]

        old_aliases = AnnotationUtils.get_annotation_from_cobra_annotation(model_metabolite)
        old_id = model_metabolite.id

        model_metabolite.annotation = {}
        compound_container = self.__compounds_ontology.get_node_by_ont_id(new_compound_ontology_id)

        new_aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(compound_container.model_seed_id)
        annotation = AnnotationUtils.get_compound_annotation_format_by_aliases(new_aliases)

        model_metabolite.annotation = annotation

        new_inchikey = ""
        if not compound_container.generic:
            model_metabolite.annotation["inchi_key"] = compound_container.inchikey
            new_inchikey = compound_container.inchikey

        model_metabolite.annotation["smiles"] = compound_container.smiles
        model_metabolite.formula = compound_container.formula
        model_metabolite.charge = compound_container.charge
        compartment = model_metabolite.compartment
        db_id = compound_container.model_seed_id
        model_metabolite.name = compound_container.name

        if db_id:

            self.__check_if_id_is_used_in_model_and_delete_it(compound_container.model_seed_id + "_" + compartment)

            self.report_material[model_metabolite.id] = compound_container.model_seed_id + "_" + compartment

            model_metabolite.id = compound_container.model_seed_id + "_" + compartment

            self.mapper.update_maps(old_inchikey, new_inchikey, old_id,
                                    model_metabolite.id, new_compound_ontology_id,
                                    old_aliases, new_aliases)
        else:
            self.__change_boimmg_format_metabolite(model_metabolite, compound_container)

    def __check_if_id_is_used_in_model_and_delete_it(self, model_id: int):
        """

        This method checks whether a compound ID is used in the model and if so, it is removed

        :param model_id: compound identifier
        :return:
        """

        inModel = self.model.metabolites.has_id(model_id)

        if inModel:
            metabolite = self.model.metabolites.get_by_id(model_id)
            self.model.remove_metabolites([metabolite])

    def __get_game_changer_conjugates(self):
        """
        This method searches for the conjugates of a given compound in model and assigns their conjugated acid or base;

        :return dict {ontology id : (metabolite already in model ? (boolean), COBRApy metabolite or ontology id)}
        """

        conjugated_acid_in_model = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
            self.__metabolite_in_model,
            "conjugated_acid_of"
        )
        conjugated_base_in_model = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
            self.__metabolite_in_model,
            "conjugated_base_of"
        )

        game_changer = {}

        if conjugated_base_in_model:
            conjugated_base_to_change = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
                self.__new_metabolite,
                "conjugated_base_of"
            )
            game_changer.update(
                self.__add_conjugate_to_game_changer("base", conjugated_base_in_model, conjugated_base_to_change))

        if conjugated_acid_in_model:
            conjugated_acid_to_change = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
                self.__new_metabolite,
                "conjugated_acid_of"
            )
            game_changer.update(
                self.__add_conjugate_to_game_changer("acid", conjugated_acid_in_model, conjugated_acid_to_change))

        return game_changer

    def __get_game_changer(self):

        """
        This method set the game changer. The game changer will be a dictionary with all of the necessary
        information to swap the compounds. The game changer will have a given compound intermediates, conjugated acid
        or/and base and the respective metabolites to replace them.

        :return dict: {ontology id : (metabolite already in model ? (boolean), COBRApy metabolite or ontology id)}
        """

        intermediates_in_model = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(
            self.__metabolite_in_model,
            "precursor_of")

        conjugated_acid_in_model = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
            self.__metabolite_in_model,
            "conjugated_acid_of"
        )
        conjugated_base_in_model = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
            self.__metabolite_in_model,
            "conjugated_base_of"
        )

        intermediates_to_change = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(
            self.__new_metabolite, "precursor_of")

        #  values   ----->  (met in model?, if true: cobraMetabolite else: ontId)
        game_changer = self.__add_to_game_changer(self.__metabolite_in_model, self.__new_metabolite)

        game_changer.update(self.__get_intermediates_game_changer(intermediates_in_model, intermediates_to_change))

        if conjugated_acid_in_model:
            conjugated_acid_to_change = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
                self.__new_metabolite,
                "conjugated_acid_of"
            )
            game_changer.update(
                self.__add_conjugate_to_game_changer("acid", conjugated_acid_in_model, conjugated_acid_to_change))

            intermediates_in_model = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(
                conjugated_acid_in_model[0],
                "precursor_of")

            intermediates_to_change = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(
                conjugated_acid_to_change[0],
                "precursor_of")

            game_changer.update(self.__get_intermediates_game_changer(intermediates_in_model, intermediates_to_change))

        if conjugated_base_in_model:
            conjugated_base_to_change = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
                self.__new_metabolite,
                "conjugated_base_of"
            )
            game_changer.update(
                self.__add_conjugate_to_game_changer("base", conjugated_base_in_model, conjugated_base_to_change))

            intermediates_in_model = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(
                conjugated_base_in_model[0],
                "precursor_of")

            intermediates_to_change = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(
                conjugated_base_to_change[0],
                "precursor_of")

            game_changer.update(self.__get_intermediates_game_changer(intermediates_in_model, intermediates_to_change))

        return game_changer

    def __add_conjugate_to_game_changer(self, change_type: str, conjugated_in_model: list, new_conjugated: list) \
            -> dict:
        """
        This method add a conjugated acid or base to the game changer

        :param string change_type: "acid" or "base"
        :param list conjugated_in_model: ontology id of the conjugate to be replaced
        :param list new_conjugated: ontology id of the conjugate that will replace

        :return dict: {ontology id : (metabolite already in model ? (boolean), COBRApy metabolite or ontology id)}
        """

        game_changer = {}

        if change_type == "base":

            metabolite_added = self.__add_to_game_changer(conjugated_in_model[0], new_conjugated[0])
            game_changer.update(metabolite_added)

        elif change_type == "acid":

            metabolite_added = self.__add_to_game_changer(conjugated_in_model[0], new_conjugated[0])
            game_changer.update(metabolite_added)

        return game_changer

    def __add_to_game_changer(self, met_in_model, new_met):
        """
        This method add metabolites to the game changer

        :param int met_in_model: ontology id of the metabolite to be replaced
        :param int new_met: ontology id of the metabolite to be replaced

        :return dict: {ontology id : (metabolite already in model ? (boolean), COBRApy metabolite or ontology id)}
        """
        game_changer = {}
        metabolite_container = self.__compounds_ontology.get_node_by_ont_id(new_met)
        inchikey = metabolite_container.inchikey
        aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(metabolite_container.model_seed_id)
        metabolites_check_in_model = self.mapper.check_metabolites_in_model(inchikey, aliases)

        if met_in_model not in self.__metabolites_not_to_change:
            if metabolites_check_in_model:
                game_changer[met_in_model] = (True, [])
                for metabolite in metabolites_check_in_model:
                    game_changer[met_in_model][1].append(
                        metabolite)  # (met in model?, if true: [cobraMetabolite] else: ontId)
            else:
                game_changer[met_in_model] = (False, new_met)

            return game_changer
        else:
            if not metabolites_check_in_model:

                if metabolite_container.model_seed_id:

                    ms_container = self.__compounds_db.get_compound_by_id(metabolite_container.model_seed_id)
                    metabolites = self.generate_new_metabolites(ms_container)

                else:

                    metabolites = self.generate_new_boimmg_metabolites(metabolite_container)

                self.model.add_metabolites(metabolites)

            return game_changer

    def _check_if_compound_is_generic(self, compound_id):
        """
        This method check if some compound is generic or not

        :param int compound_id: ontology id of the metabolite

        :return boolean: is generic or not
        """

        compound_container = self.__compounds_ontology.get_node_by_ont_id(compound_id)
        return compound_container.generic

    def __get_intermediates_game_changer(self, intermediates_in_model, intermediates_to_change):
        """
        This method finds the correspondent substitute of each intermediate in the model

        :param list intermediates_in_model: list of ontology ids of the intermediates of the metabolite to be changed
        :param list intermediates_to_change: list of ontology ids of the intermediates of the metabolite that will
        replace

        :return dict: {ontology id : (metabolite already in model? (boolean), COBRApy metabolite or ontology id) }
        """

        game_changer = {}
        if self.__type == 1:
            intermediates = self.__set_intermediates_in_game_changer_type1(intermediates_in_model,
                                                                           intermediates_to_change)
            game_changer.update(intermediates)

        elif self.__type == 2:
            intermediates = self.__set_intermediates_in_game_changer_type2(intermediates_in_model,
                                                                           intermediates_to_change)
            game_changer.update(intermediates)

        return game_changer

    def __set_intermediates_in_game_changer_type1(self, intermediates_in_model, intermediates_to_change):
        """
        This method finds the correspondent substitute of each intermediate in the model for the type 1 game change

        :param list intermediates_in_model: list of ontology ids of the intermediates of the metabolite to be changed
        :param list intermediates_to_change: list of ontology ids of the intermediates of the metabolite that will
        replace

        :return dict: {ontology id : (metabolite already in model? (boolean), COBRApy metabolite or ontology id) }
        """

        game_changer = {}
        for met in intermediates_in_model:
            parent_met_in_model = self.__compounds_ontology.get_successors_by_ont_id_rel_type(met, "is_a")
            i = 0
            found = False
            while not found and i < len(intermediates_to_change):
                parent2 = self.__compounds_ontology.get_successors_by_ont_id_rel_type(intermediates_to_change[i],
                                                                                      "is_a")
                if parent_met_in_model[0] == parent2[0]:
                    found = True

                    metabolite_added = self.__add_to_game_changer(met, intermediates_to_change[i])
                    game_changer.update(metabolite_added)

                i += 1
        return game_changer

    def __set_intermediates_in_game_changer_type2(self, intermediates_in_model, intermediates_to_change):
        """
        This method finds the correspondent substitute of each intermediate in the model for the type 2 game change

        :param list intermediates_in_model: list of ontology ids of the intermediates of the metabolite to be changed
        :param list intermediates_to_change: list of ontology ids of the intermediates of the metabolite that will
        replace

        :return dict: {ontology id : (metabolite already in model? (boolean), COBRApy metabolite or ontology id) }
        """

        game_changer = {}
        isGenericInModel = self._check_if_compound_is_generic(intermediates_in_model[0])
        self.reactions_swapper.set_type2_is_generic_in_model(isGenericInModel)

        if isGenericInModel:
            for met in intermediates_in_model:
                i = 0
                found = False
                while not found and i < len(intermediates_to_change):
                    parent_met_to_change = self.__compounds_ontology.get_successors_by_ont_id_rel_type(
                        intermediates_to_change[i], "is_a")

                    if met in parent_met_to_change:
                        game_changer.update(self.__add_to_game_changer(met, intermediates_to_change[i]))
                        found = True
                    i += 1
        else:
            for met in intermediates_in_model:
                i = 0
                found = False
                while not found and i < len(intermediates_to_change):
                    children = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(
                        intermediates_to_change[i], "is_a")

                    if met in children:
                        game_changer.update(self.__add_to_game_changer(met, intermediates_to_change[i]))
                        found = True
                    i += 1
        return game_changer

    def __get_metabolite_container_by_id(self, metabolite_id):
        """
        This method searches for a given metabolite in the Model

        :param string metabolite_id: id in the Model

        :return Model.metabolite: if the metabolite is found returns the metabolite, otherwise returns None
        """

        for met in self.model.metabolites:
            if met.id == metabolite_id:
                return met
        return None

    def change_reaction_format(self, reaction):
        """
        This method changes a given reaction format.

        :param Model.Reaction reaction: reaction to be changed
        :return:
        """

        self.reactions_swapper.swap_reactions([reaction])

    def __check_if_compound_exists_in_model_by_ontology_id(self, ontology_id):

        """
        Method to check whether a compound exists in the model using BOIMMG ID

        :param int ontology_id: boimmg id :return list<Compound>: list of compounds in model with this ID (list of
        the same compound in different compartments)

        """

        container = self.__compounds_ontology.get_node_by_ont_id(ontology_id)
        inchikey = container.inchikey
        aliases = {}
        if container.model_seed_id:
            aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(container.model_seed_id)

        aliases.update(container.aliases)
        return self.mapper.check_metabolites_in_model(inchikey, aliases, container)

    def generate_new_metabolites(self, compound_container: ModelSeedCompound):
        """
        This method will generate a new set of compounds (same metabolite in different compartments)

        :param ModelSeedCompound compound_container: internal wrappers of the ModelSEED compounds
        :return list<Metabolite>: list of newly generated compounds
        """

        model_metabolites = \
            model_utilities.generate_model_compounds_by_database_format(self.model,
                                                                        compound_container.getDbId(),
                                                                        self.__compoundsIdConverter,
                                                                        self.__compounds_db,
                                                                        self.__model_database)

        self.model.add_metabolites(model_metabolites)
        self.mapper.add_new_metabolites_to_maps(model_metabolites)
        return model_metabolites

    def generate_new_boimmg_metabolites(self, metabolite_container: CompoundNode):

        """
        This method will generate a new set of compounds (same metabolite in different compartments)

        :param CompoundNode metabolite_container: wrapper of BOIMMG compound
        :return list<Metabolite>: list of newly generated compounds
        """

        metabolites = model_utilities.generate_boimmg_metabolites(self.model,
                                                                  metabolite_container,
                                                                  self.__model_database,
                                                                  self.__compoundsIdConverter,
                                                                  self.__compoundsAnnotationConfigs)

        self.model.add_metabolites(metabolites)

        self.mapper.add_new_metabolites_to_maps(metabolites)

        return metabolites

    def setType2_isGenericInModel(self, isGeneric):
        self.reactions_swapper.set_type2_is_generic_in_model(isGeneric)

    def set_model_mapper(self, mapper):
        self.mapper = mapper
