from cobra import Model, Reaction

from boimmgpy.service.revisor.compounds_revisor import CompoundsRevisor
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.utilities import file_utilities
from boimmgpy.definitions import TOOL_CONFIG_PATH, ROOT_DIR, REACTIONS_ANNOTATION_CONFIGS_PATH, \
    COMPOUNDS_ANNOTATION_CONFIGS_PATH


class ReactionsChanger:

    def __init__(self, model: Model, swap_type: int, model_database: str,
                 universal_model: Model, compounds_converter=None):
        """
        Class constructor

        :param cobrapy.Model model: COBRApy model
        :param int swap_type: type of change (0,1,2 or 3)
        :param string model_database: database format of metabolites and reactions
        :param cobrapy.Model universal_model: universal model
        :param compounds_converter CompoundsIDConverter: object to convert metabolites identifiers
        """

        if not compounds_converter:
            self.__compoundsIDConverter = CompoundsIDConverter()

        else:
            self.__compoundsIDConverter = compounds_converter

        if not universal_model:
            self.__universal_model = Model()
        else:
            self.__universal_model = universal_model

        self.report_material = {}
        self.__type = swap_type
        self.__model_database = model_database
        self.__model = model
        self.__modelseed_hydrogen = "cpd00067"
        self.__reactionBalancer = CompoundsRevisor(self.__model)
        self.__not_found_reactions_num = 0

        self.__home_path__ = ROOT_DIR

        self.__configs = file_utilities.read_conf_file(
            TOOL_CONFIG_PATH)
        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

    @property
    def model(self):
        return self.__model

    @model.setter
    def model(self, value):
        self.__model = value

    def swap_reactions_by_compounds(self, compounds_with_reactions_to_change: list) -> list:
        """
        This method swap all the reactions associated to a list of compounds.
        Firstly, it tries to balance the reactions and then tries to change the reaction.

        :param list compounds_with_reactions_to_change: list of Model.metabolites with reactions to change.
        :return list: reactions changed
        """
        past = []
        reaction_changed = []
        for compound in compounds_with_reactions_to_change:
            reactions = self.get_reactions_by_compound(compound)
            for reaction in reactions:
                if reaction not in past and reaction.annotation.get("sbo") != "SBO:0000629" \
                        and "Biomass" not in reaction.id:
                    past.append(reaction)

                    self.change_reaction(reaction)

                    reaction_changed.append(reaction)

        return reaction_changed

    def get_reactions_by_compound(self, compound):
        """
        Return all the reactions associated with a given compound

        :param Metabolite compound:
        :return list<Reaction>: all the reactions associated with a given compound
        """
        res = []
        for reaction in self.__model.reactions:
            metabolites_ids = [m.id for m in reaction.metabolites]

            if compound.id in metabolites_ids:
                res.append(reaction)

        return res

    def swap_reactions(self, reactions_to_change) -> list:
        """
        This method swap all the reactions associated to a list of compounds.
        Firstly, it tries to balance the reactions and then tries to change the reaction.

        :param list reactions_to_change: list of Model.reactions with reactions to change.
        :return list: list of changed reactions
        """
        past = []
        reaction_changed = []
        for reaction in reactions_to_change:
            if reaction not in past and reaction.annotation.get(
                    "sbo") != "SBO:0000629" and "Biomass" not in reaction.id:
                past.append(reaction.copy())

                self.change_reaction(reaction)

                reaction_changed.append(reaction)

        return reaction_changed

    def change_reaction(self, reaction):
        """
        This method changes a given reaction's name, id, delta G, etc...
        It is worth noting that the reaction in the model is no longer the true one, since the metabolites were swapped
        1º start by searching in the ontology whether there is any entity representing the reaction in the model
        2º - if the reaction is found then other reaction in the ontology is searched to replace it
        3º - otherwise the reaction will be searched in a reaction database
        4º - if the reaction was not found, then a cannonical name and id will be set.

        :param Model.reaction reaction: reaction to be changed.

        """

        # if self.__model_database in self.__reactionsAnnotationConfigs:
        #     database = self.__reactionsAnnotationConfigs[self.__model_database]
        # else:
        #     raise Exception

        self.__change_not_found_reaction(reaction)

    def __change_not_found_reaction(self, reaction):
        """
        This method will change a not found reaction in the database

        :param Model.reaction reaction: reaction to be changed.

        """

        self.change_boimmg_reaction_id(reaction)

        reaction.name = "Changed - " + reaction.name
        reaction.annotation = {}

    def change_boimmg_reaction_id(self, reaction: Reaction):

        new_id = self.__reactionsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]
        metabolites = reaction.metabolites

        metabolites = sorted([metabolite.id for metabolite in metabolites])
        new_id += "_".join(metabolites)

        self.report_material[reaction.id] = new_id
        reaction.id = new_id

        return reaction.id

    def get_database_ids_compounds_in_model(self, compounds_in_model):
        """
        This method searchs in a given database the ids of the compounds in the model

        :param dict compounds_in_model: model compounds assigned to a model reaction

        :returns list: list with database ids of the metabolites assigned to a model reaction

        """

        res = []
        for compound in compounds_in_model.keys():
            modelseed_id_annotation = compound.annotation.get("seed.compound")

            if modelseed_id_annotation and type(modelseed_id_annotation) == str:
                modelseed_ids = [modelseed_id_annotation]

            else:
                modelseed_ids = modelseed_id_annotation

            if self.__model_database != "ModelSEED" and not modelseed_ids:
                db_id = compound.annotation.get(self.__compoundsAnnotationConfigs[self.__model_database])
                if type(db_id) == list:

                    i = 0
                    found = False

                    while not found and i < len(db_id):
                        modelseed_ids = self.__compoundsIDConverter.convert_dbID_into_modelSeedId(self.__model_database,
                                                                                                  db_id[i])

                        if modelseed_ids:
                            found = True
                else:
                    modelseed_ids = self.__compoundsIDConverter.convert_dbID_into_modelSeedId(self.__model_database,
                                                                                              db_id)

            if not modelseed_ids:
                return []

            else:
                if not res:
                    for modelseed_id in modelseed_ids:
                        res.append([modelseed_id])
                else:
                    if len(modelseed_ids) == 1:
                        for reaction in res:
                            reaction.append(modelseed_ids[0])
                    else:
                        previous = res.copy()
                        res = []
                        for modelseed_id in modelseed_ids:
                            for reaction in previous:
                                new_reaction = reaction.copy()
                                new_reaction.append(modelseed_id)
                                res.append(new_reaction)

        final_res = []
        for reaction in res:
            final_res.append(sorted(reaction))
        return final_res

    def set_type2_is_generic_in_model(self, isGenericInModel):
        """
        This method sets the state of the metabolite in the model.

        :param boolean isGenericInModel: if True -> Generic; if False -> Complete
        :return:
        """
        self.__set_type2_is_generic_in_model = isGenericInModel

    def set_type(self, new_type):

        """
        This method sets the type of swap.
        :param int new_type: Type of swap (0,1,2,3)
        :return:
        """

        self.__type = new_type
