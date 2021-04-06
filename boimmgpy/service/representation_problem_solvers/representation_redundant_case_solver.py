import re
from copy import deepcopy

from biocyc import biocyc
from cobra import Model, Reaction, Metabolite
from cobra.util import linear_reaction_coefficients


from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.service.network_modifiers.granulator import Granulator
from boimmgpy.service.revisor.compounds_revisor import CompoundsRevisor
from boimmgpy.service.model_mapper import ModelMapper
from boimmgpy.service.interfaces.representation_problem_solver import RepresentationProblemSolver
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from boimmgpy.utilities import file_utilities
from boimmgpy.definitions import ROOT_DIR, TOOL_CONFIG_PATH, COMPOUNDS_ANNOTATION_CONFIGS_PATH, REACTIONS_ANNOTATION_CONFIGS_PATH


PROGRESS_BAR = ROOT_DIR + "/service/logs/progress_bar.txt"
# PROGRESS_BAR = "/workdir/resultsWorker/progress_bar.txt"

class RedundantCaseSolver(RepresentationProblemSolver):

    def __init__(self, model, database_format):
        """
        Class constructor

        :param Model model: cobrapy model
        :param str database_format: ModelSEED, BiGG or KEGG
        """

        self.model = model

        self.objective = linear_reaction_coefficients(self.model)
        if self.objective:
            self.objective = list(self.objective.keys())[0]
        else:
            raise Exception("No objective found")

        self.__compounds_ontology = CompoundsDBAccessor()

        self.model = model
        self.__database_format = database_format
        self.__configs = file_utilities.read_conf_file(TOOL_CONFIG_PATH)
        self.__home_path__ = ROOT_DIR


        self.__modelseedCompoundsDb = ModelSeedCompoundsDB()

        self.__define_instance_variables()


        self.write_in_progress_bar("mapping model... ",1)
        self.mapper.map_model(database_format)
        self.write_in_progress_bar("model mapped ", 10)


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
    def mapper(self,value):
        self.__mapper = value

    def write_in_progress_bar(self, message, value):
        with open(PROGRESS_BAR,"w") as file:
            file.write(message + ":" + str(value))


    def __define_instance_variables(self):
        """
        Define instance variables

        :return:
        """

        self.__universal_model = Model("universal_model")
        self.__swapped = []

        self.__compounds_ontology = CompoundsDBAccessor()

        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)


        self.__compoundsIdConverter = CompoundsIDConverter()

        self.__compounds_revisor = CompoundsRevisor(self.model, self.__universal_model)

        biocyc.set_organism("meta")

        self.mapper = ModelMapper(self.model,self.__compoundsAnnotationConfigs,self.__compoundsIdConverter)

        self.__compounds_revisor.set_model_mapper(self.mapper)

        self.granulator = Granulator(self.model, self.mapper, self.__database_format, self.__compoundsIdConverter,
                                     self.__compoundsAnnotationConfigs)

    def __add_hydrogen_to_virtual_model(self):
        """
        This method set searches for the hydrogen in each compartment of the model.
        The hydrogens are used to balance reactions.
        """


        compartments = self.__model.compartments
        metabolites_in_model = self.__model.metabolites

        res = []

        i = 0
        while i<len(metabolites_in_model):

            metabolite = metabolites_in_model[i]
            if metabolite.formula == "H" and metabolite.compartment in compartments:
                res.append(metabolite)

            i+=1

        self.__virtual_model.add_metabolites(res)

    def swap_and_gap_fill(self,target_ontology_id):
        pass

    def generalize_model(self,target_ontology_id):
        pass

    def gap_fill_model_by_target(self,target_ontology_id: int, components_ont_ids: list):
        pass

    def write_reports(self,file_path):

        metabolites_report_material = self.granulator.metabolite_report_material
        reactions_report_material = self.granulator.reaction_report_material

        with open(file_path, "w") as f:

            f.write("--Metabolites--\n")
            f.write("metabolite in model,replacer metabolite\n")

            for metabolite in metabolites_report_material:
                f.write(metabolite + ",")
                f.write(metabolites_report_material[metabolite] + "\n")

            f.write("--Reactions--\n")
            f.write("reaction in model,replacer reaction\n")

            for reaction in reactions_report_material:
                f.write(reaction + ",")
                f.write(reactions_report_material[reaction] + "\n")

    def swap_from_generic(self,targets: list, components: list,same_components = False,
                          progress_bar_processes_left = 1):
        """
        Swap compounds from a generic target and components. The algorithm acts as follows:

            1st: Starts for identifying the reactions present in the requested biosynthetic pathway;

            2nd: Creates a virtual model;

            3rd: Calls the granulator to granulate the biosynthetic pathway;

        :param list<str> targets: list of lipid targets
        :param list<str> components: list of components
        :param bool same_components: boolean for same or mix of components
        :param int progress_bar_processes_left:
        :return:
        """

        targets_generic_ontology_id, components_ont_ids = \
            self.transform_targets_and_components_into_boimmg_ids(targets,components)

        self.__virtual_model = Model("virtual_model")

        self.__add_hydrogen_to_virtual_model()

        self.__virtual_model_mapper = ModelMapper(self.__virtual_model, self.__compoundsAnnotationConfigs,
                                                self.__compoundsIdConverter)

        self.__virtual_compounds_revisor = CompoundsRevisor(self.__virtual_model, self.__universal_model,
                                                            self.__compoundsAnnotationConfigs)

        self.__virtual_compounds_revisor.set_model_mapper(self.__virtual_model_mapper)



        self.granulator.set_virtual_model(self.__virtual_model,self.__virtual_model_mapper,self.__virtual_compounds_revisor)


        for target_generic_ontology_id in targets_generic_ontology_id:


            biosynthesis_reactions = self.granulator.identify_biosynthesis_pathway(target_generic_ontology_id)


            self.__biosynthesis_reactions = deepcopy(biosynthesis_reactions)

            self.__add_new_reactions_to_virtual_model(self.__biosynthesis_reactions)

            self.__model.remove_reactions(biosynthesis_reactions)

            self.granulator.granulate(target_generic_ontology_id, components_ont_ids, same_components,
                                progress_bar_processes_left)

        self.generateISAreactions()

    def __subtract_reactions_in_list(self,reactions,reactions_to_subtract):

        res = []

        for reaction in reactions:
            if reaction.id not in reactions_to_subtract:
                res.append(reaction)

        return res

    # def solve_components_reactions(self,same_components=False):
    #
    #     self.__virtual_model = Model("virtual_model")
    #     self.__universal_model = Model("universal_model")
    #
    #     self.__virtual_model_mapper = ModelMapper(self.__virtual_model, self.__compoundsAnnotationConfigs,
    #                                               self.__compoundsIdConverter)
    #
    #     self.__virtual_compounds_revisor = CompoundsRevisor(self.__virtual_model, self.__universal_model,
    #                                                         self.__compoundsAnnotationConfigs)
    #
    #     self.__virtual_compounds_revisor.set_model_mapper(self.__virtual_model_mapper)
    #
    #     granulator = Granulator(self.model, self.mapper, self.__database_format, self.__modelseedCompoundsDb,
    #                             self.__compoundsIdConverter, self.__compoundsAnnotationConfigs)
    #
    #     granulator.set_virtual_model(self.__virtual_model, self.__virtual_model_mapper,
    #                                  self.__virtual_compounds_revisor)
    #
    #     granulator.solve_components_reactions(same_components)


    def __add_new_reactions_to_virtual_model(self,new_reactions):
        self.__virtual_model_mapper.add_new_reactions_to_model(new_reactions)
        self.__virtual_model.add_reactions(new_reactions)

    def __add_new_reactions_to_model(self,new_reactions):
        self.mapper.add_new_reactions_to_model(new_reactions)
        self.model.add_reactions(new_reactions)

    def generateISAreactions(self):

        parents_and_children = self.get_parent_and_children_in_model()

        objective_reactants = [reactant.id for reactant in self.objective.reactants]

        reactions = []
        for parent in parents_and_children:

            if parent in objective_reactants:

                for child in parents_and_children[parent]:

                    name = "ISA_reaction_"+child

                    id = "ISA_reaction_"+child

                    newISAreaction = Reaction(id = id, name = name)

                    child_cobra_container = self.model.metabolites.get_by_id(child)

                    parent_cobra_container = self.model.metabolites.get_by_id(parent)

                    stoichiometry = {child_cobra_container: -1,
                                     parent_cobra_container: 1}

                    newISAreaction.add_metabolites(stoichiometry)

                    reactions.append(newISAreaction)

        self.model.add_reactions(reactions)


    def get_parent_and_children_in_model(self):

        res = {}

        for model_compound_id in self.mapper.boimmg_db_model_map.keys():

            boimmg_id = self.mapper.boimmg_db_model_map[model_compound_id]
            node = self.__compounds_ontology.get_node_by_ont_id(boimmg_id)

            if not node.generic:
                parents = self.__compounds_ontology.get_successors_by_ont_id_rel_type(boimmg_id,"is_a")

                if parents:

                    parent = parents[0]
                    if parent in self.mapper.boimmg_db_model_map_reverse.keys():
                        parents_in_model = self.mapper.boimmg_db_model_map_reverse[parent]

                        if len(parents_in_model) == 1:
                            parent_in_model = parents_in_model[0]

                            if parent_in_model in res:
                                res[parent_in_model].append(model_compound_id)

                            else:
                                res[parent_in_model] = [model_compound_id]
        return res


    def transform_targets_and_components_into_boimmg_ids(self, targets, components):

        """
        Conversion of components and lipid targets into BOIMMG's format

        :param targets:
        :param components:
        :return:
        """

        targets_res = []
        components_res = []
        accessor = CompoundsDBAccessor()
        boimmg_prefix = self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]

        for target in targets:
            new_target = self._transform_database_ids_into_boimmg_ids(target,accessor,boimmg_prefix)

            if not new_target:
                exit(2)

            targets_res.append(new_target)

        for component in components:
            new_component = self._transform_database_ids_into_boimmg_ids(component, accessor, boimmg_prefix)

            if not new_component:
                exit(2)

            components_res.append(new_component)

        return targets_res,components_res


    def _transform_database_ids_into_boimmg_ids(self, target: str, accessor: CompoundsDBAccessor, boimmg_prefix: str):

        """
        Convert target lipid identifier into BOIMMG's format

        :param str target:
        :param CompoundsDBAccessor accessor:
        :param str boimmg_prefix:
        :return:
        """

        targets_res = ''
        if boimmg_prefix in target:
            targets_res = int(target.replace(boimmg_prefix, ""))


        elif re.match("^[0-9]*$", target):

            targets_res = int(target)

        elif "cpd" not in target:
            model_seed_ids = self.__compoundsIdConverter.convert_db_id_to_model_seed_by_db_id(target)

        else:
            model_seed_ids = [target]

        if not targets_res:
            for model_seed_id in model_seed_ids:
                node1 = accessor.get_node_from_model_seed_id(model_seed_id)

                if node1:
                    targets_res = node1.id
                    break

        return targets_res






























