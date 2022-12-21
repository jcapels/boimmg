from copy import deepcopy

from cobra import Reaction
from cobra.util import linear_reaction_coefficients

from src.boimmgpy.database.accessors.compounds_rest_accessor import CompoundsRestAccessor
from src.boimmgpy.service.network_modifiers.granulator import Granulator
from src.boimmgpy.service.interfaces.representation_problem_solver import RepresentationProblemSolver

import logging
import re
import time
from logging.handlers import TimedRotatingFileHandler

import cobra
from biocyc import biocyc
from src.boimmgpy import definitions
from cobra import Model

from src.boimmgpy.database.containers.compound_node import CompoundNode
from src.boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from src.boimmgpy.database.databases_babel import AliasesTransformer
from src.boimmgpy.service.model_mapper import ModelMapper
from src.boimmgpy.service.revisor.compounds_revisor import CompoundsRevisor
from src.boimmgpy.service.network_modifiers.metabolite_swapper import MetaboliteSwapper
from src.boimmgpy.model_seed.model_seed_compound import ModelSeedCompound
from src.boimmgpy.utilities import model_utilities
from src.boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from src.boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDBRest, \
    ModelSeedCompoundsDBRaw
from src.boimmgpy.definitions import TOOL_CONFIG_PATH, ROOT_DIR, COMPOUNDS_ANNOTATION_CONFIGS_PATH, \
    REACTIONS_ANNOTATION_CONFIGS_PATH
from src.boimmgpy.utilities import file_utilities
from src.boimmgpy.utilities.rest_access_utils import RestUtils

logPath = ROOT_DIR + "/logs"

# format the log entries, backup limit and others
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
handler = TimedRotatingFileHandler(logPath, when='midnight', backupCount=20)
handler.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

PROGRESS_BAR = ROOT_DIR + "/service/logs/progress_bar.txt"


# PROGRESS_BAR = "/workdir/resultsWorker/progress_bar.txt"

def write_in_progress_bar(message, value):
    with open(PROGRESS_BAR, "w") as file:
        file.write(message + ":" + str(value))


class LipidGranulator(RepresentationProblemSolver):

    def __init__(self, model, database_format, db_accessor=CompoundsDBAccessor()):
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

        self.__compounds_ontology = db_accessor

        self.model = model
        self.__database_format = database_format
        self.__configs = file_utilities.read_conf_file(TOOL_CONFIG_PATH)
        self.__home_path__ = ROOT_DIR
        self.__define_instance_variables()

        if isinstance(db_accessor, CompoundsDBAccessor):
            self.__modelseedCompoundsDb = ModelSeedCompoundsDBRaw()

        else:
            self.__modelseedCompoundsDb = ModelSeedCompoundsDBRest()

        write_in_progress_bar("mapping model... ", 1)

        write_in_progress_bar("model mapped ", 10)

    def map_model(self):
        """
        Creates maps for all the metabolites in the model

        :return:
        """

        write_in_progress_bar("mapping the model... ", 10)
        print("############# Starting to map the model ################")

        if isinstance(self.__compounds_ontology, CompoundsDBAccessor):
            self.__mapper.map_model(self.__database_format)

        else:
            cobra.io.write_sbml_model(self.__model, definitions.ROOT_DIR + "/temp/model_to_be_submitted.xml")
            response = RestUtils.map_model(definitions.ROOT_DIR + "/temp/model_to_be_submitted.xml",
                                           self.__database_format)
            maps = response.json().get("result")

            self.__mapper.boimmg_db_model_map = maps["boimmg_db_model_map"]
            self.__mapper.boimmg_db_model_map_reverse = maps["boimmg_db_model_map_reverse"]
            self.__mapper.compounds_aliases_indexation = maps["compounds_aliases_indexation"]
            self.__mapper.compounds_aliases_indexation_reverse = maps["compounds_aliases_indexation_reverse"]
            self.__mapper.compound_inchikey_indexation = maps["compound_inchikey_indexation"]

        self.__mapper.mapped = True

        write_in_progress_bar("model mapped ", 31)

        print()
        print("################ model mapped ####################")

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

    def dump_maps(self, destination: str):
        """
        It creates a dump of all the model metabolite maps

        :param str destination: destination folder
        :return:
        """

        self.__mapper.create_map_dump(destination)

    def load_maps(self, folder):
        """
        Load all the metabolite maps in a given folder

        :param str folder: folder path
        :return:
        """

        self.__mapper.upload_maps(folder)

    def __define_instance_variables(self):
        """
        Define instance variables

        :return:
        """

        self.__universal_model = Model("universal_model")
        self.__swapped = []

        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

        self.__compoundsIdConverter = CompoundsIDConverter()

        self.__compounds_revisor = CompoundsRevisor(self.model, self.__universal_model)

        biocyc.set_organism("meta")

        self.mapper = ModelMapper(self.model, self.__compoundsIdConverter, self.__compounds_ontology)

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
        while i < len(metabolites_in_model):

            metabolite = metabolites_in_model[i]

            if metabolite.formula != "" and metabolite.formula in ["H",
                                                                   "H+"] and metabolite.compartment in compartments:
                res.append(metabolite)

            i += 1

        if res:
            self.__virtual_model.add_metabolites(res)

        else:
            raise Exception("Please insert formulas to the hydrogens or protons in your model")

    def swap_and_gap_fill(self, target_ontology_id):
        pass

    def generalize_model(self, target_ontology_id):
        pass

    def gap_fill_model_by_target(self, target_ontology_id: int, components_ont_ids: list):
        pass

    def write_report(self, file_path):

        metabolites_report_material = self.granulator.metabolite_report_material
        reactions_report_material = self.granulator.reaction_report_material

        with open(file_path, "w") as f:

            f.write("--Metabolites--\n")
            f.write("metabolite in model,replacer metabolite\n")

            for metabolite in metabolites_report_material:
                f.write(metabolite + ",")
                f.write(str(metabolites_report_material[metabolite]) + "\n")

            f.write("--Reactions--\n")
            f.write("reaction in model,replacer reaction\n")

            for reaction in reactions_report_material:
                f.write(reaction + ",")
                f.write(str(reactions_report_material[reaction]) + "\n")

    def swap_from_generic(self, targets: list, components: list, same_components=False,
                          progress_bar_processes_left=1, sources=None):
        """
        Swap compounds from a generic target and components. The algorithm acts as follows:

            1st: Starts for identifying the reactions present in the requested biosynthetic pathway;

            2nd: Creates a virtual model;

            3rd: Calls the granulator to granulate the biosynthetic pathway;

        :param list<str> targets: list of lipid targets
        :param list<str> components: list of components
        :param bool same_components: boolean for same or mix of components
        :param int progress_bar_processes_left:
        :param list sources: list with the sources of the structural defined lipids (ModelSEED, LIPID MAPS, SwissLipids)
        :return:
        """

        if sources is None:
            sources = []
        if not self.mapper.mapped:
            self.map_model()

        targets_generic_ontology_id, components_ont_ids = \
            self.transform_targets_and_components_into_boimmg_ids(targets, components)

        self.__virtual_model = Model("virtual_model")

        self.__add_hydrogen_to_virtual_model()

        self.__virtual_model_mapper = ModelMapper(self.__virtual_model,
                                                  self.__compoundsIdConverter, self.__compounds_ontology)

        self.__virtual_compounds_revisor = CompoundsRevisor(self.__virtual_model, self.__universal_model,
                                                            self.__compoundsAnnotationConfigs)

        self.__virtual_compounds_revisor.set_model_mapper(self.__virtual_model_mapper)

        self.granulator.set_virtual_model(self.__virtual_model, self.__virtual_model_mapper,
                                          self.__virtual_compounds_revisor)

        for target_generic_ontology_id in targets_generic_ontology_id:
            biosynthesis_reactions = self.granulator.identify_biosynthesis_pathway(target_generic_ontology_id)

            self.__biosynthesis_reactions = deepcopy(biosynthesis_reactions)

            self.__add_new_reactions_to_virtual_model(self.__biosynthesis_reactions)

            self.__model.remove_reactions(biosynthesis_reactions)

            write_in_progress_bar("starting granulation ", 32)
            self.granulator.granulate(target_generic_ontology_id, components_ont_ids, same_components,
                                      progress_bar_processes_left, sources)

        self.generate_isa_reactions()

    @staticmethod
    def __subtract_reactions_in_list(reactions, reactions_to_subtract):

        res = []

        for reaction in reactions:
            if reaction.id not in reactions_to_subtract:
                res.append(reaction)

        return res

    def __add_new_reactions_to_virtual_model(self, new_reactions):
        self.__virtual_model_mapper.add_new_reactions_to_model(new_reactions)
        self.__virtual_model.add_reactions(new_reactions)

    def __add_new_reactions_to_model(self, new_reactions):
        self.mapper.add_new_reactions_to_model(new_reactions)
        self.model.add_reactions(new_reactions)

    def generate_isa_reactions(self):

        parents_and_children = self.get_parent_and_children_in_model()

        objective_reactants = [reactant.id for reactant in self.objective.reactants]

        reactions = []
        for parent in parents_and_children:

            if parent in objective_reactants:

                for child in parents_and_children[parent]:
                    name = "ISA_reaction_" + child

                    reaction_id = "ISA_reaction_" + child

                    newISAreaction = Reaction(id=reaction_id, name=name)

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
                parents = self.__compounds_ontology.get_successors_by_ont_id_rel_type(boimmg_id, "is_a")

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
            new_target = self._transform_database_ids_into_boimmg_ids(target, accessor, boimmg_prefix)

            if not new_target:
                exit(2)

            targets_res.append(new_target)

        for component in components:
            new_component = self._transform_database_ids_into_boimmg_ids(component, accessor, boimmg_prefix)

            if not new_component:
                exit(2)

            components_res.append(new_component)

        return targets_res, components_res

    def _transform_database_ids_into_boimmg_ids(self, target: str, accessor: CompoundsDBAccessor, boimmg_prefix: str):

        """
        Convert target lipid identifier into BOIMMG's format

        :param str target:
        :param CompoundsDBAccessor accessor:
        :param str boimmg_prefix:
        :return:
        """

        targets_res = ''

        model_seed_ids = []

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


# TODO: implement interface
class CofactorSwapper:

    def __init__(self, model: Model, database_format: str, db_accessor=CompoundsRestAccessor()):
        """
        Class constructor

        :param Model model: model being analysed
        :param database_format: database format (ModelSEED, BiGG or KEGG)
        :param modelseed_compoundsdb:
        """

        self.report_material = {}
        self.__universal_model = Model("universal_model")
        self.model = model
        self.__database_format = database_format
        self.__configs = file_utilities.read_conf_file(TOOL_CONFIG_PATH)
        self.__home_path__ = ROOT_DIR
        self.__compounds_ontology = db_accessor

        if isinstance(db_accessor, CompoundsDBAccessor):
            self.__modelseedCompoundsDb = ModelSeedCompoundsDBRaw()

        else:
            self.__modelseedCompoundsDb = ModelSeedCompoundsDBRest()

        self.__define_instance_variables()

    def dump_maps(self, destination: str):
        """
        It creates a dump of all the model metabolite maps

        :param str destination: destination folder
        :return:
        """

        self.__mapper.create_map_dump(destination)

    def load_maps(self, folder):
        """
        Load all the metabolite maps in a given folder

        :param str folder: folder path
        :return:
        """

        self.__mapper.upload_maps(folder)

    def map_model(self):
        """
        Creates maps for all the metabolites in the model

        :return:
        """

        start = time.time()
        self.write_in_progress_bar("mapping the model... ", 10)
        logger.info("mapping model")

        if isinstance(self.__compounds_ontology, CompoundsDBAccessor):
            self.__mapper.map_model(self.__database_format)

        else:
            cobra.io.write_sbml_model(self.__model, definitions.ROOT_DIR + "/temp/model_to_be_submitted.xml")
            response = RestUtils.map_model(definitions.ROOT_DIR + "/temp/model_to_be_submitted.xml",
                                           self.__database_format)
            maps = response.json().get("result")

            self.__mapper.boimmg_db_model_map = maps["boimmg_db_model_map"]
            self.__mapper.boimmg_db_model_map_reverse = {int(k): maps["boimmg_db_model_map_reverse"][k] for k, v in
                                                         maps["boimmg_db_model_map_reverse"].items()}
            self.__mapper.compounds_aliases_indexation = maps["compounds_aliases_indexation"]
            self.__mapper.compounds_aliases_indexation_reverse = maps["compounds_aliases_indexation_reverse"]
            self.__mapper.compound_inchikey_indexation = maps["compound_inchikey_indexation"]

        self.__mapper.mapped = True

        logger.info("model mapped")

        self.write_in_progress_bar("model mapped ", 31)
        finished = time.time()
        logger.info("mapping finished: %d" % (finished - start))

    @staticmethod
    def write_in_progress_bar(message: str, value: float):

        """
        Write information related to the progress bar

        :param str message: message of what is being done
        :param float value: value to assign in the progress bar
        :return:
        """

        with open(PROGRESS_BAR, "w") as file:
            file.write(message + ":" + str(value))

    @staticmethod
    def get_progress_bar_state() -> float:
        """
        This method reads the progress bar message and the respective value

        :return float state: value to assign in the progress bar
        """
        with open(PROGRESS_BAR, "r") as file:
            line = file.read()
            line = line.strip()
            state = line.split(":")[1]
            return float(state)

    @staticmethod
    def calculate_division_for_progress_bar(state, total_processes_left):
        left = 100 - state
        for_each_part = left / total_processes_left

        return for_each_part

    @property
    def model(self):
        return self.__model

    @model.setter
    def model(self, value):
        if isinstance(value, Model):
            self.__model = value
        else:
            raise ValueError("introduce a cobrapy model")

    def __define_instance_variables(self):
        """
        Method to define a set of instance values

        :return:
        """

        self.__swapped = []
        self.lineage = []

        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)

        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(REACTIONS_ANNOTATION_CONFIGS_PATH)
        self.__compoundsIdConverter = CompoundsIDConverter()

        self.__not_to_change_classes = ["cpd03476"]
        self.__set_not_to_change_compounds()

        self.__compounds_revisor = CompoundsRevisor(self.model, self.__universal_model)
        biocyc.set_organism("meta")

        self.__swapper = MetaboliteSwapper(self.model, None, None, 0,
                                           self.__modelseedCompoundsDb,
                                           self.__compounds_ontology,
                                           model_database=self.__database_format,
                                           compoundsIdConverter=self.__compoundsIdConverter,
                                           universal_model=self.__universal_model,
                                           not_to_change_compounds=self.__not_to_change_compounds
                                           )

        self.__mapper = ModelMapper(self.model,
                                    self.__compoundsIdConverter, self.__compounds_ontology)

        self.__swapper.set_model_mapper(self.__mapper)
        self.__compounds_revisor.set_model_mapper(self.__mapper)

    def __set_not_to_change_compounds(self):
        """
        One would want to protect few compounds of being changed.

        :return:
        """
        self.__not_to_change_compounds = []
        for compound_class in self.__not_to_change_classes:
            parent_id = self.__compounds_ontology.get_node_id_from_model_seed_id(compound_class)
            children = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(parent_id, "is_a")
            self.__not_to_change_compounds.extend(children)

    def swap_compound(self, compound_in_model: str, compound_to_change: str, compounds_left=1):
        """
        This method aims at swapping all the biosynthetic precursors, and the conjugated acid and base of a given
        target. It accepts all types of changes.

        :param str compound_in_model: compound in model
        :param str compound_to_change: compound to replace the one in the model
        :param str compounds_left: compounds left in the swapping operation
        :return:
        """

        if not self.__mapper.mapped:
            self.map_model()

        compound_in_model = self.transform_input_into_boimmg_ids(compound_in_model)
        if not compound_to_change:
            exit(2)

        compound_to_change = self.transform_input_into_boimmg_ids(compound_to_change)
        if not compound_to_change:
            exit(2)

        compound_in_model_node = self.__compounds_ontology.get_node_by_ont_id(compound_in_model)
        compound_to_change_node = self.__compounds_ontology.get_node_by_ont_id(compound_to_change)

        state = self.get_progress_bar_state()
        per_iteration = self.calculate_division_for_progress_bar(state, compounds_left)

        self.write_in_progress_bar("Swapping " + compound_in_model_node.name
                                   + " to " + compound_to_change_node.name, state)

        logger.info("Swapping " + compound_in_model_node.name + " to " + compound_to_change_node.name)

        swapping_type = 0

        if compound_in_model_node.generic and not compound_to_change_node.generic:
            parents = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_to_change, "is_a")
            if compound_in_model in parents:
                swapping_type = 2
            else:
                swapping_type = 3

        elif compound_to_change_node.generic and not compound_in_model_node.generic:
            parents = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_in_model, "is_a")
            if compound_to_change in parents:
                swapping_type = 2
            else:
                swapping_type = 3

        elif not compound_to_change_node.generic and not compound_in_model_node.generic:

            parents1 = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_in_model, "is_a")
            parents2 = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_to_change, "is_a")
            if parents1 and parents2:
                if parents1[0] == parents2[0]:
                    swapping_type = 1
                else:
                    swapping_type = 2

        if swapping_type != 0:
            self.__swapper.set_type(swapping_type)
            self.__swapper.set_compound_in_model(compound_in_model)
            self.__swapper.set_new_compound(compound_to_change)

            self.__swapper.swap_metabolites()

        else:
            self.__swapper.set_type(swapping_type)
            self.__swapper.set_compound_in_model(compound_in_model)
            self.__swapper.set_new_compound(compound_to_change)

            self.__swapper.swap_only_metabolites_and_conjugates()

        self.write_in_progress_bar("Swap performed", per_iteration)
        logger.info("Swap performed")

    def write_report(self, file_path):

        metabolites_report_material = self.__swapper.report_material
        reactions_report_material = self.__swapper.reactions_swapper.report_material

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

    def __check_if_compound_exists_in_model_by_ontology_id(self, ontology_id: int) -> list:
        """
        This method will check whether a BOIMMG compound exists in the model.

        :param int ontology_id: BOIMMG id
        :return list: compounds in model
        """
        res = []

        for model_id, boimmg_id in self.__mapper.boimmg_db_model_map.items():
            if boimmg_id == ontology_id:
                model_compound_container = self.model.metabolites.get_by_id(model_id)
                res.append(model_compound_container)

        return res

    def __replace_electron_donors_and_acceptors(self):
        """
        This method will swap all the generic electron donors and acceptors eventually added when introducing
        MetaCyc reactions.

        :return:
        """

        donor_model_seed_id = "cpd26978"
        electron_transfer_quinol = 749301

        ontology_id = self.__compounds_ontology.get_node_id_from_model_seed_id(donor_model_seed_id)
        container = self.__compounds_ontology.get_node_by_ont_id(ontology_id)
        inchikey = container.inchikey
        aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(donor_model_seed_id)
        metabolites_in_model = self.__mapper.check_metabolites_in_model(inchikey, aliases)
        self.__swapper.set_compound_in_model(ontology_id)
        self.__swapper.set_type(0)

        if metabolites_in_model:

            leaves = self.__compounds_ontology.get_leaves_from_ont_id(ontology_id)

            for metabolite in metabolites_in_model:
                compartment = metabolite.compartment
                found_leaf = False
                i = 0
                while not found_leaf and i < len(leaves):
                    # TODO: test this
                    parents = self.__compounds_ontology.get_all_parents(leaves[i])
                    if electron_transfer_quinol not in parents:

                        leaf_container = self.__compounds_ontology.get_node_by_ont_id(leaves[i])
                        modelseedid = leaf_container.model_seed_id
                        inchikey = leaf_container.inchikey
                        aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(modelseedid)

                        metabolites_in_model_to_change = \
                            self.__mapper.check_metabolites_in_model(inchikey, aliases,
                                                                     leaf_container)

                        found_compartment_leaf = False
                        j = 0
                        while not found_compartment_leaf and j < len(metabolites_in_model_to_change):
                            if metabolites_in_model_to_change[j].compartment == compartment:
                                found_compartment_leaf = True
                                found_leaf = True

                                self.__swapper.set_new_compound(leaves[i])
                                self.__swapper.swap_only_metabolites_and_conjugates(metabolite,
                                                                                    metabolites_in_model_to_change[j])
                            else:
                                j += 1
                    i += 1

    def __swap_metabolites(self, metabolite_in_model: int, new_metabolite: int, swapping_type: int,
                           target=False) -> list:
        """
        This method aims at calling the MetaboliteSwapper to swap a given metabolite, their biosynthetic precursors
        and their conjugated base and acid, eventually.

        :param int metabolite_in_model: BOIMMG id of the metabolite in model
        :param int new_metabolite: BOIMMG id of the new metabolite
        :param int swapping_type: type of swap
        :param boolean target: whether the new metabolite is the target metabolite or not
        :return:
        """

        self.__swapper.set_compound_in_model(metabolite_in_model)
        self.__swapper.set_new_compound(new_metabolite)
        self.__swapper.set_type(swapping_type)

        if swapping_type == 2:
            self.__type_2_is_generic_setter(metabolite_in_model)

        if not target:
            changed_reactions = self.__swapper.swap_only_metabolites_and_conjugates()
            self.__swapped.append(new_metabolite)
        else:
            changed_reactions = self.__swapper.swap_metabolites()

        self.__universal_model = self.__swapper.get_universal_model()

        return changed_reactions

    def __type_2_is_generic_setter(self, metabolite_in_model):

        container = self.__compounds_ontology.get_node_by_ont_id(metabolite_in_model)
        if container.generic:
            self.__swapper.setType2_isGenericInModel(True)
        else:
            self.__swapper.setType2_isGenericInModel(False)

    def __check_and_change_all_children(self, parent: int, new_metabolite: int, target=False) -> list:
        """
        This method absorbs all the compounds of a specific chemical type and swap them into the :param new_metabolite

        :param int parent: BOIMMG id of the structural parent
        :param int new_metabolite: BOIMMG id of the metabolite that will prevail
        :param boolean target: whether the new metabolite is the target or not
        :return list: changed reactions
        """

        children = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(parent, "is_a")
        changed_reactions = []

        for child in children:
            if child != new_metabolite:

                child_container = self.__compounds_ontology.get_node_by_ont_id(child)

                aliases = AliasesTransformer.transform_boimmg_aliases_into_model_seed(child_container.aliases)

                model_seed_id = child_container.model_seed_id

                if model_seed_id:

                    converted_aliases = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(model_seed_id)

                    merged_aliases = self.merge_dictionaries(converted_aliases, aliases)

                    child_in_model = self.__mapper.check_metabolites_in_model(child_container.inchikey,
                                                                              merged_aliases)

                else:

                    child_in_model = self.__mapper.check_metabolites_in_model(child_container.inchikey,
                                                                              child_container.aliases,
                                                                              child_container)

                if child_in_model:
                    changed_reactions = self.__swap_metabolites(child, new_metabolite, 1, target)

        return changed_reactions

    def convert_aliases_into_mapper_format(self, aliases):
        """
        Such method can convert the aliases from retrieved from ModelSEED into the mapper format

        :param dict aliases: aliases from ModelSEED
        :return dict new_aliases: the converted dictionary
        """
        new_aliases = {}
        for alias in aliases:
            if alias in self.__compoundsAnnotationConfigs.keys():
                new_key = self.__compoundsAnnotationConfigs[alias]
                new_aliases[new_key] = aliases[alias]

        return new_aliases

    @staticmethod
    def merge_dictionaries(dict1, dict2):
        for key in dict1:
            if key in dict2:
                dict1[key].extend(dict2[key])
        return dict1

    def generate_new_metabolite(self, compound_container):
        """
        Such method generates a whole new metabolite using a ModelSEED compound wrapper or a CompoundNode

        :param compound_container: ModelSEED compound wrapper
        :return list: list of newly generated metabolites
        """

        if isinstance(compound_container, ModelSeedCompound):
            model_metabolites = \
                model_utilities.generate_model_compounds_by_database_format(self.model,
                                                                            compound_container.getDbId(),
                                                                            self.__compoundsIdConverter,
                                                                            self.__modelseedCompoundsDb,
                                                                            self.__database_format)
        elif isinstance(compound_container, CompoundNode):
            model_metabolites = \
                model_utilities.generate_model_compounds_by_database_format(self.model,
                                                                            compound_container.model_seed_id,
                                                                            self.__compoundsIdConverter,
                                                                            self.__modelseedCompoundsDb,
                                                                            self.__database_format)

        else:
            raise ValueError("Not expected type")

        self.__model.add_metabolites(model_metabolites)
        self.__mapper.add_new_metabolites_to_maps(model_metabolites)
        return model_metabolites

    def generate_new_boimmg_metabolite(self, boimmg_container: CompoundNode) -> list:

        """
        Such method generates new BOIMMG compounds and adds them to the model (same compound in different compartments)

        :param CompoundNode boimmg_container: BOIMMG container for model addition
        :return list<Metabolite>: list of model metabolites in different compartments
        """

        model_metabolites = model_utilities.generate_boimmg_metabolites(self.model,
                                                                        boimmg_container,
                                                                        self.__database_format,
                                                                        self.__compoundsIdConverter,
                                                                        self.__compoundsAnnotationConfigs,
                                                                        )

        self.__model.add_metabolites(model_metabolites)
        self.__mapper.add_new_metabolites_to_maps(model_metabolites)
        return model_metabolites

    def add_metabolites_to_model(self, metabolites):

        for metabolite in metabolites:

            container = self.__compounds_ontology.get_node_by_ont_id(metabolite)
            metabolites_in_model = self.__mapper.check_metabolites_in_model(container.inchikey, container.aliases,
                                                                            container)

            if not metabolites_in_model:

                if container.model_seed_id:
                    ms_container = self.__modelseedCompoundsDb.get_compound_by_id(container.model_seed_id)
                    self.generate_new_metabolite(ms_container)

                else:
                    self.generate_new_boimmg_metabolite(container)

    def transform_input_into_boimmg_ids(self, target: str) -> int:
        """
        Method that converts a target ID from an external database into BOIMMG's

        :param str target: target identifier of an external database
        :return int: boimmg identifier
        """

        accessor = CompoundsDBAccessor()

        boimmg_prefix = self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]
        target_res = self._transform_database_ids_into_boimmg_ids(target, accessor, boimmg_prefix)

        return target_res

    def _transform_database_ids_into_boimmg_ids(self, target: str, accessor: CompoundsDBAccessor,
                                                boimmg_prefix: str) -> int:

        """
        Such method transforms a target ID from an external database into BOIMMG's

        :param str target: target identifier from an external database
        :param CompoundsDBAccessor accessor:
        :param str boimmg_prefix: BOIMMG prefix
        :return:
        """

        targets_res = ''
        model_seed_ids = []
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
