import json
import logging
import re
import time
from logging.handlers import TimedRotatingFileHandler

import cobra
from biocyc import biocyc
from boimmgpy import definitions
from cobra import Model, Metabolite, Reaction
from cobra.flux_analysis.gapfilling import GapFiller

from boimmgpy.database.accessors.compounds_rest_accessor import CompoundsRestAccessor
from boimmgpy.database.containers.compound_node import CompoundNode
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.database.databases_babel import AliasesTransformer
from boimmgpy.service.model_mapper import ModelMapper
from boimmgpy.service.network_handlers.pathway_handler import PathwayHandler
from boimmgpy.service.revisor.compounds_revisor import CompoundsRevisor
from boimmgpy.service.network_modifiers.metabolite_swapper import MetaboliteSwapper
from boimmgpy.service.interfaces.representation_problem_solver import RepresentationProblemSolver
from boimmgpy.model_seed.model_seed_compound import ModelSeedCompound
from boimmgpy.utilities import model_utilities
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB, ModelSeedCompoundsDBRest, \
    ModelSeedCompoundsDBRaw
from boimmgpy.definitions import TOOL_CONFIG_PATH, ROOT_DIR, COMPOUNDS_ANNOTATION_CONFIGS_PATH, REACTIONS_ANNOTATION_CONFIGS_PATH
from boimmgpy.utilities import file_utilities
from boimmgpy.utilities.rest_access_utils import RestUtils

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

class SimpleCaseSolver():

    def __init__(self, model: Model, database_format: str, db_accessor = CompoundsDBAccessor()):
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


        if isinstance(db_accessor,CompoundsDBAccessor):
            self.__modelseedCompoundsDb = ModelSeedCompoundsDBRaw()
            self.__define_instance_variables()
            self.__mapper.map_model(database_format)

        else:
            self.__modelseedCompoundsDb = ModelSeedCompoundsDBRest()
            self.__define_instance_variables()

            cobra.io.write_sbml_model(model,definitions.ROOT_DIR + "/temp/model_to_be_submitted.xml")
            response = RestUtils.map_model(definitions.ROOT_DIR + "/temp/model_to_be_submitted.xml", database_format)
            maps = json.load(response)

            self.__mapper.boimmg_db_model_map = maps["boimmg_db_model_map"]
            self.__mapper.boimmg_db_model_map_reverse = maps["boimmg_db_model_map_reverse"]
            self.__mapper.compounds_aliases_indexation = maps["compounds_aliases_indexation"]
            self.__mapper.compounds_aliases_indexation_reverse = maps["compounds_aliases_indexation_reverse"]
            self.__mapper.compound_inchikey_indexation = maps["compound_inchikey_indexation"]



    def dump_maps(self,destination:str):
        """
        It creates a dump of all the model metabolite maps

        :param str destination: destination folder
        :return:
        """

        self.__mapper.create_map_dump(destination)

    def load_maps(self,folder):
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
        self.__mapper.map_model(self.__database_format)
        logger.info("model mapped")

        self.write_in_progress_bar("model mapped ", 31)
        finished = time.time()
        logger.info("mapping finished: %d" % (finished - start))

    def write_in_progress_bar(self, message: str, value: float):

        """
        Write information related to the progress bar

        :param str message: message of what is being done
        :param float value: value to assign in the progress bar
        :return:
        """

        with open(PROGRESS_BAR,"w") as file:
            file.write(message + ":" + str(value))

    def get_progress_bar_state(self) -> float:
        """
        This method reads the progress bar message and the respective value

        :return float state: value to assign in the progress bar
        """
        with open(PROGRESS_BAR,"r") as file:
            line = file.read()
            line = line.strip()
            state = line.split(":")[1]
            return float(state)

    def calculate_division_for_progress_bar(self,state,total_processes_left):
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

        # self.__compounds_ontology = CompoundsDBAccessor()

        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)

        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(REACTIONS_ANNOTATION_CONFIGS_PATH)
        self.__compoundsIdConverter = CompoundsIDConverter()

        # self.__ptw_handler_quinones = PathwayHandler(self.__modelseedCompoundsDb, self.__compoundsIdConverter)
        # self.__ptw_handler_quinones.load_from_file(self.__home_path__ + self.__configs["quinones_pathways"])

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
                                    self.__compoundsIdConverter,self.__compounds_ontology)

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


    def swap_compound(self, compound_in_model:str, compound_to_change: str, compounds_left=1):
        """
        This method aims at swapping all the biosynthetic precursors, and the conjugated acid and base of a given target.
        It accepts all types of changes.

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
        per_iteration = self.calculate_division_for_progress_bar(state,compounds_left)


        self.write_in_progress_bar("Swapping "+compound_in_model_node.name
                                   + " to " + compound_to_change_node.name,state)

        logger.info("Swapping "+compound_in_model_node.name + " to " + compound_to_change_node.name)

        type = 0

        if compound_in_model_node.generic and not compound_to_change_node.generic:
            parents = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_to_change, "is_a")
            if compound_in_model in parents:
                type = 2
            else:
                type = 3

        elif compound_to_change_node.generic and not compound_in_model_node.generic:
            parents = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_in_model, "is_a")
            if compound_to_change in parents:
                type = 2
            else:
                type = 3

        elif not compound_to_change_node.generic and not compound_in_model_node.generic:

            parents1 = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_in_model, "is_a")
            parents2 = self.__compounds_ontology.get_successors_by_ont_id_rel_type(compound_to_change, "is_a")
            if parents1 and parents2:
                if parents1[0] == parents2[0]:
                    type = 1
                else:
                    type = 2

        if type != 0:
            self.__swapper.set_type(type)
            self.__swapper.set_compound_in_model(compound_in_model)
            self.__swapper.set_new_compound(compound_to_change)

            self.__swapper.swap_metabolites()

        else:
            self.__swapper.set_type(type)
            self.__swapper.set_compound_in_model(compound_in_model)
            self.__swapper.set_new_compound(compound_to_change)

            self.__swapper.swap_only_metabolites_and_conjugates()

        self.write_in_progress_bar("Swap performed", per_iteration)
        logger.info("Swap performed")


    def write_report(self,file_path):

        metabolites_report_material = self.__swapper.report_material
        reactions_report_material = self.__swapper.reactions_swapper.report_material

        with open(file_path,"w") as f:

            f.write("--Metabolites--\n")
            f.write("metabolite in model,replacer metabolite\n")

            for metabolite in metabolites_report_material:
                f.write(metabolite+",")
                f.write(metabolites_report_material[metabolite]+"\n")

            f.write("--Reactions--\n")
            f.write("reaction in model,replacer reaction\n")

            for reaction in reactions_report_material:
                f.write(reaction+",")
                f.write(reactions_report_material[reaction]+"\n")


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



    def __swap_metabolites(self, metabolite_in_model: int, new_metabolite: int, type: int, target=False) -> list:
        """
        This method aims at calling the MetaboliteSwapper to swap a given metabolite, their biosynthetic precursors
        and their conjugated base and acid, eventually.

        :param int metabolite_in_model: BOIMMG id of the metabolite in model
        :param int new_metabolite: BOIMMG id of the new metabolite
        :param int type: type of swap
        :param boolean target: whether the new metabolite is the target metabolite or not
        :return:
        """

        self.__swapper.set_compound_in_model(metabolite_in_model)
        self.__swapper.set_new_compound(new_metabolite)
        self.__swapper.set_type(type)

        if type == 2:
            self.__type_2_isGeneric_setter(metabolite_in_model)

        if not target:
            changed_reactions = self.__swapper.swap_only_metabolites_and_conjugates()
            self.__swapped.append(new_metabolite)
        else:
            changed_reactions = self.__swapper.swap_metabolites()

        self.__universal_model = self.__swapper.get_universal_model()

        return changed_reactions

    def __type_2_isGeneric_setter(self, metabolite_in_model):

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

    def convert_aliases_into_mapper_format(self,aliases):
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
            model_metabolites = model_utilities.generate_model_compounds_by_database_format(self.model,
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

    def transform_input_into_boimmg_ids(self, target: str) -> str:
        """
        Method that converts a target ID from an external database into BOIMMG's

        :param str target: target identifier of an external database
        :return str: boimmg identifier
        """

        accessor = CompoundsDBAccessor()

        boimmg_prefix = self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]
        target_res = self._transform_database_ids_into_boimmg_ids(target,accessor,boimmg_prefix)


        return target_res

    def _transform_database_ids_into_boimmg_ids(self, target: str, accessor: CompoundsDBAccessor, boimmg_prefix: str)-> str:

        """
        Such method transforms a target ID from an external database into BOIMMG's

        :param str target: target identifier from an external database
        :param CompoundsDBAccessor accessor:
        :param str boimmg_prefix: BOIMMG prefix
        :return:
        """

        targets_res=''
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