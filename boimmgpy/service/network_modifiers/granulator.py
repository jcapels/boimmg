import re
import time
from copy import deepcopy
from typing import Tuple, List

from cobra import Reaction, Model, Metabolite

from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from boimmgpy.service.model_mapper import ModelMapper
from boimmgpy.utilities import chemo_utilities, model_utilities, file_utilities
from boimmgpy.definitions import ROOT_DIR, REACTIONS_ANNOTATION_CONFIGS_PATH

PROGRESS_BAR = ROOT_DIR + "/service/logs/progress_bar.txt"
# PROGRESS_BAR = "/workdir/resultsWorker/progress_bar.txt"

def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='|', printEnd="\r"):
    """
  Call in a loop to create terminal progress bar
  @params:
      iteration   - Required  : current iteration (Int)
      total       - Required  : total iterations (Int)
      prefix      - Optional  : prefix string (Str)
      suffix      - Optional  : suffix string (Str)
      decimals    - Optional  : positive number of decimals in percent complete (Int)
      length      - Optional  : character length of bar (Int)
      fill        - Optional  : bar fill character (Str)
      printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
  """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="",flush=True)
    # Print New Line on Complete
    if iteration == total:
        print()

class Granulator:

    def __init__(self, model: Model, mapper: ModelMapper, database_format: str, id_converter: CompoundsIDConverter,
                 compoundsAnnotationConfigs: dict):

        """
        Class constructor

        :param Model model: Cobrapy model
        :param ModelMapper mapper: model mapper
        :param str database_format: ModelSEED, KEGG or BiGG
        :param CompoundsIDConverter id_converter:
        :param dict compoundsAnnotationConfigs: compounds annotation configurations
        """

        self.__compoundsAnnotationConfigs = compoundsAnnotationConfigs
        self.__model = model
        self.__compoundsIdConverter = id_converter
        self.__modelseedCompoundsDb = ModelSeedCompoundsDB()
        self.__database_format = database_format
        self.mapper = mapper

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
    def virtual_model_mapper(self):
        return self.__virtual_model_mapper

    @virtual_model_mapper.setter
    def virtual_model_mapper(self, value):
        self.__virtual_model_mapper = value

    def define_instance_variables(self):
        """
        Instance variables definition

        :return:
        """

        self.__compounds_ontology  = CompoundsDBAccessor()

        self.__added_reactions_and_respective_target = {}

        self.__reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

        self.metabolite_report_material = {}
        self.reaction_report_material = {}
        self.__set_components()

    def set_virtual_model(self,virtual_model, virtual_model_mapper, virtual_compounds_revisor):

        self.__virtual_model = virtual_model
        self.__virtual_compounds_revisor = virtual_compounds_revisor
        self.__virtual_model_mapper = virtual_model_mapper

    def write_in_progress_bar(self, message, value):
        with open(PROGRESS_BAR,"w") as file:
            file.write(message + ":" + str(value))

    def get_progress_bar_state(self):
        with open(PROGRESS_BAR,"r") as file:
            line = file.read()
            line = line.strip()
            state = line.split(":")[1]
            return float(state)

    def calculate_division_for_progress_bar(self,state,total,total_processes_left):
        left = 100 - state
        for_each_part = left / total_processes_left
        for_each_iteration = for_each_part / total

        return for_each_iteration

    def granulate(self, target_generic_ontology_id, components, same_components, progress_bar_processes_left):
        """
        Operation to granulate a given compound:

            1st: Build up structurally defined compounds using the information on the components list. It builds up
            compounds by mixing :param components if :param same_components is False, otherwise it builds up compounds
            without a mix of :param components.

            2nd: Calls a function to generate the whole biosynthetic pathway of all the built up compounds.

        :param int target_generic_ontology_id: generic compound database identifier
        :param list<int> components: list of BOIMMG identifiers of all the requested components
        :param bool same_components:
        :param int progress_bar_processes_left:
        :return:
        """


        if same_components:
            targets_to_replace = \
                self.__compounds_ontology.get_compounds_with_only_one_component(target_generic_ontology_id, components)

        else:
            targets_to_replace = \
                self.__compounds_ontology.get_compounds_with_specific_parent_within_set_of_components(
                    target_generic_ontology_id,
                    components)

        targets_to_replace = self.filter_targets_to_replace(targets_to_replace,target_generic_ontology_id)

        state = self.get_progress_bar_state()
        for_each_iteration = self.calculate_division_for_progress_bar(
                                        state,len(targets_to_replace),progress_bar_processes_left)

        for target in targets_to_replace:

            self.write_in_progress_bar("generating network for BMGC"+ str(target),state)
            self.handle_target_in_virtual_model(target, target_generic_ontology_id)
            state += for_each_iteration


    def handle_target_in_virtual_model(self, target: int, parent: int):
        """
        Method to granulate a given lipid target in the virtual model.

            1st: Processes the target biosynthesis pathway;
            2nd: Balance the granulated reactions in the previously processed biosynthesis pathway;
            3rd: Adds them to the model

        :param int target: lipid target BOIMMG identifier (structurally defined)
        :param int parent: lipid parent BOIMMG identifier (generic lipid)
        :return:
        """

        granulated_reactions = self.__process_biosynthesis_pathway(target, parent)


        granulated_reactions = self.__virtual_compounds_revisor.balance_reactions(
            granulated_reactions)

        reactions_to_add = deepcopy(granulated_reactions)

        self.__add_new_reactions_to_model(reactions_to_add)


    def identify_biosynthesis_pathway(self, target_generic_ontology_id: int) -> list:
        """
        Identifies the biosynthesis pathway of a given generic lipid.

        :param int target_generic_ontology_id:
        :return list<Reaction>: list of reactions
        """

        reactions_to_delete = self.get_target_associated_network(target_generic_ontology_id)

        return reactions_to_delete


    def __process_biosynthesis_pathway(self, target: int, parent: int) -> list:
        """
        Granulates the biosynthesis pathway of a given target lipid. It uses a "processing" stack to process each reaction
        and each new structurally defined lipid. When the stack is empty it is assumed the pathway to be fully processed.

        :param int target: lipid target BOIMMG identifier (structurally defined)
        :param int parent: lipid parent BOIMMG identifier (generic lipid)
        :return list<Reaction>: list of granulated reactions
        """

        processed_compounds = []
        new_compounds = [(target, parent)]
        granulated_reactions = []

        while new_compounds:
            new_compound_and_parent = new_compounds.pop()

            if new_compound_and_parent[0] not in processed_compounds:


                new_reactions, added_compounds = self.__process_new_compound(new_compound_and_parent)

                processed_compounds.append(new_compound_and_parent[0])

                self.__add_new_reactions_to_virtual_model(new_reactions)

                granulated_reactions.extend(new_reactions)

                new_compounds.extend(added_compounds)

        return granulated_reactions


    def __process_new_compound(self, new_compound_and_parent: tuple):
        """
        Processes a new compound

        :param tuple new_compound_and_parent: (structurally defined compound, its generic parent)
        :return:
        """

        new_compound = new_compound_and_parent[0]

        parent = new_compound_and_parent[1]
        model_boimmg_parents = self.virtual_model_mapper.check_if_boimmg_metabolite_in_model(parent)

        if not model_boimmg_parents:
            exit(2)
            raise Exception("please introduce a mapped biosynthetic target")

        parent_and_precursors = self.get_targets_to_replace_by_type(new_compound)

        new_reactions = []

        added_compounds = []

        processed = []
        for model_parent in model_boimmg_parents:

            if model_parent not in processed:
                cobra_model_parent = self.__virtual_model.metabolites.get_by_id(model_parent)
                reactions_to_process = self.get_reactions_by_compound(cobra_model_parent)

                for reaction_to_process in reactions_to_process:

                    if reaction_to_process.id not in self.__added_reactions_and_respective_target.keys() or \
                            new_compound not in self.__added_reactions_and_respective_target[reaction_to_process.id]:

                        go = True
                        compounds_ids = [comp.id for comp in reaction_to_process.products]

                        if go and model_parent in compounds_ids:
                            new_reaction, added = self._process_reaction(reaction_to_process, new_compound,
                                                                         cobra_model_parent, parent_and_precursors)


                            if reaction_to_process.id in self.reaction_report_material:
                                self.reaction_report_material[reaction_to_process.id].append(new_reaction.id)
                            else:
                                self.reaction_report_material[reaction_to_process.id] = [new_reaction.id]

                            new_reactions.append(new_reaction)

                            for compound in added:
                                if compound not in added_compounds:
                                    added_compounds.append(compound)

                            if reaction_to_process.id not in self.__added_reactions_and_respective_target.keys():
                                self.__added_reactions_and_respective_target[reaction_to_process.id] = [new_compound]

                            else:
                                self.__added_reactions_and_respective_target[reaction_to_process.id].append(new_compound)

                processed.append(model_parent)

        return new_reactions, added_compounds


    def get_reactions_by_compound(self, compound):


        res = []
        for reaction in self.__virtual_model.reactions:
            metabolites_ids = [m.id for m in reaction.metabolites]

            if compound.id in metabolites_ids:
                res.append(reaction)

        return res


    def _process_reaction(self, reaction_to_process: Reaction, new_compound: int, model_parent: Metabolite,
                          parent_and_precursors: dict) -> Tuple[Reaction, List[Tuple[Metabolite, int]]]:
        """
        Granulates the :param reaction_to_process

        :param Reaction reaction_to_process: reaction in model to be granulated
        :param int new_compound: structurally defined lipid BOIMMG identifier
        :param Metabolite model_parent: model generic lipid
        :param dict parent_and_precursors: {parent: [children]} -> children of the generic precursors
        :return tuple: the granulated reaction in the cobrapy Reaction object and list of tuples: (added compound, BOIMMG ID)
        """

        new_reaction = deepcopy(reaction_to_process)

        new_reaction_products = list(new_reaction.products)

        new_reaction, added = self.__process_reaction(new_reaction, new_compound, model_parent, parent_and_precursors,
                                                          new_reaction_products)

        reaction_id = self.change_boimmg_reaction_id(new_reaction)

        new_reaction.id = reaction_id

        return new_reaction, added

    def change_boimmg_reaction_id(self, reaction: Reaction) -> str:

        """
        Changes the reaction ID into a cannonical BOIMMG reaction identifier

        :param Reaction reaction:
        :return str: new identifier
        """

        new_id = self.__reactionsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]

        metabolites = reaction.metabolites

        metabolites = sorted([metabolite.id for metabolite in metabolites])
        new_id += "_".join(metabolites)

        if len(new_id) > 250:
            new_id = new_id[:230]

        reaction.id = new_id

        return reaction.id



    def __process_reaction(self, new_reaction: Reaction, new_compound: int, model_parent: Metabolite,
                           parent_and_precursors: dict,
                           new_reaction_metabolites: List[Metabolite]) -> Tuple[Reaction, List[Tuple[Metabolite, int]]]:

        """
        Granulates a given reaction. Replaces the generic metabolite by its structurally defined child and corrects the
        precursors.

        :param Reaction new_reaction: the new reaction
        :param int new_compound: the reaction structurally defined product
        :param Metabolite model_parent: the model parent of the structurally defined product
        :param dict parent_and_precursors: dictionary with lipid parents and their children
        :param List[Metabolite] new_reaction_metabolites: list of possibly new metabolites
        :return Tuple[Reaction, List[Tuple[Metabolite, int]]]: returns a tuple with the granulated model reaction
        and the new added compounds that will be processed further.
        """

        i = 0
        found = False

        while not found and i < len(new_reaction_metabolites):
            if new_reaction_metabolites[i].id == model_parent.id:
                new_model_parent = new_reaction_metabolites[i]
                found = True

            i += 1

        coef = new_reaction.get_coefficient(model_parent.id)

        new_metabolite = self.add_boimmg_metabolites_to_reaction(new_compound, coef, new_reaction, model_parent.compartment)

        if model_parent.id in self.metabolite_report_material:
            self.metabolite_report_material[model_parent.id].append(new_metabolite.id)
        else:
            self.metabolite_report_material[model_parent.id] = [new_metabolite.id]

        met_to_subtract = {}

        met_to_subtract[new_model_parent] = coef
        new_reaction.subtract_metabolites(met_to_subtract)

        added = self._correct_precursors(new_reaction, new_compound, parent_and_precursors)

        return new_reaction, added


    def _correct_precursors(self, new_reaction: Reaction, new_compound:int, parent_and_precursors:dict) -> \
            List[Tuple[Metabolite, int]]:

        """
        Correct the precursors in the new reaction.

        :param Reaction new_reaction: the new reaction
        :param int new_compound: the new structurally defined lipid species BOIMMG identifier
        :param dict parent_and_precursors: dictionary with lipid parents and their children
        :return List[Tuple[Metabolite, int]]: returns a tuple with the granulated model reaction
        and the new added compounds that will be processed further.
        """

        new_compound_components = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(new_compound,
                                                                                                "component_of")

        reactants = list(new_reaction.reactants)
        added = self.__correct_precursors(new_reaction, parent_and_precursors, reactants,
                                              new_compound_components)

        return added


    def __correct_precursors(self, new_reaction: Reaction, parent_and_precursors: dict,
                             metabolites: List[Metabolite], new_compound_components: List[int]) -> \
                            List[Tuple[Metabolite, int]]:

        """
        Correct the precursors in the new reaction.

        :param Reaction new_reaction: the new reaction
        :param dict parent_and_precursors: dictionary with lipid parents and their children
        :param List[Metabolite] metabolites: number of reactants in the new reaction to be granulated
        :param List[int] new_compound_components: components of the new product in the reaction
        :return:
        """

        added = []
        while metabolites:

            reactant = metabolites.pop()

            boimmg_id = self.virtual_model_mapper.get_boimmg_id_from_model_compound_id(reactant.id)
            if boimmg_id:

                coef = new_reaction.metabolites[reactant]
                if boimmg_id in parent_and_precursors:

                    precursors = parent_and_precursors[boimmg_id]
                    for precursor in precursors:
                        self.add_boimmg_metabolites_to_reaction(precursor, coef, new_reaction, reactant.compartment)

                        added.append((precursor, boimmg_id))

                    met_to_subtract = {}

                    met_to_subtract[reactant] = coef
                    new_reaction.subtract_metabolites(met_to_subtract)



                else:
                    to_replace = self.__compounds_ontology.get_compounds_with_specific_parent_set_of_components(
                        boimmg_id,
                        new_compound_components)

                    if to_replace:
                        to_replace = self.filter_targets_to_replace(to_replace,boimmg_id)

                        self.add_boimmg_metabolites_to_reaction(to_replace[0], coef, new_reaction, reactant.compartment)

                        met_to_subtract = {}

                        met_to_subtract[reactant] = coef
                        new_reaction.subtract_metabolites(met_to_subtract)

                        added.append((to_replace[0], boimmg_id))
        return added


    def add_boimmg_metabolites_to_reaction(self, boimmg_id: int, coef: int, reaction: Reaction, compartment: str):
        """
        Add a metabolite to a given reaction receiving BOIMMG identifier as input.

        :param int boimmg_id: BOIMMG identifier of the metabolite being inserted in the reaction
        :param int coef: stoichimetric coefficient in the reaction
        :param Reaction reaction: model reaction
        :param str compartment: compartment
        :return:
        """

        model_seed_ids = self.__compounds_ontology.get_all_model_seed_ids(boimmg_id)
        node = self.__compounds_ontology.get_node_by_ont_id(boimmg_id)

        aliases = node.aliases

        found = False
        i = 0
        while not found and i<len(model_seed_ids):

            aliases2 = self.__compoundsIdConverter.get_all_aliases_by_modelSeedID(model_seed_ids[i])

            aliases.update(aliases2)

            model_child = self.mapper.check_if_boimmg_metabolite_in_model(boimmg_id,aliases)

            if model_child:
                found = True

            i+=1

        if not model_seed_ids:
            model_child = self.mapper.check_if_boimmg_metabolite_in_model(boimmg_id)


        if not model_child:

            new_model_compound = self.__create_boimmg_metabolite(boimmg_id, compartment)

            met_to_add = {}

            met_to_add[new_model_compound] = coef
            reaction.add_metabolites(met_to_add)


        else:
            found = False
            i = 0
            while not found and i < len(model_child):
                cobra_child_container = self.model.metabolites.get_by_id(model_child[i])
                if cobra_child_container.compartment == compartment:
                    found = True

                    if self.__virtual_model.metabolites.has_id(model_child[i]):
                        found_model_compound = cobra_child_container

                    else:
                        virtual_compound = deepcopy(cobra_child_container)
                        self.__virtual_model.add_metabolites([virtual_compound])
                        found_model_compound = virtual_compound

                i += 1

            if not found:
                new_model_compound = self.__create_boimmg_metabolite(boimmg_id, compartment)

                met_to_add = {}

                met_to_add[new_model_compound] = coef
                reaction.add_metabolites(met_to_add)

            else:
                met_to_add = {}

                new_model_compound = found_model_compound
                met_to_add[found_model_compound] = coef
                reaction.add_metabolites(met_to_add)

        return new_model_compound

    def __create_boimmg_metabolite(self, boimmg_id:int, compartment:str) -> Metabolite:

        """
        Generates a new model compound using the internal database information.

        :param int boimmg_id: BOIMMG identifier
        :param str compartment: compartment
        :return Metabolite: the new metabolite
        """

        new_compound_node = self.__compounds_ontology.get_node_by_ont_id(boimmg_id)

        model_seed_ids = self.__compounds_ontology.get_all_model_seed_ids(boimmg_id)

        for model_seed_id in model_seed_ids:

            translation = self.__compoundsIdConverter.get_modelSeedIdToDb().get(model_seed_id)

            definitive_model_seed_id = model_seed_id
            if self.__database_format == "BiGG" and \
                    ("BiGG" in translation.keys() or "BiGG1" in translation.keys()):

                definitive_model_seed_id = model_seed_id
                break


            elif self.__database_format == "KEGG" and \
                    "KEGG" in translation.keys():

                definitive_model_seed_id = model_seed_id
                break

        new_model_compound = model_utilities.generate_boimmg_metabolites_in_compartment(new_compound_node, compartment,
                                                                                       self.__compoundsIdConverter,
                                                                                       self.__compoundsAnnotationConfigs)
        self.add_metabolites_to_the_virtual_model([new_model_compound])

        return new_model_compound


    def add_metabolites_to_the_virtual_model(self, metabolites):
        self.__virtual_model.add_metabolites(metabolites)
        self.virtual_model_mapper.add_new_metabolites_to_maps(metabolites)


    @staticmethod
    def merge_dictionaries(dict1, dict2):
        for key in dict1:
            if key in dict2:
                dict1[key].extend(dict2[key])
        return dict1


    def get_target_associated_network(self, target_generic_ontology_id: int) -> List[Reaction]:
        """
        Identifies the target's biosynthesis network in the model

        :param int target_generic_ontology_id:
        :return List[Reaction]: list of reactions in the network
        """

        res = []
        res_ids = []

        predecessors = \
            self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(target_generic_ontology_id, "precursor_of")

        predecessors_and_target = predecessors.copy()
        predecessors_and_target.append(target_generic_ontology_id)

        network_metabolites_in_model = []

        for token in predecessors_and_target:
            metabolites_in_model = self.__check_if_compound_exists_in_model_by_ontology_id(token)
            network_metabolites_in_model.extend(metabolites_in_model)

        network_metabolites_ids = [m.id for m in network_metabolites_in_model]

        for metabolites_in_model in network_metabolites_in_model:
            reactions = metabolites_in_model.reactions
            for reaction in reactions:

                if reaction.id not in res_ids:

                    reactants = reaction.reactants
                    products = reaction.products

                    reactant_found = False
                    product_found = False
                    i = 0
                    while not reactant_found and i < len(reactants):
                        reactant = reactants[i]
                        reactant_id = reactant.id
                        if reactant_id in network_metabolites_ids:
                            reactant_found = True

                        i += 1

                    j = 0
                    while not product_found and j < len(products):
                        product = products[j]
                        product_id = product.id
                        if product_id in network_metabolites_ids:
                            product_found = True
                        j += 1

                    if product_found and reactant_found:
                        res.append(reaction)
                        res_ids.append(reaction.id)


        return res


    def __check_if_compound_exists_in_model_by_ontology_id(self, ontology_id: int) -> List[Metabolite]:
        """
        This method will check whether a BOIMMG compound exists in the model.

        :param int ontology_id: BOIMMG id
        :return List[Metabolite]: compounds in model
        """
        res = []

        for model_id, boimmg_id in self.mapper.boimmg_db_model_map.items():
            if boimmg_id == ontology_id:
                model_compound_container = self.model.metabolites.get_by_id(model_id)
                res.append(model_compound_container)

        return res

    def __add_new_reactions_to_virtual_model(self,new_reactions):
        self.virtual_model_mapper.add_new_reactions_to_model(new_reactions)
        self.__virtual_model.add_reactions(new_reactions)

    def __add_new_reactions_to_model(self,new_reactions):
        self.mapper.add_new_reactions_to_model(new_reactions)
        self.model.add_reactions(new_reactions)

    def get_targets_to_replace_by_type(self, target: int) -> dict:
        """
        Identifies targets components, precursors and successors in the biosynthesis network. Returns a dictionary
        with the precursors' generic parent as keys and a list of structurally defined precursors as values.

        :param int target: BOIMMG identifier
        :return dict: dictionary with the precursors' generic parent as keys and a list
        of structurally defined precursors as values
        """

        res={}
        precursors = self.__compounds_ontology.get_all_predecessors_by_ont_id_rel_type(target,"precursor_of")

        successors = self.__compounds_ontology.get_successors_by_ont_id_rel_type(target,"precursor_of")

        components = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(target,"component_of")

        for component in components:
            if component in self.fatty_acids:
                if self.generic_fatty_acid in res:
                    res.get(self.generic_fatty_acid).append(component)
                else:
                    res[self.generic_fatty_acid] = [component]

        for precursor in precursors:
            parent = self.__compounds_ontology.get_successors_by_ont_id_rel_type(precursor,"is_a")

            if parent[0] in res:
                res.get(parent[0]).append(precursor)
            else:
                res[parent[0]] = [precursor]

        for successor in successors:
            parent = self.__compounds_ontology.get_successors_by_ont_id_rel_type(successor,"is_a")

            if parent[0] in res:
                res.get(parent[0]).append(successor)
            else:
                res[parent[0]] = [successor]


        return res

    # def find_components_reactions_and_generalize(self):
    #
    #     candidate_reactions = []
    #
    #     res = {}
    #
    #     generic_compounds = []
    #     pseudo_components = []
    #     components = []
    #     components_boimmg_id = []
    #
    #     for compound in self.mapper.boimmg_db_model_map:
    #
    #         boimmg_id =self.mapper.boimmg_db_model_map[compound]
    #         if boimmg_id in self.__components:
    #             model_compound = self.model.metabolites.get_by_id(compound)
    #             pseudo_components.append(compound)
    #             candidate_reactions.extend(list(model_compound.reactions))
    #
    #             if boimmg_id in self.fatty_acids:
    #                 components.append(compound)
    #
    #         node = self.__compounds_ontology.get_node_by_ont_id(self.mapper.boimmg_db_model_map[compound])
    #
    #         if node.generic:
    #             generic_compounds.append(compound)
    #
    #     for reaction in candidate_reactions:
    #
    #         regex = re.compile("[Bb][Ii][Oo][Mm][Aa][Ss][Ss]")
    #
    #         if reaction.annotation != "SBO:0000629" and not regex.match(reaction.name):
    #
    #             reactants = reaction.reactants
    #             products = reaction.products
    #
    #             reactant_found = False
    #             product_found = False
    #
    #             pseudo_component_in_reactants = False
    #
    #             isComponent = False
    #
    #             for reactant in reactants:
    #                 reactant_id = reactant.id
    #
    #                 if reactant_id in generic_compounds:
    #                     reactant_found = True
    #
    #                 elif reactant_id in pseudo_components:
    #                     pseudo_component_in_reactants = True
    #
    #                     if reactant_id in components:
    #                         isComponent = True
    #
    #             for product in products:
    #                 product_id = product.id
    #                 if product_id in generic_compounds:
    #                     product_found = True
    #
    #                 elif product_id in components:
    #                     isComponent = True
    #
    #             if reactant_found or product_found:
    #
    #                 new_reaction, subtracted = \
    #                     self.remove_components_from_reaction(reaction,pseudo_components,pseudo_component_in_reactants,isComponent)
    #
    #                 for sub in subtracted:
    #                     if sub not in components_boimmg_id:
    #                         components_boimmg_id.append(sub)
    #
    #
    #                 res[reaction.id] = (new_reaction,isComponent,pseudo_component_in_reactants,subtracted)
    #
    #     return res,components_boimmg_id

    # def remove_components_from_reaction(self,reaction,pseudo_components,pseudo_component_in_reactants,isComponent):
    #
    #     new_reaction = deepcopy(reaction)
    #
    #     component_parent = None
    #     if pseudo_component_in_reactants:
    #
    #         reactants = list(new_reaction.reactants)
    #
    #         to_subtract = {}
    #
    #         while reactants:
    #             reactant = reactants.pop()
    #             if reactant.id in self.mapper.boimmg_db_model_map.keys():
    #                 if reactant.id in pseudo_components:
    #
    #                     if not component_parent:
    #                         child = self.mapper.get_boimmg_id_from_model_compound_id(reactant.id)
    #                         parent = self.__compounds_ontology.get_successors_by_ont_id_rel_type(child, "is_a")
    #                         component_parent = parent[0]
    #
    #                     coef = new_reaction.metabolites[reactant]
    #                     to_subtract[reactant] = coef
    #
    #         new_reaction.subtract_metabolites(to_subtract)
    #
    #         components = []
    #         for subtracted in to_subtract:
    #             if isComponent:
    #                 boimmg_id = self.mapper.get_boimmg_id_from_model_compound_id(subtracted.id)
    #                 components.append(boimmg_id)
    #
    #             else:
    #                 boimmg_id = self.mapper.get_boimmg_id_from_model_compound_id(subtracted.id)
    #                 new_components = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(boimmg_id,
    #                                                                                                "component_of")
    #                 components.extend(new_components)
    #
    #         if isComponent:
    #             self.add_boimmg_metabolites_to_reaction(self.generic_fatty_acid, -1, new_reaction, reactant.compartment)
    #
    #         else:
    #
    #             self.add_boimmg_metabolites_to_reaction(component_parent,-1,new_reaction,reactant.compartment)
    #
    #         return new_reaction, components
    #
    #     else:
    #         products = list(new_reaction.products)
    #
    #         to_subtract = {}
    #
    #         while products:
    #             product = products.pop()
    #
    #             if product.id in self.mapper.boimmg_db_model_map.keys():
    #                 if product.id in pseudo_components:
    #
    #                     if not component_parent:
    #
    #                         child = self.mapper.get_boimmg_id_from_model_compound_id(product.id)
    #                         parent = self.__compounds_ontology.get_successors_by_ont_id_rel_type(child, "is_a")
    #
    #                         component_parent = parent[0]
    #
    #                     coef = new_reaction.metabolites[product]
    #                     to_subtract[product] = coef
    #
    #         if isComponent:
    #             self.add_boimmg_metabolites_to_reaction(self.generic_fatty_acid, 1, new_reaction, product.compartment)
    #
    #         else:
    #
    #             self.add_boimmg_metabolites_to_reaction(component_parent, 1, reaction, new_reaction.compartment)
    #
    #         new_reaction.subtract_metabolites(to_subtract)
    #
    #         components = []
    #         for subtracted in to_subtract:
    #             if isComponent:
    #                 boimmg_id = self.mapper.get_boimmg_id_from_model_compound_id(subtracted.id)
    #                 components.append(boimmg_id)
    #
    #             else:
    #                 boimmg_id = self.mapper.get_boimmg_id_from_model_compound_id(subtracted.id)
    #                 new_components = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(boimmg_id,
    #                                                                                                "component_of")
    #                 components.extend(new_components)
    #
    #
    #         return new_reaction, components


    # def __get_children_in_model(self,product):
    #
    #     children_in_model = {}
    #
    #     parent_boimmg_id = self.virtual_model_mapper.get_boimmg_id_from_model_compound_id(product.id)
    #
    #     if parent_boimmg_id:
    #         children = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(parent_boimmg_id,"is_a")
    #
    #         for mapped_compound in self.virtual_model_mapper.boimmg_db_model_map:
    #             eventual_child = self.virtual_model_mapper.boimmg_db_model_map[mapped_compound]
    #
    #             if eventual_child in children:
    #
    #                 children_in_model[eventual_child] = self.virtual_model_mapper.boimmg_db_model_map[mapped_compound]
    #
    #     return children_in_model

    # def check_if_decomposition_reaction(self,reaction:Reaction):
    #
    #     reactants = reaction.reactants
    #     products = reaction.products
    #
    #     inReactants = self.__check_component_in_reaction_side(reactants)
    #
    #     inProducts = self.__check_component_in_reaction_side(products)
    #
    #     if inReactants or inProducts:
    #         return True,inReactants
    #
    #     else:
    #         return False, inReactants



    # def __check_component_in_reaction_side(self,metabolites):
    #     for metabolite in metabolites:
    #
    #         boimmg_id = self.virtual_model_mapper.get_boimmg_id_from_model_compound_id(metabolite.id)
    #
    #         if boimmg_id:
    #             if boimmg_id == self.generic_fatty_acid or boimmg_id == self.__generic_acps or \
    #                     boimmg_id == self.__generic_acyl_coa:
    #                 return True
    #
    #     else:
    #         return False

    def __set_components(self):

        self.fatty_acids = self.__compounds_ontology.get_leaves_from_ont_id(527)
        self.generic_fatty_acid = 527

        # self.__acyl_coas = self.__compounds_ontology.get_leaves_from_ont_id(43)
        # self.__generic_acyl_coa = 43
        #
        # self.__acyl_acps = self.__compounds_ontology.get_leaves_from_ont_id(704832)
        # self.__generic_acps = 704832

        # self.__components = self.fatty_acids + self.__acyl_acps + self.__acyl_coas
        # self.__components = self.fatty_acids

    # def solve_components_reactions(self,same_components=False,components_reactions = {}, components = []):
    #
    #     if not components_reactions and not components:
    #         components_reactions,components = self.find_components_reactions_and_generalize()
    #
    #     new_reactions = [components_reactions[reaction][0] for reaction in components_reactions.keys()]
    #
    #     self.__add_new_reactions_to_virtual_model(new_reactions)
    #
    #     self.__solve_components_reactions(components_reactions,components,same_components)

    def __solve_components_reactions(self,components_reactions,components,same_components):


        for reaction in components_reactions:

            self.__correct_component_reaction(components_reactions[reaction],components,same_components)

    # def __correct_component_reaction(self,reaction_features,components,same_components):
    #
    #     new_reaction = reaction_features[0]
    #     isReactant = reaction_features[2]
    #     deleted = reaction_features[3]
    #
    #     reactants = new_reaction.reactants
    #     products = new_reaction.products
    #
    #     # r_groups_reactants,r_groups_products = self.get_r_groups_from_reaction_compounds(reactants,products)
    #
    #     if isReactant:
    #
    #         boimmg_compound = self.find_boimmg_compound(products)
    #
    #         if boimmg_compound:
    #
    #             self.granulate(boimmg_compound, components, same_components)
    #
    #     else:
    #         boimmg_compound = self.find_boimmg_compound(reactants)
    #
    #         if boimmg_compound:
    #             self.granulate(boimmg_compound, components, same_components,decomposition=True)


    def find_boimmg_compound(self,compounds):
        """
        It finds only one boimmg compound within a list of model compounds
        :param compounds:
        :return:
        """

        for compound in compounds:

            boimmg_id = self.mapper.get_boimmg_id_from_model_compound_id(compound.id)

            if boimmg_id:
                return boimmg_id

        return None


    # def get_r_groups_from_reaction_compounds(self,reactants,products):
    #
    #     r_groups_reactants = 0
    #     r_groups_products = 0
    #
    #     for reactant in reactants:
    #         r_groups_reactants += self.get_r_groups_from_compound(reactant)
    #
    #     for product in products:
    #         r_groups_products += self.get_r_groups_from_compound(product)
    #
    #     return r_groups_reactants,r_groups_products

    # def get_r_groups_from_compound(self,model_compound):
    #
    #     boimmg_id = self.virtual_model_mapper.get_boimmg_id_from_model_compound_id(model_compound.id)
    #     if boimmg_id:
    #         node = self.__compounds_ontology.get_node_by_ont_id(boimmg_id)
    #         smiles = node.smiles
    #
    #         if smiles:
    #             return smiles.count("*")
    #
    #         return 0
    #     return 0

    def filter_targets_to_replace(self,targets: List[int],parent: int) -> List[int]:
        """
        Filter the predicted precursors by side chains.

        :param List[int] targets:
        :param int parent:
        :return List[int]: filtered list
        """

        res =  []

        for target in targets:
            components = self.__compounds_ontology.get_predecessors_by_ont_id_rel_type(target,"component_of")

            node = self.__compounds_ontology.get_node_by_ont_id(target)

            child_smiles = node.smiles

            parent_node = self.__compounds_ontology.get_node_by_ont_id(parent)

            parent_smiles = parent_node.smiles

            molecules = chemo_utilities.get_side_chains(parent_smiles, child_smiles)

            if molecules:
                new_molecules = []
                new_molecules_n = 0
                for molecule in molecules:

                    if molecule not in new_molecules:
                        new_molecules_n+=1
                        new_molecules.append(molecule)

                if new_molecules_n == len(components):
                    res.append(target)

        return res