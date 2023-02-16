import itertools

from biocyc import biocyc
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MolToInchiKey, MolToSmiles, rdmolops
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles

from boimmgpy.database.databases_babel import ModelSEEDDatabases
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.service.network_handlers.pathway_handler import PathwayHandler
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.id_converters.reactions_id_converter import ReactionsIDConverter
from boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from boimmgpy.model_seed.model_seed_reactions_database import ModelSeedReactionsDB
from boimmgpy.utilities import chemo_utilities, file_utilities
from boimmgpy.definitions import EXCEPTIONS


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
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="", flush=True)
    # Print New Line on Complete
    if iteration == total:
        print()


class TransformationsHandler:

    def __init__(self):
        biocyc.set_organism('meta')
        self.exceptions_conf = file_utilities.read_conf_file(EXCEPTIONS)
        self.reactions_transformations = {}

        self.compounds_simplifier = {"cpd24471": "cpd12300",
                                     "cpd27426": "cpd27914",
                                     "cpd11611": "cpd29330",
                                     "cpd00487": "cpd19000"}

    @property
    def reactions_transformations(self):
        return self.__reactions_transformations

    @reactions_transformations.setter
    def reactions_transformations(self, reactions_transformations):
        self.__reactions_transformations = reactions_transformations

    @property
    def modelSeedCompoundsDB(self):
        return self.__modelSeedCompDB

    @modelSeedCompoundsDB.setter
    def modelSeedCompoundsDB(self, db):
        self.__modelSeedCompDB = db

    @property
    def compoundsIDConverter(self):
        return self.__compoundsIDConverter

    @compoundsIDConverter.setter
    def compoundsIDConverter(self, compoundsIDConverter):
        self.__compoundsIDConverter = compoundsIDConverter

    def extract_metacyc_pathways(self, parent, generic=False):

        biocyc.set_organism('meta')
        ptw_class = biocyc.get(parent)
        l = [ptw_class]
        res = []
        while len(l) > 0:
            try:
                node = l.pop(0)

                if not node.subclasses and not node.instances:
                    res.append(node)
                    # print(node.name)

                following = node.subclasses
                following.extend(node.instances)
                for elem in following:
                    if elem not in res and elem not in l:
                        l.insert(0, elem)

            except:
                pass

        i = 0
        for instance in res:
            try:
                # printProgressBar(i, len(res))
                self.__pathway_handler.add_metacyc_pathway(instance.id, generic)
                i += 1

            except:
                pass

    def __simplify_compounds(self, modelseed_compounds):
        res = []
        for compound in modelseed_compounds:
            if compound in self.compounds_simplifier.keys():
                res.append(self.compounds_simplifier.get(compound))
            else:
                res.append(compound)

        return res

    def choose_transformations_from_pathway(self, pathway_parent, compounds_to_avoid=[], accept_cycles=False):

        self.modelSeedCompoundsDB = ModelSeedCompoundsDB()
        self.modelSeedReactionsDB = ModelSeedReactionsDB()
        self.reactionsIDConverter = ReactionsIDConverter()

        self.__pathway_handler = PathwayHandler(self.modelSeedCompoundsDB)

        self.extract_metacyc_pathways(pathway_parent, True)
        target_pathways = self.__pathway_handler.get_target_pathways()

        for pathway in target_pathways:

            paths_within_pathway = self.__pathway_handler.get_all_paths_using_pathway_id(pathway, accept_cycles)
            for path in paths_within_pathway:

                reverse_pathway = path[::-1]

                for reaction in reverse_pathway:

                    print(reaction)

                    if reaction not in self.reactions_transformations.keys():

                        reactions = self.reactionsIDConverter.convert_db_id_to_model_seed_by_db_id(reaction)

                        found = False
                        for react in reactions:
                            temp_reaction = self.modelSeedReactionsDB.getReaction(react)
                            compounds = temp_reaction.getCompounds()
                            for compound in compounds:
                                temp_compound = self.modelSeedCompoundsDB.get_compound_by_id(compound)
                                if "R" in temp_compound.getFormula() or "*" in temp_compound.getSmiles():
                                    model_seed_reaction = temp_reaction
                                    found = True
                                    break
                            if found:
                                break

                        # model_seed_reaction = self.modelSeedReactionsDB.getReaction(reaction_id[0])
                        if found:
                            reactants = model_seed_reaction.getReactants()
                            if "cpd26685" in reactants:
                                reactants.remove("cpd26685")
                                reactants.append("cpd00167")

                            if "cpd28189" in reactants:
                                reactants.remove("cpd28189")
                                reactants.append("cpd00431")

                            products = model_seed_reaction.getProducts()
                            if "cpd26685" in products:
                                products.remove("cpd26685")
                                products.append("cpd00167")

                            if "cpd28189" in products:
                                products.remove("cpd28189")
                                products.append("cpd00431")

                            reactants = self.__simplify_compounds(reactants)
                            products = self.__simplify_compounds(products)

                            model_seed_reactants = []

                            for reac in reactants:
                                model_seed_reactant = self.modelSeedCompoundsDB.get_compound_by_id(reac)

                                if "R" in model_seed_reactant.getFormula() or "*" in model_seed_reactant.getSmiles():

                                    if (model_seed_reactant.getSmiles() or
                                        model_seed_reactant.getDbId() in self.exceptions_conf.keys()) and \
                                            model_seed_reactant.getDbId() not in compounds_to_avoid:
                                        model_seed_reactants.append(model_seed_reactant)

                            model_seed_products = []
                            for prod in products:

                                model_seed_product = self.modelSeedCompoundsDB.get_compound_by_id(prod)

                                if "R" in model_seed_product.getFormula() or "*" in model_seed_product.getSmiles():

                                    if (model_seed_product.getSmiles() or
                                        model_seed_product.getDbId() in self.exceptions_conf.keys()) and \
                                            model_seed_product.getDbId() not in compounds_to_avoid:
                                        model_seed_products.append(model_seed_product)

                            # products_model_seed = self.__convert_metacyc_compound_into_model_seed(products)
                            # reactants_model_seed = self.__convert_metacyc_compound_into_model_seed(reactants)

                            if model_seed_products:
                                transformation, _ = self.__convert_reaction_into_transformation(model_seed_reactants,
                                                                                                model_seed_products)

                                self.reactions_transformations[reaction] = transformation

    def save_transformations(self, filename):

        with open(filename, "w") as file:
            for reaction in self.reactions_transformations.keys():
                file.write(reaction + "!!!" + self.reactions_transformations[reaction] + "\n")

    def load_transformations(self, filepath):

        with open(filepath, "r") as file:

            lines = file.readlines()

            for line in lines:
                if line:
                    line_lst = line.split("!!!")
                    reaction = line_lst[0]
                    transformation = line_lst[1].strip()

                    self.reactions_transformations[reaction] = transformation

    def __convert_metacyc_compound_into_model_seed(self, metacyc_compounds):
        res = []
        for metacyc_compound in metacyc_compounds:
            modelseedids = self.compoundsIDConverter.convert_dbID_into_modelSeedId("MetaCyc", metacyc_compound.id)
            if modelseedids:
                model_seed_compound = self.__check_if_generic(modelseedids)
                if model_seed_compound:
                    res.append(model_seed_compound)

        return res

    def __neutralize_molecules(self, molecules):
        molecules_smiles = []

        for molecule in molecules:
            model_seed_id = molecule.getDbId()
            if model_seed_id in self.exceptions_conf.keys():
                new_smiles = self.exceptions_conf.get(model_seed_id)
            else:
                new_smiles, neutralised = chemo_utilities.NeutraliseCharges(molecule.getSmiles())

            molecules_smiles.append(new_smiles)

        return molecules_smiles

    def __convert_reaction_into_transformation(self, reactants, products):
        """
        This method aims at transforming a MetaCyc reaction into a SMIRK transformation

        :param reactants: list of modelseed generic reactants
        :param products: list of modelseed generic products
        :return: a string representing the transformation
        """

        reactants_smiles = self.__neutralize_molecules(reactants)
        products_smiles = self.__neutralize_molecules(products)

        react = ".".join(reactants_smiles)
        prod = ".".join(products_smiles)

        r_groups_reactants = react.count("*")
        r_groups_products = prod.count("*")

        n_primary_precursors = 0
        while r_groups_reactants < r_groups_products:
            reactants_smiles.append(reactants_smiles[0])  # martelado

            react = ".".join(reactants_smiles)
            r_groups_reactants = react.count("*")
            n_primary_precursors += 1

        reactant_smiles_combinations = chemo_utilities.get_combinations_of_dummy_atoms(react)
        product_smiles = chemo_utilities.convert_model_seed_dummy_atoms(prod)

        transformations_combinations = []
        for reactant_smiles in reactant_smiles_combinations:
            reaction_smart = product_smiles + ">>" + reactant_smiles
            transformations_combinations.append(reaction_smart)

        user_interaction_string = ""
        i = 1
        dict_for_user_to_choose = {}

        if len(transformations_combinations) > 1:

            for transformation in transformations_combinations:
                user_interaction_string += str(i) + ": "
                user_interaction_string += transformation + "\n"
                dict_for_user_to_choose[i] = transformation
                i += 1

            transformation_i = input(user_interaction_string)
            chosen_transformation = dict_for_user_to_choose[int(transformation_i)]

        else:
            chosen_transformation = transformations_combinations[0]

        return chosen_transformation, n_primary_precursors

    def __check_if_generic(self, modelseedids):
        for modelseedid in modelseedids:

            compound = self.modelSeedCompoundsDB.get_compound_by_id(modelseedid)
            formula = compound.getFormula()
            if "R" in formula:
                return compound
        return None


class OntologyGeneratorDA:

    def __init__(self):

        self.cores = None
        self.modelSeedCompoundsDB = ModelSeedCompoundsDB()
        self.modelSeedReactionsDB = ModelSeedReactionsDB()
        biocyc.set_organism("meta")
        # self.__modelSeedCompoundsDatabase = ModelSeedCompoundsDB()
        self.accessor = CompoundsDBAccessor()

        self.compoundsIDConverter = CompoundsIDConverter()
        self.reactionsIDConverter = ReactionsIDConverter()

        self.compounds_simplifier = {"cpd24471": "cpd12300", "cpd27426": "cpd27914", "cpd11611": "cpd29330",
                                     "cpd00487": "cpd19000", "cpd26685": "cpd00167", "cpd28189": "cpd00431"}

        # ############################# exceptionally for biosynthetic pathways in the sphingolipids pathway
        # ##################

        # #######################################
        # ###############################################################################

        self.exceptions_conf = file_utilities.read_conf_file(EXCEPTIONS)

    @property
    def modelSeedCompoundsDB(self):
        return self.__modelSeedCompDB

    @modelSeedCompoundsDB.setter
    def modelSeedCompoundsDB(self, db):
        self.__modelSeedCompDB = db

    @property
    def modelSeedReactionsDB(self):
        return self.__modelSeedReactionsDB

    @modelSeedReactionsDB.setter
    def modelSeedReactionsDB(self, db):
        self.__modelSeedReactionsDB = db

    @property
    def transformationsHandler(self):
        return self.__transformationsHandler

    @transformationsHandler.setter
    def transformationsHandler(self, transformationsHandler):
        self.__transformationsHandler = transformationsHandler

    @property
    def cores(self):
        return self.__cores

    @cores.setter
    def cores(self, core):
        self.__cores = core

    @property
    def compoundsIDConverter(self):
        return self.__compoundsIDConverter

    @compoundsIDConverter.setter
    def compoundsIDConverter(self, compoundsIDConverter):
        self.__compoundsIDConverter = compoundsIDConverter

    @property
    def reactionsIDConverter(self):
        return self.__reactionsIDConverter

    @reactionsIDConverter.setter
    def reactionsIDConverter(self, reactionsIDConverter):
        self.__reactionsIDConverter = reactionsIDConverter

    @property
    def accessor(self):
        return self.__accessor

    @accessor.setter
    def accessor(self, accessor):
        self.__accessor = accessor

    @property
    def compounds_simplifier(self):
        return self.__compounds_simplifier

    @compounds_simplifier.setter
    def compounds_simplifier(self, value):
        self.__compounds_simplifier = value

    def get_combinations_with_the_same_components(self, combinations):
        res = []
        i = 0

        comp_and_components = {}
        print("Getting combinations... ")
        for combination in combinations:
            i += 1

            printProgressBar(i, len(combinations))

            if combination[0] not in comp_and_components:
                temp_components = self.accessor.get_predecessors_by_ont_id_rel_type(combination[0], "component_of")
                comp_and_components[combination[0]] = temp_components

            else:
                temp_components = comp_and_components[combination[0]]

            found = False
            for compound in combination[1:]:
                if compound not in comp_and_components.keys():
                    components = self.accessor.get_predecessors_by_ont_id_rel_type(compound, "component_of")
                    comp_and_components[compound] = components

                else:
                    components = comp_and_components[compound]

                if components:
                    component = components[0]
                    if component in temp_components:
                        found = True
                    else:
                        found = False
                        break
            if found:
                res.append(combination)
        return res

    def create_new_compounds_using_transformation(self, transformation, reactants, products, reaction_id="",
                                                  source="", add_rel=False, same_components=False,
                                                  reverse=True):
        children = []
        for reactants_ont_id in reactants:

            if same_components:
                reactant_children = self.accessor.get_predecessors_with_same_component(reactants_ont_id, "is_a")

            else:
                reactant_children = self.accessor.get_predecessors_by_ont_id_rel_type(reactants_ont_id, "is_a")

            children.append(reactant_children.copy())

        combinations = list(itertools.product(*children))
        if same_components and len(children) > 1:
            combinations = self.get_combinations_with_the_same_components(combinations)

        i = 0
        for combination in combinations:
            printProgressBar(i + 1, len(combinations), prefix='Progress:', suffix='Complete')
            i += 1

            new_molecules_list = self.generate_new_molecule_with_transformation(list(combination),
                                                                                transformation)

            for new_molecule in new_molecules_list:
                self.__add_molecule_and_relationships(products,
                                                      new_molecule,
                                                      reaction_id,
                                                      source=source,
                                                      add_biosynthetic_rel=add_rel, reverse=reverse)

    def handle_manual_reaction(self, reaction_transformation, reactants, products, reaction_id, source):

        model_seed_products = []
        for product in products:
            model_seed_compound = self.modelSeedCompoundsDB.get_compound_by_id(product)
            model_seed_products.append(model_seed_compound)

        # products_model_seed = self.__convert_metacyc_compounds_into_model_seed(products)

        if model_seed_products:

            products_ont_ids = self.__get_generic_compounds_from_ontology(model_seed_products)

            children = []
            for products_ont_id in products_ont_ids:
                product_children = self.accessor.get_predecessors_by_ont_id_rel_type(products_ont_id, "is_a")
                children.append(product_children)

            combinations = list(itertools.product(*children))

            i = 0
            for combination in combinations:

                printProgressBar(i + 1, len(children), prefix='Progress:', suffix='Complete')
                i += 1

                reactants = [self.modelSeedCompoundsDB.get_compound_by_id(reac) for reac in reactants]
                products = [self.modelSeedCompoundsDB.get_compound_by_id(prod) for prod in products]

                new_reaction_products = self.__check_if_product_is_compatible_with_ont_container(products,
                                                                                                 combination)

                reactants_ont_ids = self.__get_generic_compounds_from_ontology(reactants)

                for reactants_ont_id in reactants_ont_ids:
                    for products_ont_id in products_ont_ids:
                        self.accessor.add_biosynthetic_relationship(reactants_ont_id, products_ont_id,
                                                                    reaction=reaction_id)

                if reactants_ont_ids and new_reaction_products:

                    new_molecules_list = self.generate_new_molecule_with_transformation(new_reaction_products,
                                                                                        reaction_transformation)

                    for new_molecules in new_molecules_list:

                        if new_molecules:
                            self.__add_molecule_and_relationships(reactants_ont_ids,
                                                                  new_molecules,
                                                                  reaction_id,
                                                                  source)

    def handle_reaction(self, reaction, pathway_id, transformations_file):
        self.transformationsHandler = TransformationsHandler()
        self.transformationsHandler.load_transformations(transformations_file)

        modelseedid = self.reactionsIDConverter.convert_db_id_to_model_seed_by_db_id(reaction)
        model_seed_reaction = self.modelSeedReactionsDB.getReaction(modelseedid[0])

        # products = metacyc_reaction.compounds_right
        products = model_seed_reaction.getProducts()
        model_seed_products = []

        for product in products:
            model_seed_compound = self.modelSeedCompoundsDB.get_compound_by_id(product)
            model_seed_products.append(model_seed_compound)

        # products_model_seed = self.__convert_metacyc_compounds_into_model_seed(products)

        if model_seed_products:

            products_ont_ids = self.__get_generic_compounds_from_ontology(model_seed_products)

            children = []
            for products_ont_id in products_ont_ids:
                product_children = self.accessor.get_predecessors_by_ont_id_rel_type(products_ont_id, "is_a")

                children.append(product_children)

            combinations = list(itertools.product(*children))

            i = 0
            for combination in combinations:

                printProgressBar(i + 1, len(children), prefix='Progress:', suffix='Complete')
                i += 1

                reaction_id = self.reactionsIDConverter.convert_db_id_to_model_seed_by_db_id(reaction)
                model_seed_reaction = self.modelSeedReactionsDB.getReaction(reaction_id[0])
                model_seed_reactants = model_seed_reaction.getReactants()
                model_seed_products = model_seed_reaction.getProducts()

                model_seed_reactants = self.__simplify_compounds(model_seed_reactants)
                model_seed_products = self.__simplify_compounds(model_seed_products)

                reactants = [self.modelSeedCompoundsDB.get_compound_by_id(reac) for reac in model_seed_reactants]
                products = [self.modelSeedCompoundsDB.get_compound_by_id(prod) for prod in model_seed_products]

                new_reaction_products = self.__check_if_product_is_compatible_with_ont_container(products,
                                                                                                 combination)

                reactants_ont_ids = self.__get_generic_compounds_from_ontology(reactants)

                for reactants_ont_id in reactants_ont_ids:
                    for products_ont_id in products_ont_ids:
                        self.accessor.add_biosynthetic_relationship(reactants_ont_id, products_ont_id,
                                                                    reaction=reaction,
                                                                    pathway=pathway_id)

                if reactants_ont_ids and new_reaction_products:
                    transformation = self.transformationsHandler.reactions_transformations[reaction]

                    new_molecules_list = self.generate_new_molecule_with_transformation(new_reaction_products,
                                                                                        transformation)

                    for new_molecules in new_molecules_list:

                        if new_molecules:
                            self.__add_molecule_and_relationships(reactants_ont_ids,
                                                                  new_molecules,
                                                                  reaction,
                                                                  pathway_id)

    def __call__(self, pathway_parent, cores, transformations_file, pathway_to_avoid=None, reactions_to_avoid=None,
                 accept_cycles=False, list_of_targets=None, all=True):

        if reactions_to_avoid is None:
            reactions_to_avoid = []
        if pathway_to_avoid is None:
            pathway_to_avoid = []
        if list_of_targets is None:
            list_of_targets = []
        self.transformationsHandler = TransformationsHandler()
        self.transformationsHandler.load_transformations(transformations_file)
        self.seen = pathway_to_avoid

        self.reactions_to_avoid = reactions_to_avoid
        self.cores = cores
        self.__pathway_handler = PathwayHandler(self.modelSeedCompoundsDB)
        self.extract_metacyc_pathways(pathway_parent, True)
        target_pathways = self.__pathway_handler.get_target_pathways()

        for pathway in target_pathways:
            print("Trying to handle this pathway: %s" % pathway)
            if pathway not in self.seen:
                paths_within_pathways = self.__pathway_handler.get_all_paths_using_pathway_id(pathway, accept_cycles)

                for path in paths_within_pathways:
                    to_continue = True
                    for reaction in self.reactions_to_avoid:
                        if reaction in path:
                            to_continue = False
                            break

                    if to_continue:
                        reverse_path = path[::-1]
                        print("Target: %s" % reverse_path[0])
                        self.__handle_pathway(reverse_path, pathway, list_of_targets, all)

                # self.seen.extend(reverse_pathway)

    def __handle_pathway(self, reverse_pathway, pathway_id, list_of_targets=None, all=True):

        if list_of_targets is None:
            list_of_targets = []
        print("Handling the following pathway: %s" % pathway_id)
        first_reaction = reverse_pathway[0]

        # metacyc_reaction = biocyc.get(first_reaction)
        modelseedid = self.reactionsIDConverter.convert_db_id_to_model_seed_by_db_id(first_reaction)
        model_seed_reaction = self.modelSeedReactionsDB.getReaction(modelseedid[0])

        # products = metacyc_reaction.compounds_right
        products = model_seed_reaction.getProducts()
        model_seed_products = []

        for product in products:
            model_seed_compound = self.modelSeedCompoundsDB.get_compound_by_id(product)
            model_seed_products.append(model_seed_compound)

        # products_model_seed = self.__convert_metacyc_compounds_into_model_seed(products)

        if model_seed_products:

            products_ont_ids = self.__get_generic_compounds_from_ontology(model_seed_products)

            i = 0

            if not list_of_targets:
                children = []

                for products_ont_id in products_ont_ids:

                    if not all:
                        product_children = self.accessor.get_predecessors_only_from_lipid_maps_and_ms(products_ont_id)

                    else:
                        product_children = self.accessor.get_predecessors_by_ont_id(products_ont_id)

                    children.append(product_children)

            else:
                children = [list_of_targets]

            combinations = list(itertools.product(*children))

            for combination in combinations:

                printProgressBar(i + 1, len(combinations), prefix='Progress:', suffix='Complete')
                i += 1

                children_products_ont_ids = list(combination)
                for reaction in reverse_pathway:

                    # metacyc_reaction = biocyc.get(reaction)
                    # print(metacyc_reaction.id)
                    # reactants = metacyc_reaction.compounds_left
                    # products = metacyc_reaction.compounds_right

                    reaction_id = self.reactionsIDConverter.convert_db_id_to_model_seed_by_db_id(reaction)
                    model_seed_reaction = self.modelSeedReactionsDB.getReaction(reaction_id[0])
                    model_seed_reactants = model_seed_reaction.getReactants()
                    model_seed_products = model_seed_reaction.getProducts()

                    model_seed_reactants = self.__simplify_compounds(model_seed_reactants)
                    model_seed_products = self.__simplify_compounds(model_seed_products)

                    reactants = [self.modelSeedCompoundsDB.get_compound_by_id(reac) for reac in model_seed_reactants]
                    products = [self.modelSeedCompoundsDB.get_compound_by_id(prod) for prod in model_seed_products]

                    if children_products_ont_ids:

                        for child in children_products_ont_ids:
                            new_reaction_products = self.__check_if_product_is_compatible_with_ont_container(products,
                                                                                                             [child])

                            # reactants_model_seed = self.__convert_metacyc_compounds_into_model_seed(reactants)
                            reactants_ont_ids = self.__get_generic_compounds_from_ontology(reactants)

                            if reaction not in self.seen:
                                for reactants_ont_id in reactants_ont_ids:
                                    for products_ont_id in products_ont_ids:
                                        self.accessor.add_biosynthetic_relationship(reactants_ont_id, products_ont_id,
                                                                                    reaction=reaction,
                                                                                    pathway=pathway_id)

                                self.seen.append(reaction)

                            if reactants_ont_ids and new_reaction_products:

                                if reaction in self.transformationsHandler.reactions_transformations:
                                    transformation = self.transformationsHandler.reactions_transformations[reaction]

                                    # the following describe the reaction in a reverse way,
                                    # the products of the actual reaction will be the reactants in the transformation

                                    new_molecules_list = self.generate_new_molecule_with_transformation(
                                        new_reaction_products, transformation)

                                    new_products = []

                                    for new_molecules in new_molecules_list:

                                        if new_molecules:
                                            new_molecules_ont_ids = self.__add_molecule_and_relationships(
                                                reactants_ont_ids,
                                                new_molecules, reaction, pathway_id)
                                            new_products.extend(new_molecules_ont_ids)

                                    children_products_ont_ids = new_products.copy()

    def generate_new_molecule_with_transformation(self, reactants, transformation):
        """

        :param reactants: ontology ids of the different possible reactants
        :param transformation: SMART transformation
        :return list: list of tuples (ontology_id of downstream compound, new molecules smiles)
        """
        res = []
        reactants_mol = []
        for reactant in reactants:

            reactant = self.accessor.get_node_by_ont_id(reactant)

            if not reactant.smiles:
                return []

            neutralized_smiles, _ = chemo_utilities.NeutraliseCharges(reactant.smiles)
            reactant_mol = MolFromSmiles(neutralized_smiles)

            reactants_mol.append(reactant_mol)

        reactants_tuple = tuple(reactants_mol)

        transformation_smarts = AllChem.ReactionFromSmarts(transformation)

        try:
            new_molecules = [Chem.MolToSmiles(x, 1) for x in transformation_smarts.RunReactants(
                reactants_tuple)[0]]

            res.append((reactants, new_molecules))
        except:
            res.append(None)

        return res

    def extract_metacyc_pathways(self, parent, generic=False):

        biocyc.set_organism('meta')
        ptw_class = biocyc.get(parent)
        l = [ptw_class]
        res = []
        while len(l) > 0:
            try:
                node = l.pop(0)

                if not node.subclasses and not node.instances:
                    res.append(node)
                    # print(node.name)

                following = node.subclasses
                following.extend(node.instances)
                for elem in following:
                    if elem not in res and elem not in l:
                        l.insert(0, elem)

            except:
                pass

        i = 0
        for instance in res:
            try:
                # printProgressBar(i, len(res))
                self.__pathway_handler.add_metacyc_pathway(instance.id, generic)
                i += 1

            except:
                pass

    @staticmethod
    def get_similarity(molecules, reactant):
        similarities = []
        for molecule in molecules:
            fgp1 = AllChem.GetMorganFingerprint(MolFromSmiles(molecule), 3, useFeatures=True, useChirality=True)
            fgp_reactant = AllChem.GetMorganFingerprint(reactant, 1, useFeatures=True, useChirality=True)
            similarity = DataStructs.TanimotoSimilarity(fgp_reactant, fgp1, useFeatures=True, useChirality=True)

            similarities.append(similarity)

        average = sum(similarities) / len(similarities)
        return average

    @staticmethod
    def __get_filtered_compounds(smiles_keys, smarts, smiles_database):
        new_smiles_list = smiles_keys.copy()
        filtered_compounds = []
        for smile in smiles_keys:
            molecule = MolFromSmiles(smile)
            if molecule:
                match = molecule.GetSubstructMatch(smarts)

                if match:
                    filtered_compounds.append(smiles_database[smile])
                    new_smiles_list.remove(smile)

        smiles_keys = new_smiles_list.copy()

        return (filtered_compounds, smiles_keys)

    def __convert_metacyc_compounds_into_model_seed(self, metacyc_compounds):
        res = []
        for metacyc_compound in metacyc_compounds:
            modelseedids = self.__compoundsIDConverter.convert_dbID_into_modelSeedId("MetaCyc", metacyc_compound.id)
            if modelseedids:
                model_seed_compound = self.__check_if_generic(modelseedids)
                if model_seed_compound:
                    res.append(model_seed_compound)

        return res

    def __check_if_generic(self, modelseedids):
        for modelseedid in modelseedids:

            compound = self.modelSeedCompoundsDB.get_compound_by_id(modelseedid)
            formula = compound.getFormula()
            if "R" in formula:
                return compound
        return None

    def __neutralize_molecules(self, molecules):
        molecules_smiles = []
        for molecule in molecules:
            # mol = MolFromSmiles(molecule.getSmiles())
            new_smiles, neutralised = chemo_utilities.NeutraliseCharges(molecule.getSmiles())
            molecules_smiles.append(new_smiles)

        return molecules_smiles

    def __get_generic_compounds_from_ontology(self, model_seed_molecules):

        res = []

        for compound in model_seed_molecules:
            db_id = compound.getDbId()
            smiles = compound.getSmiles()
            smiles, _ = chemo_utilities.NeutraliseCharges(smiles)
            if "*" in smiles or "R" in compound.getFormula():
                ont_id = self.accessor.get_node_from_model_seed_id(db_id).id

                res.append(ont_id)

        return res

    def __add_molecule_and_relationships(self, possible_structural_parents, precursor_successor_and_smiles, reaction_id,
                                         source, pathway_id=None, add_biosynthetic_rel=True, reverse=True):

        new_compounds_ont_ids = []
        precursor_successors = precursor_successor_and_smiles[0]
        smiles_lst = precursor_successor_and_smiles[1]

        for smiles in smiles_lst:
            mol = MolFromSmiles(smiles)

            found = False
            i = 0
            while not found and i < len(possible_structural_parents):
                structuralParentOntId = possible_structural_parents[i]
                structuralParentContainer = self.accessor.get_node_by_ont_id(structuralParentOntId)
                generic_smiles = structuralParentContainer.smiles
                isStructuralParent = \
                    chemo_utilities.check_similarity_between_generic_and_complete_representation(generic_smiles, smiles)

                if isStructuralParent:
                    found = True

                i += 1

            if found:

                newCompound_ont_id = self.__generate_complete_compound_and_establish_relationships(mol,
                                                                                                   structuralParentContainer,
                                                                                                   precursor_successors,
                                                                                                   reaction_id,
                                                                                                   pathway_id, source,
                                                                                                   add_biosynthetic_rel,
                                                                                                   reverse)

                if newCompound_ont_id not in new_compounds_ont_ids:
                    new_compounds_ont_ids.append(newCompound_ont_id)

        return new_compounds_ont_ids

    def __generate_complete_compound_and_establish_relationships(self,
                                                                 mol,
                                                                 structuralParentContainer,
                                                                 precursor_sucessors, reaction_id,
                                                                 pathway_id, source, add_biosynthetic_rel=True,
                                                                 reverse=True):

        inchikey = MolToInchiKey(mol)
        smiles = MolToSmiles(mol)
        if inchikey:
            node = self.accessor.get_compound_by_inchikey(inchikey)
        else:
            node = self.accessor.get_compound_by_smiles(smiles)

        if node and precursor_sucessors:
            new_molecule_ont_id = node.id

            if add_biosynthetic_rel:
                for successor in precursor_sucessors:
                    if pathway_id:
                        if reverse:
                            self.accessor.add_biosynthetic_relationship(new_molecule_ont_id, successor,
                                                                        reaction=reaction_id, pathway=pathway_id,
                                                                        source=source)
                        else:
                            self.accessor.add_biosynthetic_relationship(successor, new_molecule_ont_id,
                                                                        reaction=reaction_id, pathway=pathway_id,
                                                                        source=source)
                    elif reverse:
                        self.accessor.add_biosynthetic_relationship(new_molecule_ont_id, successor,
                                                                    reaction=reaction_id, source=source)
                    else:
                        self.accessor.add_biosynthetic_relationship(successor, new_molecule_ont_id,
                                                                    reaction=reaction_id, source=source)

            self.accessor.establish_structural_relationship(new_molecule_ont_id, structuralParentContainer.id)

        else:
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

            if "*" in formula:
                number_r_groups = formula.count("*")
                formula = formula.replace("*", "") + "R" + str(number_r_groups)

            inchikey = MolToInchiKey(mol)
            smiles = MolToSmiles(mol)

            side_chain_providers = self.__find_side_chain_providers(mol, structuralParentContainer, precursor_sucessors)

            name = self.__generate_new_molecule_name(structuralParentContainer, side_chain_providers)

            charge = rdmolops.GetFormalCharge(mol)

            new_molecule_ont_container = self.accessor.create_compound(name=name, formula=formula,
                                                                       inchikey=inchikey,
                                                                       smiles=smiles, annotated=False, charge=charge)
            new_molecule_ont_id = new_molecule_ont_container.id

            self.accessor.establish_structural_relationship(new_molecule_ont_container.id, structuralParentContainer.id)

            for side_chain_provider in side_chain_providers:
                self.accessor.add_relationship(side_chain_provider.id, new_molecule_ont_id, "component_of")

            if add_biosynthetic_rel:
                for successor in precursor_sucessors:
                    if pathway_id:
                        if reverse:
                            self.accessor.add_biosynthetic_relationship(new_molecule_ont_id, successor,
                                                                        reaction=reaction_id, pathway=pathway_id,
                                                                        source=source)
                        else:
                            self.accessor.add_biosynthetic_relationship(successor, new_molecule_ont_id,
                                                                        reaction=reaction_id, pathway=pathway_id,
                                                                        source=source)
                    elif reverse:
                        self.accessor.add_biosynthetic_relationship(new_molecule_ont_id, successor,
                                                                    reaction=reaction_id, source=source)
                    else:
                        self.accessor.add_biosynthetic_relationship(successor, new_molecule_ont_id,
                                                                    reaction=reaction_id, source=source)

            self.accessor.establish_structural_relationship(new_molecule_ont_id, structuralParentContainer.id)

        return new_molecule_ont_id

    @staticmethod
    def __generate_new_molecule_name(product_container, combination_list):

        names = []

        for container in combination_list:
            name = container.name
            if name not in names:
                names.append(name)

        new_name = product_container.name
        for unique_name in names:
            new_name += " (" + unique_name + ") "

        return new_name

    def __generate_complete_compound(self, lipid_compound,
                                     model_seed_compound,
                                     mol,
                                     structuralParentContainer=None,
                                     side_chain_providers=None):

        if lipid_compound:
            aliases = lipid_compound.getAliases()
            formula = lipid_compound.getFormula()
            inchikey = lipid_compound.getInchiKey()
            smiles = lipid_compound.getSmiles()
            name = lipid_compound.getName()

        if lipid_compound and model_seed_compound:
            id = model_seed_compound.getDbId()
            model_seed_aliases = self.__compoundsIDConverter.get_all_aliases_by_modelSeedID(id)
            aliases.update(model_seed_aliases)

        elif lipid_compound and not model_seed_compound:
            id = lipid_compound.getDbId()


        elif model_seed_compound:
            id = model_seed_compound.getDbId()
            model_seed_aliases = self.__compoundsIDConverter.get_all_aliases_by_modelSeedID(id)
            aliases = model_seed_aliases
            formula = model_seed_compound.getFormula()
            inchikey = model_seed_compound.getInchikey()
            smiles = model_seed_compound.getSmiles()
            name = model_seed_compound.getName()

        # elif structuralParentContainer and side_chain_providers:
        else:
            id = ""
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            inchikey = MolToInchiKey(mol)
            smiles = MolToSmiles(mol)
            aliases = {}
            name = self.__generate_new_molecule_name(structuralParentContainer, side_chain_providers)

        return id, name, formula, inchikey, smiles, aliases

    def __find_side_chain_providers(self, mol, structuralParentContainer, precursor_sucessors):

        r_groups_n1 = structuralParentContainer.smiles.count("*")
        r_groups_n2 = 0

        side_chains = []
        for precursor in precursor_sucessors:
            side_chains_lst = \
                self.accessor.get_predecessors_by_ont_id_rel_type(precursor, "component_of")

            # try:
            precursor_parent = self.accessor.get_successors_by_ont_id_rel_type(precursor, "is_a")
            node_container = self.accessor.get_node_by_ont_id(precursor_parent[0])
            r_groups_n2 += node_container.smiles.count("*")
            # except:
            #     pass

            side_chains.extend(side_chains_lst)

        model_seed_ids = structuralParentContainer.aliases[ModelSEEDDatabases.MODEL_SEED.value]

        if r_groups_n1 == r_groups_n2 and all(ms_id not in self.exceptions_conf for ms_id in model_seed_ids):
            res = []

            for side_chain in side_chains:
                container = self.accessor.get_node_by_ont_id(side_chain)
                res.append(container)

            return res

        parent_neutralized_smiles, _ = chemo_utilities.NeutraliseCharges(structuralParentContainer.smiles)
        backup_parent = parent_neutralized_smiles
        parent_neutralized_smiles = parent_neutralized_smiles.replace("(*)", "").replace("*", "")

        parent_smarts = MolFromSmarts(parent_neutralized_smiles)
        if not parent_smarts:
            parent_smarts = MolFromSmarts(backup_parent)

        sidechain_compound = AllChem.DeleteSubstructs(mol, parent_smarts)
        sidechain2, sidechains_smiles_list2 = chemo_utilities.retrieve_fragments(sidechain_compound)
        components = []

        for sidechain in side_chains:
            container = self.accessor.get_node_by_ont_id(sidechain)

            neutralized_smiles, _ = chemo_utilities.NeutraliseCharges(container.smiles)

            for core in self.cores:
                sidechain1 = AllChem.DeleteSubstructs(MolFromSmiles(neutralized_smiles), MolFromSmiles(core))
                sidechain1, sidechains_smiles_list1 = chemo_utilities.retrieve_fragments(sidechain1)

                for side_chain_smiles in sidechains_smiles_list2:
                    if sidechains_smiles_list1:
                        equal_side_chains = \
                            chemo_utilities.check_if_molecules_are_equal(MolFromSmiles(side_chain_smiles),
                                                                         MolFromSmiles(sidechains_smiles_list1[0]))

                        if equal_side_chains:
                            if container not in components:
                                components.append(container)

        return components

    def __check_if_product_is_compatible_with_ont_container(self, model_seed_compounds, ont_ids):

        res = []

        for ont_id in ont_ids:

            for compound in model_seed_compounds:

                parent = self.accessor.get_successors_by_ont_id_rel_type(ont_id, "is_a")
                if parent:

                    ms_id = compound.getDbId()
                    container = self.accessor.get_node_by_ont_id(parent[0])
                    aliases = container.aliases
                    if ModelSEEDDatabases.MODEL_SEED.value in aliases.keys():
                        for alias in aliases[ModelSEEDDatabases.MODEL_SEED.value]:
                            if ms_id == alias:
                                res.append(ont_id)

        return res

    def __simplify_compounds(self, modelseed_compounds):
        res = []
        for compound in modelseed_compounds:
            if compound in self.compounds_simplifier.keys():
                res.append(self.compounds_simplifier.get(compound))
            else:
                res.append(compound)

        return res


if __name__ == "__main__":
    coiso = OntologyGeneratorDA()
    core_fatty_acids = \
        "C(=O)O"
    core_alcohol = "*CO"

    reaction_transformation = "[*:1]C(=O)OC[C@@H](O)COP(=O)(O)OC[C@@H](O)CO.[*:2]C(=O)OC[C@H](COP(=O)(O)OC[C@@H](O)CO)OC([*:3])=O>>[*:1]C(=O)OC[C@@H](COP(O)(=O)OC[C@@H](COC([*:2])=O)OC([*:3])=O)O"
    coiso.create_new_compounds_using_transformation(reaction_transformation, [110, 41], [765625], "LPLIPA6", "BiGG",
                                                    True, True, reverse=False)

    # coiso("TRIGLSYN-PWY", [core_fatty_acids,core_alcohol] , "transformations_TRIGLSYN-PWY", ["PWY-7277","PWY-7835"],["RXN-12383","RXN-1641"])

    # coiso.handle_reaction("RXN-15211","PWY-7277","transformations_Sphingolipid-Biosynthesis")
    # coiso.handle_reaction("RXN-15212", "PWY-7277", "transformations_Sphingolipid-Biosynthesis")

    # transformations_cardiolipin = TransformationsHandler()
    # transformations_cardiolipin.choose_transformations_from_pathway("TRIGLSYN-PWY", ["cpd27031", "cpd27029", "cpd00487"])
    # transformations_cardiolipin.save_transformations("transformations_TRIGLSYN-PWY")
