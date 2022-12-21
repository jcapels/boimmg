import re

from src.boimmgpy.kegg.kegg_compound import KeggCompound
from src.boimmgpy.kegg.kegg_pathway import KeggPathway
from src.boimmgpy.kegg.kegg_reaction import KeggReaction
from src.boimmgpy.service.network_handlers.metabolic_network import MetabolicNetwork
from src.boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from src.boimmgpy.model_seed.model_seed_compounds_database import ModelSeedCompoundsDB
from biocyc import biocyc
from src.boimmgpy.definitions import ROOT_DIR


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█'):
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
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end='', flush=True)
    # Print New Line on Complete
    if iteration == total:
        print()


class PathwayHandler(object):

    def __init__(self, compounds_database=None, compounds_converter=None):
        self.target_pathways_map = {}
        self.__pathways_map = {}
        self.taxID_map = {}
        self.__home_path__ = ROOT_DIR

        if not compounds_converter:
            self.__compoundsIdConverter = CompoundsIDConverter()
        else:
            self.__compoundsIdConverter = compounds_converter

        if not compounds_database:
            self.__compounds_database = ModelSeedCompoundsDB()
        else:
            self.__compounds_database = compounds_database

        biocyc.set_organism('meta')

    @property
    def target_pathways_map(self):
        return self.__target_pathways_map

    @target_pathways_map.setter
    def target_pathways_map(self, target_pathways_map):
        self.__target_pathways_map = target_pathways_map

    def get_targets(self, pathwayid):
        return self.target_pathways_map[pathwayid]

    def save_in_file(self, filename):
        with open(filename, "w") as file:
            for pathway in self.target_pathways_map:
                lst = self.target_pathways_map[pathway]
                new_lst = []
                temp_lst = []
                for elem in lst:
                    ms_number = elem.split("cpd")[1]
                    if int(ms_number) > 100:
                        new_lst.append(elem)

                if not new_lst:
                    new_lst = temp_lst.copy()

                if new_lst:
                    file.write(">" + ",".join(new_lst) + "\n")
                    path = self.get_path_from_pathway(pathway)

                    file.write(pathway + ": " + ">".join(path) + "\n")

    def load_from_file(self, filename):
        with open(filename, "r") as file:
            line = file.readline()
            while line:
                targets = re.sub(">", "", line).strip()
                targets = targets.split(",")
                line = file.readline()
                while not re.search("^>", line) and line:

                    pathwayid_path = line.split(": ")
                    path = pathwayid_path[1].strip().split(">")
                    self.target_pathways_map[pathwayid_path[0]] = targets
                    network = MetabolicNetwork()
                    if len(path) == 1:
                        network.add_vertex(path[0])
                    else:
                        for i in range(0, len(path) - 1):
                            network.add_edge(path[i], path[i + 1])

                    self.__pathways_map[pathwayid_path[0]] = network
                    line = file.readline()

    def get_networks_from_target(self, target_modelseedid):
        return self.target_pathways_map[target_modelseedid]

    def get_target_pathways(self):
        return self.target_pathways_map

    def get_network_from_target_and_pathway_id(self, target_modelseedid, metacyc_pathway_id):
        return self.target_pathways_map.get(target_modelseedid).get(metacyc_pathway_id)

    def add_metacyc_pathway(self, metacyc_pathway_id, generic=False):
        biocyc.set_organism('meta')
        ptw = biocyc.get(metacyc_pathway_id)

        network = self.scrap_metacyc_path(ptw, generic)
        self.convert_compounds_into_modelseedid(network)

        pathway = network.get_pathway()
        if pathway:
            target_reaction = pathway[-1]
            last_reaction = biocyc.get(target_reaction)
            products = last_reaction.compounds_right
            targets = []

            for product in products:
                modelseedids = self.__compoundsIdConverter.convert_dbID_into_modelSeedId("MetaCyc", product.id)
                if modelseedids is not None:
                    targets.extend(modelseedids)

            self.__pathways_map[ptw.id] = network
            self.target_pathways_map[ptw.id] = targets

    def get_pathways(self):
        return self.__pathways_map

    def get_network_by_metacyc_pathway_id(self, pathway_id):
        return self.__pathways_map[pathway_id]

    def scrap_metacyc_path(self, ptw, generic):
        if generic:
            return self.__scrap_generic_metacyc_path(ptw)
        else:
            return self.__scrap_metacyc_path(ptw)

    def scrap_kegg_pathway(self, kegg_pathway):

        pathway = KeggPathway(kegg_pathway, False)
        return self.__scrap_generic_kegg_path(pathway)

    def __scrap_generic_kegg_path(self, pathway):

        network = MetabolicNetwork()

        for reaction in pathway.get_reactions():
            reaction = KeggReaction(reaction)
            reaction_id = reaction.id

            reactants = reaction.get_reactants()
            products = reaction.get_products()

            for reactant in reactants:
                reactant = KeggCompound(reactant)
                ids = self.__compoundsIdConverter.convert_dbID_into_modelSeedId("KEGG", reactant.id)
                for reactant_id in ids:
                    compound_container = self.__compounds_database.get_compound_by_id(reactant_id)
                    if "R" in compound_container.getFormula():
                        network.add_vertex(reaction_id)
                        network.node_types["reaction"].append(reaction_id)

                        network.node_types["metabolite"].append(reactant.id)
                        network.add_vertex(reactant.id)
                        network.add_edge(reactant.id, reaction_id)
                        network.add_edge(reaction_id, reactant.id)
                        break

            for product in products:
                product = KeggCompound(product)
                ids = self.__compoundsIdConverter.convert_dbID_into_modelSeedId("KEGG", product.id)
                for product_id in ids:
                    compound_container = self.__compounds_database.get_compound_by_id(product_id)
                    if "R" in compound_container.getFormula():
                        network.node_types["metabolite"].append(product.id)
                        network.add_vertex(product.id)
                        network.add_edge(reaction_id, product.id)
                        network.add_edge(product.id, reaction_id)
                        break

        metabolite_network = MetabolicNetwork()
        metabolite_network.convert_reaction_graph(network)
        return metabolite_network

    @staticmethod
    def __get_reactants_and_products_by_r_groups(reaction):

        products = reaction.get_products()
        reactants = reaction.get_reactants()

        products_r_groups = 0
        reactants_r_groups = 0

        for product in products:
            product = KeggCompound(product)
            formula = product.get_formula()
            if "R" in formula:
                match = re.match(formula, "R")
                start = match.start()
                number = formula[start + 1]

                if re.match(number, "[2-9]"):
                    number = int(number)
                    products_r_groups += number
                else:
                    products_r_groups += 1

        for reactant in reactants:
            reactant = KeggCompound(reactant)
            formula = reactant.get_formula()
            if "R" in formula:
                match = re.match(formula, "R")
                start = match.start()
                number = formula[start + 1]

                if re.match(number, "[2-9]"):
                    number = int(number)
                    reactants_r_groups += number
                else:
                    reactants_r_groups += 1

        if reactants_r_groups >= products_r_groups:

            return reactants, products

        else:
            return products, reactants

    @staticmethod
    def __scrap_metacyc_path(ptw):
        network = MetabolicNetwork()

        for reaction in ptw.reactions:
            reaction_id = reaction.id
            network.add_vertex(reaction_id)
            network.node_types["reaction"].append(reaction_id)

            reactants = reaction.compounds_left
            products = reaction.compounds_right

            for reactant in reactants:
                if reactant.id != "PROTON":
                    network.node_types["metabolite"].append(reactant.id)
                    network.add_vertex(reactant.id)
                    network.add_edge(reactant.id, reaction_id)

            for product in products:
                if product.id != "PROTON":
                    network.node_types["metabolite"].append(product.id)
                    network.add_vertex(product.id)
                    network.add_edge(reaction_id, product.id)

        metabolite_network = MetabolicNetwork()
        metabolite_network.convert_reaction_graph(network)
        return metabolite_network

    def __scrap_generic_metacyc_path(self, ptw):
        network = MetabolicNetwork()

        for reaction in ptw.reactions:
            reaction_id = reaction.id

            network.add_vertex(reaction_id)
            network.node_types["reaction"].append(reaction_id)

            reactants = reaction.compounds_left
            products = reaction.compounds_right

            if reaction_id == "RXN-17728":
                reactants = reaction.compounds_right
                products = reaction.compounds_left

            for reactant in reactants:
                ids = self.__compoundsIdConverter.convert_dbID_into_modelSeedId("MetaCyc", reactant.id)
                for reactant_id in ids:
                    compound_container = self.__compounds_database.get_compound_by_id(reactant_id)
                    if "R" in compound_container.getFormula():
                        network.node_types["metabolite"].append(reactant.id)
                        network.add_vertex(reactant.id)
                        network.add_edge(reactant.id, reaction_id)
                        break

            for product in products:
                ids = self.__compoundsIdConverter.convert_dbID_into_modelSeedId("MetaCyc", product.id)
                for product_id in ids:
                    compound_container = self.__compounds_database.get_compound_by_id(product_id)
                    if "R" in compound_container.getFormula():
                        network.node_types["metabolite"].append(product.id)
                        network.add_vertex(product.id)
                        network.add_edge(reaction_id, product.id)

        metabolite_network = MetabolicNetwork()
        metabolite_network.convert_reaction_graph(network)
        return metabolite_network

    def convert_compounds_into_modelseedid(self, network, db="MetaCyc"):
        res = {}
        for compound in network.get_nodes():
            modelseedids = self.__compoundsIdConverter.convert_dbID_into_modelSeedId(db, compound)
            if modelseedids is not None:
                res[compound] = modelseedids
        return res

    def get_path_from_pathway(self, pathway_id):
        network = self.__pathways_map[pathway_id]
        return network.get_complete_pathway()

    @staticmethod
    def get_unique(lst):
        res = []
        for elem in lst:
            if elem not in res:
                res.append(elem)

        return res

    def get_pathway_by_target(self, target):
        pathways = []
        for path in self.target_pathways_map:
            if target in self.target_pathways_map[path]:
                pathways.append(path)
        return pathways

    def get_all_paths_using_pathway_id(self, pathway_id, accept_cycles=False):
        network = self.get_network_by_metacyc_pathway_id(pathway_id)
        if not accept_cycles:
            network.prune_redundant_cycles()
        pathways = network.get_all_target_pathways()
        return pathways

    def get_path_model_seed_ids(self, ptw):
        network = self.__pathways_map[ptw]
        dict_res = self.convert_compounds_into_modelseedid(network)
        return dict_res
