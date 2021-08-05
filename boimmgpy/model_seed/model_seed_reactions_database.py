import re

from boimmgpy.model_seed.model_seed_reaction import ModelSeedReaction
from boimmgpy.utilities import file_utilities
from boimmgpy.definitions import ROOT_DIR, DATABASE_CONFIGS


class ModelSeedReactionsDB:

    def __init__(self):

        self.model_seed_reactions_database = {}
        self.__database_config = file_utilities.read_conf_file(DATABASE_CONFIGS)
        path = ROOT_DIR + self.__database_config["path_model_seed_db_reactions"]
        self.read_model_seed_reactions_database(path)

    def getReaction(self, reaction_id: str) -> ModelSeedReaction:
        return self.model_seed_reactions_database.get(reaction_id)

    def getDb(self):
        return self.model_seed_reactions_database

    def getReactionsProducingCompound(self, compound):
        reactions = self.getReactionsByCompound(compound)

        res = []
        for reaction in reactions:
            products = reaction.getProducts()

            if compound in products:
                res.append(reaction)

        return res

    def getReactionsByCompounds(self, reactants, products):
        res = []
        if reactants and products:
            for reaction in self.model_seed_reactions_database:
                reaction_container = self.model_seed_reactions_database[reaction]
                reaction_stoich = reaction_container.getStoichiometry()
                reaction_reactants = []
                reaction_products = []
                for comp in reaction_stoich:
                    if reaction_stoich[comp] < 0:
                        reaction_reactants.append(comp)
                    else:
                        reaction_products.append(comp)

                found_reactants = 0
                found_products = 0
                for reactant in reactants:
                    if reactant in reaction_reactants:
                        found_reactants += 1
                    else:
                        break

                for product in products:
                    if product in reaction_products:
                        found_products += 1

                if found_reactants == len(reactants) and found_products == len(products):
                    res.append(reaction)

            return res

    def getReactionsByCompound(self, compound):
        res = []
        for reaction in self.model_seed_reactions_database:
            reaction_container = self.model_seed_reactions_database[reaction]
            compounds = reaction_container.getCompounds()
            if compound in compounds:
                res.append(reaction_container)
        return res

    def read_model_seed_reactions_database(self, path):

        with open(path, "r") as reaction_file:
            reaction_file.readline()
            lines = reaction_file.readlines()
            for line in lines:
                line_lst = line.split("\t")

                is_transport = int(line_lst[5])
                is_obsolete = int(line_lst[18])
                if is_transport == 0 and is_obsolete == 0:

                    reaction_id = line_lst[0]
                    compounds_list = line_lst[16].split(";")
                    compounds = []
                    for compound in compounds_list:
                        compound = re.sub("\[.\]", "", compound)
                        compounds.append(compound)
                    if line_lst[12] != "null":
                        aliases = line_lst[12].split("|")
                        aliases_dict = {}
                        for alias in aliases:
                            clean_alias = re.sub("\s+", "", alias)
                            alias_lst = clean_alias.split(":")
                            aliases_dict[alias_lst[0]] = alias_lst[1].split(";")

                    stoichiometry = line_lst[4]
                    stoichiometry_lst = stoichiometry.split(";")
                    stoich_dict = {}
                    for compound in stoichiometry_lst:
                        if compound != "":
                            coeff_compound_id = compound.split(":")
                            coeff = float(coeff_compound_id[0])
                            compound_id = coeff_compound_id[1]
                            stoich_dict[compound_id] = coeff

                    pathways = line_lst[11]
                    pathways_dict = {}
                    if pathways != "null":
                        database_pathways_lst = pathways.split("|")

                        for database in database_pathways_lst:
                            database_name_and_path = database.split(":")
                            database_name = database_name_and_path[0].strip()
                            paths = database_name_and_path[1]
                            pathways_dict[database_name] = paths

                    equation = re.sub("\[[0-9]\]", "", line_lst[6])
                    name = line_lst[2]
                    direction = line_lst[9]
                    ec_numbers = line_lst[13].split(" | ")
                    deltag = float(line_lst[14])
                    pathways = pathways_dict
                    aliases = aliases_dict
                    stoichiometry = stoich_dict

                    self.model_seed_reactions_database[reaction_id] = \
                        ModelSeedReaction(reaction_id, name, equation, direction, ec_numbers, deltag, pathways,
                                          compounds, aliases, stoichiometry)
