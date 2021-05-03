import math
import re
from collections import defaultdict

from cobra import Model, Reaction
from sympy import solve, Symbol, Integer

from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.utilities import file_utilities
from boimmgpy.utilities.annotation_utils import AnnotationUtils
from boimmgpy.definitions import COMPOUNDS_ANNOTATION_CONFIGS_PATH


class CompoundsRevisor:

    def __init__(self,model,universal_model = None,compoundsAnnotationConfigs = None):
        if not universal_model:
            self.__universal_model = Model()
        else:
            self.__universal_model=universal_model

        if not compoundsAnnotationConfigs:
            self.compoundsAnnotationConfigs = file_utilities.read_conf_file(COMPOUNDS_ANNOTATION_CONFIGS_PATH)

        else:
            self.compoundsAnnotationConfigs = compoundsAnnotationConfigs

        self.__model = model
        self.__get_hydrogen_from_model()
        self.__changed_reaction = None

    def __lp_to_balance_reaction(self,eq):

        Ls = list('abcdefghijklmnopqrstuvwxyz')

        Ss, Os, Es, a, i = defaultdict(list), Ls[:], [], 1, 1

        for p in eq.split('->'):
            for k in p.split('+'):
                c = [Ls.pop(0), 1]
                for e, m in re.findall('([A-Z][a-z]?)([0-9]*)', k):
                    m = 1 if m == '' else int(m)
                    a *= m
                    d = [c[0], c[1] * m * i]
                    Ss[e][:0], Es[:0] = [d], [[e, d]]
            i = -1


        Ys = dict((s, eval('Symbol("' + s + '")')) for s in Os if s not in Ls)

        Qs = [eval('+'.join('%d*%s' % (c[1], c[0]) for c in Ss[s]), {}, Ys) for s in Ss] + [Ys['a'] - a]

        k = solve(Qs, *Ys)

        # k = lcm(Qs, *Ys)
        new_k = k.copy()
        defined_variable_value = 0
        symbols_list = []
        list_of_equivalent_symbols = []
        for symbol in k:

            if isinstance(k[symbol], Symbol):
                new_k[k[symbol]] = symbol
                list_of_equivalent_symbols.append(symbol)

            t = k[symbol]
            present_symbols = t.free_symbols
            for s in present_symbols:
                if s not in symbols_list:
                    symbols_list.append(s)

            if isinstance(k[symbol], Integer):
                defined_variable_value = k[symbol]

        if len(symbols_list) == 1 and defined_variable_value != 0:

            to_process_symbol = symbols_list[0]
            for symbol in new_k:

                if not isinstance(new_k[symbol], Integer):
                    t = new_k[symbol]

                    new = t.subs({to_process_symbol: defined_variable_value})

                    new_k[symbol] = new

                if new_k[symbol] in list_of_equivalent_symbols:
                    t = new_k[symbol]

                    new = t.subs({new_k[symbol]: defined_variable_value})

                    new_k[symbol] = new

        return new_k, Ys

    def balance_reaction(self,reaction: Reaction):

        eq = self.__get_reaction_equation(reaction)

        try:
            k,Ys = self.__lp_to_balance_reaction(eq)
        except:
            return False

        isHydrogenInReaction = True

        if len(reaction.compartments) == 1:
            hydrogen = self.__hydrogens_in_model.get(list(reaction.compartments)[0])
            isHydrogenInReaction = self.__check_if_compound_in_reaction(reaction, hydrogen)

        if len(k) == len(Ys) and 0 not in k.values():

            try:
                stoich = self.__process_result_of_new_balanced_reaction(k,Ys,eq,reaction)
                reaction.add_metabolites(stoich,combine=False)
                return True

            except:
                return False
            # print('->'.join('+'.join(pM(N.pop(0)) + str(t) for t in p.split('+')) for p in eq.split('->')))

        elif not isHydrogenInReaction:

            reactants_products = eq.split("->")
            reactants_products[0] += "+H"

            new_eq = "->".join(reactants_products)
            k, Ys = self.__lp_to_balance_reaction(new_eq)

            if len(k) == len(Ys) and 0 not in k.values():

                new_reaction = reaction.copy()

                new_reaction.add_metabolites({hydrogen: -1})

                try:
                    stoich = self.__process_result_of_new_balanced_reaction(k, Ys, new_eq, new_reaction)

                    reaction.add_metabolites(stoich, combine=False)
                    return True

                except:
                    return False
            else:
                reactants_products = eq.split("->")
                reactants_products[1] += "+H"

                new_eq = "->".join(reactants_products)
                k, Ys = self.__lp_to_balance_reaction(new_eq)

                if len(k) == len(Ys) and 0 not in k.values():

                    new_reaction = reaction.copy()

                    new_reaction.add_metabolites({hydrogen: -1})

                    try:
                        stoich = self.__process_result_of_new_balanced_reaction(k, Ys, new_eq, new_reaction)

                        reaction.add_metabolites(stoich, combine=False)
                        return True

                    except:
                        return False

        elif len(reaction.compartments) == 1:

            reactants_products = eq.split("->")

            reactants = reactants_products[0].split("+")
            products = reactants_products[1].split("+")

            if "H" in reactants:
                reactants.remove("H")

                new_reactants = "+".join(reactants)
                new_eq = "->".join([new_reactants,reactants_products[1]])

                k, Ys = self.__lp_to_balance_reaction(new_eq)

                if len(k) == len(Ys):

                    new_reaction = reaction.copy()

                    new_reaction.subtract_metabolites({hydrogen: -1})
                    try:
                        stoich = self.__process_result_of_new_balanced_reaction(k, Ys, new_eq, new_reaction)

                        reaction.subtract_metabolites({hydrogen: -1})
                        reaction.add_metabolites(stoich, combine=False)
                        return True

                    except:

                        return False

            else:
                products.remove("H")

                new_products = "+".join(products)
                new_eq = "->".join([reactants_products[0], new_products])

                k, Ys = self.__lp_to_balance_reaction(new_eq)

                if len(k) == len(Ys):
                    new_reaction = reaction.copy()

                    new_reaction.subtract_metabolites({hydrogen: 1})
                    try:
                        stoich = self.__process_result_of_new_balanced_reaction(k, Ys, new_eq, new_reaction)

                        reaction.subtract_metabolites({hydrogen: 1})
                        reaction.add_metabolites(stoich, combine=False)
                        return True
                    except:
                        return False


        return False


    def __process_result_of_new_balanced_reaction(self,k,Ys,eq,reaction):

        N = []
        for s in sorted(Ys):
            value = k[Ys[s]]
            N.append(value)

        g = N[0]

        for a1, a2 in zip(N[0::2], N[1::2]):
            g = math.gcd(g, a2)

        N = [i / g for i in N]

        pM = lambda c: str(c) if c != 1 else ''

        res = {}

        equation_reactants_products = eq.split('->')
        reactants_in_equation = equation_reactants_products[0]
        products_in_equation = equation_reactants_products[1]

        for reactant in reactants_in_equation.split('+'):

            for reactant2 in reaction.reactants:
                if reactant2.formula == reactant:
                    stoich = pM(N.pop(0))
                    if stoich:
                        coef = int(stoich)
                    else:
                        coef = 1

                    res[reactant2] = coef * -1
                    break

        for product in products_in_equation.split('+'):

            for product2 in reaction.products:
                if product2.formula == product:
                    stoich = pM(N.pop(0))
                    if stoich:
                        coef = int(stoich)
                    else:
                        coef = 1

                    res[product2] = coef
                    break

        return res


    def __get_reaction_equation(self,reaction):
        reactants = reaction.reactants
        products = reaction.products
        eq = ""

        i=0
        while i<len(reactants)-1:

            reactant_formula = reactants[i].formula.replace("+","").replace("-","")
            eq+=reactant_formula+"+"

            i+=1


        eq += reactants[i].formula + "->"

        j=0
        while j < len(products) - 1:
            product_formula = products[j].formula.replace("+","").replace("-","")
            eq += product_formula + "+"

            j += 1

        eq += products[j].formula
        return eq


    def __balance_reaction(self, reaction, diff):
        """
        This method tries to balance a given reaction by using hydrogens

        :param Model.reaction reaction: reaction to balance.
        :param int diff: hydrogen difference between reactants and products

        :returns boolean: whether a reaction is balanced or not
        """

        compartment = list(reaction.compartments)
        if len(compartment) == 1:
            stoich_add = {}
            stoich_subtr = {}
            hydrogen = self.__hydrogens_in_model.get(compartment[0])
            new_reaction = reaction.copy()
            if diff < 0:

                isHydrogenInReaction = self.__check_if_compound_in_reaction(reaction,hydrogen)
                if isHydrogenInReaction:
                    stoich_subtr[hydrogen] = diff
                    new_reaction.subtract_metabolites(stoich_subtr)
                else:
                    stoich_add[hydrogen] = -diff
                    new_reaction.add_metabolites(stoich_add)

                self.__changed_reaction = new_reaction
                return True
            else:
                isHydrogenInReaction = self.__check_if_compound_in_reaction(reaction, hydrogen)
                if isHydrogenInReaction:
                    stoich_subtr[hydrogen] = -diff
                    new_reaction.subtract_metabolites(stoich_subtr)
                else:
                    stoich_add[hydrogen] = diff
                    new_reaction.add_metabolites(stoich_add)

                self.__changed_reaction = new_reaction
                return True


        return False



    def __check_if_compound_in_reaction(self,reaction,compound):
        for metabolite_in_reaction in reaction.metabolites:
            if metabolite_in_reaction.id == compound.id:
                return True
        return False

    def get_changed_reaction(self):
        return self.__changed_reaction

    def __get_hydrogen_from_model(self):
        """
        This method set searches for the hydrogen in each compartment of the model.
        The hydrogens are used to balance reactions.
        """

        self.__hydrogens_in_model = {}

        compartments = self.__model.compartments
        metabolites_in_model = self.__model.metabolites

        i = 0
        while i<len(metabolites_in_model):

            metabolite = metabolites_in_model[i]
            if metabolite.formula == "H" and metabolite.compartment in compartments:
                self.__hydrogens_in_model[metabolite.compartment] = metabolite

            i+=1

    def check_balanced_reaction(self, reaction):

        """
        This method checks if a given reaction is balanced or not, if not it tries to balance it.

        :param Model.reaction reaction: reaction to check whether is balanced.

        :returns boolean: whether a reaction is balanced or not
        """


        isBalanced = self.__check_if_reaction_is_balanced(reaction)

        if not isBalanced:

            return self.balance_reaction(reaction)

        else:
            return True


    def balance_and_clean_unbalanced_reactions(self):

        for reaction in self.__model.reactions:

            balanced = self.check_balanced_reaction(reaction)

            go = (not balanced and reaction not in self.__model.exchanges and
                  reaction not in self.__model.demands and reaction.annotation.get("sbo") != "SBO:0000629")
            if go:
                # self.__universal_model.add_reactions([reaction.copy()])
                self.__model.remove_reactions([reaction])


    def __check_if_reaction_is_balanced(self,reaction):
        (reactants,products) = self.__extract_elements_of_reactions_compounds(reaction)
        for element in reactants:
            if element in products.keys():
                diff = reactants[element] - products[element]
                if diff != 0:
                    return False
            else:
                return False

        for element in products:
            if element not in reactants.keys():
                return False

        return True


    def __extract_elements_of_reactions_compounds(self, reaction):
        """
       This method extracts the atomic elements of each compound in reaction

       :param Model.reaction reaction: reaction to extract atomic elements.

       :returns tuple: (reactants atomic elements, products atomic elements)
       """

        reactants_elements = {}
        products_elements = {}

        stoich = reaction.metabolites

        for metabolite in stoich.keys():
            formula = metabolite.formula
            coef = stoich[metabolite]
            match = re.findall("[A-Z][a-z]*[0-9]*", formula)
            for m in match:
                element = re.findall("[A-Z][a-z]*", m)[0]
                stoi = re.findall("[0-9]*", m)
                while "" in stoi:
                    stoi.remove("")

                if stoi:
                    if coef < 0:
                        reactants_elements[element[0]] = int(stoi[0]) * coef * -1 + reactants_elements.get(element[0], 0)
                    else:
                        products_elements[element[0]] = int(stoi[0]) * coef + products_elements.get(element[0], 0)
                else:
                    if coef < 0:
                        reactants_elements[element[0]] = coef * -1 + reactants_elements.get(element[0], 0)
                    else:
                        products_elements[element[0]] = coef + products_elements.get(element[0], 0)

        return (reactants_elements,products_elements)


    def check_compounds_representation_and_balance_reactions(self, reactions:list):

        res = []
        self.__get_hydrogen_from_model()

        compounds_ontology = CompoundsDBAccessor()

        to_remove = []
        for reaction in reactions:

            reaction_ontology_reactants = []
            reaction_ontology_products = []

            for compound in reaction.reactants:
                aliases = AnnotationUtils.get_annotation_from_cobra_annotation(compound)
                found = False
                for key in aliases:
                    for alias in aliases[key]:
                        new_key = self.compoundsAnnotationConfigs.get(key)

                        if key == "BOIMMG":
                            id = self.get_boimmg_id_from_annotation(aliases[key])
                            container = compounds_ontology.get_node_by_ont_id(id)

                            reaction_ontology_reactants.append(container)
                            found = True
                            break

                        elif new_key in self.__mapper.compounds_aliases_indexation.keys():
                            entities = self.__mapper.compounds_aliases_indexation.get(new_key).get(alias)
                            if entities:
                                for entity in entities:
                                    entity_id = entity.id

                                    if entity_id in self.__mapper.boimmg_db_model_map:

                                        ont_id = self.__mapper.boimmg_db_model_map.get(entity_id)
                                        entity_container = compounds_ontology.get_node_by_ont_id(ont_id)
                                        reaction_ontology_reactants.append(entity_container)
                                        found = True
                                        break
                    if found:
                        break


            for compound in reaction.products:

                found = False
                aliases = AnnotationUtils.get_annotation_from_cobra_annotation(compound)
                for key in aliases:
                    new_key = self.compoundsAnnotationConfigs.get(key)
                    for alias in aliases[key]:

                        if key == "BOIMMG":

                            id = self.get_boimmg_id_from_annotation(aliases[key])

                            container = compounds_ontology.get_node_by_ont_id(id)

                            reaction_ontology_products.append(container)
                            found = True
                            break


                        elif new_key in self.__mapper.compounds_aliases_indexation.keys():

                            entities = self.__mapper.compounds_aliases_indexation.get(new_key).get(alias)
                            if entities:
                                for entity in entities:
                                    entity_id = entity.id

                                    if entity_id in self.__mapper.boimmg_db_model_map:

                                        ont_id = self.__mapper.boimmg_db_model_map.get(entity_id)
                                        entity_container = compounds_ontology.get_node_by_ont_id(ont_id)
                                        reaction_ontology_products.append(entity_container)
                                        found = True
                                        break

                    if found:
                        break

            # self.__check_generic_complete_relationship(reaction_ontology_reactants,reaction_ontology_products,reaction)

            # if reaction_ontology_reactants and reaction_ontology_products:
            balanced = self.check_balanced_reaction(reaction)

            if not balanced:
                to_remove.append(reaction)

            else:
                res.append(reaction)

                # if generalization:
                #     self.__check_if_generic_complete_relationship()

        # print("reactions to remove: ")
        # print(to_remove)
        self.__model.remove_reactions(to_remove)
        return res

    def balance_reactions(self,reactions):
        to_remove = []
        res = []
        for reaction in reactions:
            balanced = self.check_balanced_reaction(reaction)

            if not balanced:
                to_remove.append(reaction)

            else:
                res.append(reaction)
        self.__model.remove_reactions(to_remove)
        return res

    def get_boimmg_id_from_annotation(self,boimmg_id):

        if isinstance(boimmg_id,list):
            boimmg_id = boimmg_id[0]

        boimmg_construction = self.compoundsAnnotationConfigs.get("BOIMMG_ID_CONSTRUCTION")

        id = int(boimmg_id.replace(boimmg_construction,""))

        return id

    def __check_generic_complete_relationship(self,reactants,products,reaction):

        generic_reactants = []
        complete_reactants = []
        generic_products = []
        complete_products = []

        compounds_ontology = CompoundsDBAccessor()

        for reactant in reactants:
            if reactant.generic:
                generic_reactants.append(reactant)
            else:
                complete_reactants.append(reactant)

        for product in products:
            if product.generic:
                generic_products.append(product)
            else:
                complete_products.append(product)

        if len(complete_reactants) == len(generic_products) or len(generic_reactants) == len(complete_products):
            return True

        elif not generic_reactants and not generic_products:
            reactants_components = []
            products_components = []

            reactants_components_dict = {}
            products_components_dict = {}


            for complete_reactant in complete_reactants:
                reactant_components = compounds_ontology.get_predecessors_by_ont_id_rel_type(complete_reactant.id,"component_of")
                reactants_components.append(reactant_components)
                reactants_components_dict[complete_reactant.id] = reactant_components


            for complete_product in complete_products:
                product_components = compounds_ontology.get_predecessors_by_ont_id_rel_type(complete_product.id,"component_of")
                products_components.append(product_components)
                products_components_dict[complete_product.id] = products_components

            reactants_components = sorted(reactants_components)
            products_components = sorted(products_components)

            if reactants_components!=products_components:
                deleted = self.__delete_components_from_reaction(reactants_components,
                                                       products_components,reactants_components_dict,
                                                       products_components_dict,reaction)

                if deleted:
                    return True

                else:
                    return False

        return True

    def __delete_components_from_reaction(self,reactants_components,
                                                products_components,reactants_components_dict,
                                                products_components_dict,reaction):

        to_delete = []

        if len(reactants_components) < len(products_components):
            reaction_products = reaction.products
            for product in reaction_products:
                product_components = products_components_dict[product.id]

                componentInReactants = []
                for component in product_components:
                    componentInReactants.append(component in reactants_components)

                if not any(componentInReactants):
                    for component in product_components:
                        products_components.remove(component)

                    to_delete.append(product)
                else:
                    # if it is necessary to delete a required component
                    return False

            if reactants_components != products_components:
                return False

            for metabolite in to_delete:
                reaction.subtract_metabolites({metabolite:1})

        else:
            reaction_reactants = reaction.reactants
            for reactant in reaction_reactants:
                reactant_components = reactants_components_dict[reactant.id]

                componentInProducts = []
                for component in reactant_components:
                    componentInProducts.append(component in reactants_components)

                if not any(componentInProducts):
                    for component in reactant_components:
                        reactants_components.remove(component)

                    to_delete.append(reactant)
                else:
                    # if it is necessary to delete a required component
                    return False

            if reactants_components != products_components:
                return False

            for metabolite in to_delete:
                reaction.subtract_metabolites({metabolite: 1})

        return True


    def set_model_mapper(self, mapper):
        self.__mapper = mapper
