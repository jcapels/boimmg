import re

from boimmgpy.kegg.kegg_entity import KeggEntity


class KeggReaction(KeggEntity):

    def __init__(self, entry, byModules=True):
        super().__init__(entry)
        self.compounds = {}
        self.products = {}
        self.reactants = {}
        self.byModules = byModules
        self.__get_all_information()

    def __get_all_information(self):
        info = self.get_raw_data()

        for text in info.split("\n"):
            new_line = re.sub('\s{2,}', "\t", text)

            if re.search("\ANAME", new_line):
                self.name = new_line.split("\t")[1]

            elif re.search("\AEQUATION", new_line):
                self.equation = new_line.split("\t")[1]
                self.get_reactants_and_products()

            elif re.search("\AENZYME", new_line):
                self.enzymes = new_line.split("\t")[1:]

            elif re.search("\ARCLASS", new_line):
                self.rclass = {}
                compounds = new_line.split("\t")[2]
                compound_pair = compounds.split("_")
                compound1 = compound_pair[0]
                compound2 = compound_pair[1]
                self.rclass[compound1] = compound2
                self.rclass[compound2] = compound1

    def get_name(self):
        try:
            name = self.name
            return name
        except:
            return None

    def get_reactants_and_products(self):
        match = re.search("<*=>*", self.equation)
        start = match.start()
        end = match.end()

        self.direction = self.equation[start:end]
        if self.direction == "<=>":
            self.reversibility = True
        else:
            self.reversibility = False

        sides = re.split("<*=>*", self.equation)
        reactants = sides[0].split("+")

        compound = None
        for reactant in reactants:
            coef = 1
            r = reactant.split(" ")
            for l in r:
                if re.search("^[0-9][0-9]*$", l):
                    coef = int(re.sub("\s+", "", l))
                elif re.search("C[0-9]*", l):
                    compound = re.sub("\s+", "", l)
            self.reactants[compound] = coef * -1

        products = sides[1].split("+")
        for product in products:
            r = product.split(" ")
            coef = 1
            for l in r:
                if re.search("^[0-9][0-9]*$", l):
                    coef = int(re.sub("\s+", "", l))
                elif re.search("C[0-9]*", l):
                    compound = re.sub("\s+", "", l)

            self.products[compound] = coef

        self.compounds.update(self.products)
        self.compounds.update(self.reactants)

    def get_id(self):
        return self.id

    def get_stoichiometry(self):
        return self.compounds

    def get_compounds(self):
        return self.compounds

    def get_reactants(self):
        return self.reactants

    def get_products(self):
        return self.products

    def get_reversibility(self):
        return self.reversibility

    def get_direction(self):
        return self.direction

    def get_equation(self):
        return self.equation

    def get_rclass(self):
        return self.rclass
