


class ModelSeedReaction:

    def __init__(self,reaction_id, name, equation, direction, ec_numbers, deltag, pathways, compounds, aliases, stoichiometry):
        self.id = reaction_id
        self.name = name
        self.equation = equation
        self.set_direction(direction)
        self.ec_numbers = ec_numbers
        self.deltag = deltag
        self.pathways = pathways
        self.compounds = compounds
        self.aliases = aliases
        self.__setReversibility()
        if "Name" in self.aliases.keys():
            del self.aliases["Name"]
        self.stoichiometry = stoichiometry

    def set_direction(self,direction):
        if direction == "=":
            self.direction = "<=>"
        elif direction == ">":
            self.direction = "=>"
        elif direction == "<":
            self.direction = "<="

    def getEcNumbers(self):
        return self.ec_numbers

    def getDirection(self):
        return self.direction

    def getDbId(self):
        return self.id

    def getReactants(self):
        reactants = []
        for compound in self.stoichiometry:
            if self.stoichiometry[compound]<0:
                reactants.append(compound)

        return reactants

    def getProducts(self):
        products = []
        for compound in self.stoichiometry:
            if self.stoichiometry[compound] > 0:
                products.append(compound)

        return products

    def getCompounds(self):
        return self.compounds

    def getDeltaG(self):
        return self.deltag

    def getName(self):
        return self.name

    def getStoichiometry(self):
        return self.stoichiometry

    def getAliases(self):
        return self.aliases

    def getAliasesByDatabase(self,db):
        if db == "ModelSEED":
            return [self.getDbId()]

        if db in self.aliases.keys():
            return self.aliases.get(db)
        else:
            return None

    def __setReversibility(self):
        if self.direction == "=":
            self.reversibility = True

        else:
            self.reversibility= False

    def getReversibility(self):
        return self.reversibility