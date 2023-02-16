class ModelSeedCompound:

    def __init__(self, properties):

        self.__id = properties.get("modelseed_id")
        self.__name = properties.get("name")
        self.__smiles = properties.get("smiles")
        self.__inchikey = properties.get("inchikey")
        self.__formula = properties.get("formula")
        self.__charge = properties.get("charge")
        self.__inchi = properties.get("inchi")

    def setStructures(self, smiles=None, inchikey=None, inchi=None):
        if smiles is not None:
            self.__smiles = smiles
        if inchikey is not None:
            self.__inchikey = inchikey

        if inchi:
            self.__inchi = inchi

    def getDbId(self):
        return self.__id

    def getName(self):
        return self.__name

    def getSmiles(self):
        return self.__smiles

    def getInchikey(self):
        return self.__inchikey

    def getFormula(self):
        return self.__formula

    def getCharge(self):
        return self.__charge

    def getInchi(self):
        return self.__inchi
