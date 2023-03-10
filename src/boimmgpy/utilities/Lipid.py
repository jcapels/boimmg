
class Lipid:

    def __init__(self ,db_id ,name ,sys_name ,exactmass ,formula ,inchi ,inchikey ,smiles ,aliases ,charge=0):
        self.__db_id = db_id
        self.__charge = charge
        self.__name = name
        self.__mass = exactmass
        self.__sys_name = sys_name
        self.__formula = formula
        self.__inchi = inchi
        self.__inchikey = inchikey
        self.__smiles = smiles
        self.__aliases = aliases


    def getDbId(self):
        return self.__db_id

    def getName(self):
        return self.__name

    def getSmiles(self):
        return self.__smiles

    def getInchiKey(self):
        return self.__inchikey

    def getFormula(self):
        return self.__formula

    def getCharge(self):
        return 0

    def getAliases(self):
        return self.__aliases

    def getInchi(self):
        return self.__inchi
