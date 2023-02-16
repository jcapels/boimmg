import copy


class SplitMolecule(object):
    def __init__(self, complete_molecule, core, sidechains, ontid=None, dbid=None):
        self.__complete_molecule = complete_molecule
        self.__core__ = core
        self.__sidechains = sidechains
        self.__db_id = dbid
        self.__ontology_id__ = ontid

    def getSidechains(self):
        return self.__sidechains

    def getDbId(self):
        return self.__db_id

    def getOntID(self):
        return self.__ontology_id__

    def __copy__(self):
        return SplitMolecule(self.__complete_molecule, self.__core__,self.__sidechains,self.__ontology_id__,self.__db_id)

    def __deepcopy__(self, memo):
        return SplitMolecule(copy.deepcopy(self.__complete_molecule),
                             copy.deepcopy(self.__core__),copy.deepcopy(self.__sidechains),
                             copy.deepcopy(self.__ontology_id__),copy.deepcopy(self.__db_id))