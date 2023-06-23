class ModelSeedCompound:
    """ 
    Class representing a Model Seed compound.
    """
    def __init__(self, id:str, name:str, charge:int, formula:str):
        """
        Initializes a ModelSeedCompound instance.

        :param id: The ID of the compound.
        :type id: str
        :param name: The name of the compound.
        :type name: str
        :param charge: The charge of the compound.
        :type charge: int
        :param formula: The chemical formula of the compound.
        :type formula: str
        """
        self.__id = id
        self.__name = name
        self.__charge = charge
        self.__formula = formula
        self.__inchi = ""

    def setStructures(self, smiles:str=None, inchikey:str=None, inchi:str=None):
        """
        Sets the structures of the compound.

        :param smiles: The SMILES representation of the compound.
        :type smiles: str, optional
        :param inchikey: The InChI key of the compound.
        :type inchikey: str, optional
        :param inchi: The InChI representation of the compound.
        :type inchi: str, optional
        """
        if smiles is not None:
            self.__smiles = smiles
        if inchikey is not None:
            self.__inchikey = inchikey

        if inchi:
            self.__inchi = inchi

    def getDbId(self):
        """
        Returns the ID of the compound.

        :return: The ID of the compound.
        :rtype: str
        """
        return self.__id

    def getName(self):
        """
        Returns the name of the compound.

        :return: The name of the compound.
        :rtype: str
        """
        return self.__name

    def getSmiles(self):
        """
        Returns the smiles of the compound.

        :return: The smiles of the compound.
        :rtype: str
        """
        return self.__smiles

    def getInchikey(self):
        """
        Returns the inchikey of the compound.

        :return: The inchikey of the compound.
        :rtype: str
        """
        return self.__inchikey

    def getFormula(self):
        """
        Returns the formula of the compound.

        :return: The formula of the compound.
        :rtype: str
        """
        return self.__formula

    def getCharge(self):
        """
        Returns the charge of the compound.

        :return: The charge of the compound.
        :rtype: int
        """
        return self.__charge

    def getInchi(self):
        """
        Returns the inchi of the compound.

        :return: The inchi of the compound.
        :rtype: str
        """
        return self.__inchi
