from boimmgpy.database.databases_babel import BOIMMGDatabases, AliasesTransformer
from boimmgpy.database.boimmg_properties import BOIMMGProperties
from boimmgpy.database.interfaces.node import OntologyNode


class CompoundNode(OntologyNode):

    def __init__(self,node_id,node_properties,other_aliases):
        self.id = node_id
        self.extract_node_properties(node_properties,other_aliases)

    def __getattribute__(self, attr):
        try:
            return object.__getattribute__(self, attr)
        except AttributeError:
            return None

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self,id):
        self._id = id


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,name):
        self._name = name


    @property
    def inchi(self):
        return self._inchi

    @inchi.setter
    def inchi(self,inchi):
        self._inchi = inchi


    @property
    def inchikey(self):
        return self._inchikey

    @inchikey.setter
    def inchikey(self,inchikey):
        self._inchikey = inchikey

    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, smiles):
        self._smiles = smiles


    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, mass):
        self._mass = mass


    @property
    def generic(self):
        return self._generic

    @generic.setter
    def generic(self, generic):
        self._generic = generic


    @property
    def charge(self):
        return self._charge

    @property
    def aliases(self):
        return self._aliases

    @charge.setter
    def charge(self, value):
        self._charge = value


    @aliases.setter
    def aliases(self, aliases):
        self._aliases = aliases

    @property
    def formula(self):
        return self.__formula

    @formula.setter
    def formula(self,value):
        self.__formula = value

    @property
    def model_seed_id(self):
        return self.__model_seed_id

    @model_seed_id.setter
    def model_seed_id(self,value):
        self.__model_seed_id = value



    def extract_node_properties(self, node_properties,other_aliases):

        boimmg_databases = [db.value for db in BOIMMGDatabases]
        aliases = {}
        for property in node_properties:

            if property in boimmg_databases:
                aliases[property] = [node_properties.get(property)]


            if property == BOIMMGDatabases.MODEL_SEED.value:
                self.model_seed_id = node_properties[BOIMMGDatabases.MODEL_SEED.value]

        for aliases_node in other_aliases:
            for property in aliases_node:
                if property in aliases:
                    if aliases_node[property] not in aliases[property]:
                        aliases[property].append(aliases_node[property])
                else:
                    aliases[property] = [aliases_node[property]]

        if BOIMMGProperties.FORMULA.value in node_properties:
            self.formula = node_properties[BOIMMGProperties.FORMULA.value]

        if BOIMMGProperties.MASS.value in node_properties:
            self.mass = node_properties[BOIMMGProperties.MASS.value]

        if BOIMMGProperties.SMILES.value in node_properties:
            self.smiles = node_properties[BOIMMGProperties.SMILES.value]

        if BOIMMGProperties.INCHI.value in node_properties:
            self.inchi = node_properties[BOIMMGProperties.INCHI.value]

        if BOIMMGProperties.INCHIKEY.value in node_properties:
            self.inchikey = node_properties[BOIMMGProperties.INCHIKEY.value]

        if BOIMMGProperties.CHARGE.value in node_properties:
            self.charge = node_properties[BOIMMGProperties.CHARGE.value]

        self.generic = node_properties[BOIMMGProperties.GENERIC.value]

        if BOIMMGProperties.NAME.value in node_properties:
            self.name = node_properties[BOIMMGProperties.NAME.value]

        self.aliases = AliasesTransformer.transform_boimmg_aliases_into_model_seed(aliases)



