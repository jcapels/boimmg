
from neo4j import GraphDatabase

from boimmgpy import definitions
from boimmgpy.database.containers.compound_node import CompoundNode
from boimmgpy.database.databases_babel import AliasesTransformer
from boimmgpy.database.interfaces.boimmg_database_accessor import BOIMMGDatabaseAccessor
from boimmgpy.utilities import file_utilities


class CompoundsDBAccessor(BOIMMGDatabaseAccessor):

    def __init__(self,uri="",user="",password=""):


        if not uri or not user or not password:
            uri, user, password = self.read_config_file()

        self.__uri = uri
        self.__user = user
        self.__password = password

    def read_config_file(self):
        configs = file_utilities.read_conf_file(definitions.BOIMMG_DATABASE)

        uri = configs["uri"]
        user = configs["user"]
        password = configs["password"]

        return uri, user, password

    def login(self):

        driver = GraphDatabase.driver(self.__uri, auth=(self.__user, self.__password), encrypted=False)
        self.tx = driver

    @property
    def tx(self):
        return self._tx

    @tx.setter
    def tx(self,tx):
        self._tx = tx


    def add_relationship(self,origin : int, target : int, type : str):
        self.login()
        with self.tx.session() as session:
            session.run("MATCH (c:Compound), (p:Compound)"
                        "WHERE c.boimmg_id = $target and p.boimmg_id = $origin "
                        "MERGE (c)<-[r:" + type+ "]-(p) "
                        "ON CREATE SET r.timestamp = timestamp() ",
                        origin = origin, target= target, type = type)

    def add_biosynthetic_relationship(self,origin : int, target : int, **kwargs):
        self.login()
        properties = []
        for key in kwargs:
            properties.append( key + ":'"+ kwargs.get(key)+"'")

        properties_cypher = ",".join(properties)

        with self.tx.session() as session:
            session.run("MATCH (c:Compound), (p:Compound)"
                    "WHERE c.boimmg_id = $target and p.boimmg_id = $origin "
                    "MERGE (c)<-[r:precursor_of {" + properties_cypher + "} ]-(p) "
                    "ON CREATE SET r.timestamp = timestamp()", origin = origin,target=target
                    )



    def get_predecessors_by_swiss_lipids_id(self, db_id) -> list:
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r]-(p:Compound) "
                        "WHERE c.swiss_lipids_id = $db_id "
                        "RETURN p, p.boimmg_id as id",
                        db_id=db_id)


            data = result.data()
            for node in data:
                node_properties = node.get("p")
                node_id = node.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id,node_properties,other_aliases)
                predecessors.append(node)

        return predecessors

    def get_predecessors_by_model_seed_id(self, db_id:str) -> list:
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r]-(p:Compound) "
                        "WHERE c.model_seed_id = $db_id "
                        "RETURN p, p.boimmg_id as id",
                        db_id=db_id)


            data = result.data()
            for node in data:
                node_properties = node.get("p")
                node_id = node.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id,node_properties,other_aliases)
                predecessors.append(node)

        return predecessors

    def get_node_from_model_seed_id(self, model_seed_id:str) -> CompoundNode:
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:ModelSeedCompound)-[:is_db_link_of]->(d:Compound) "
                                 "USING INDEX c:ModelSeedCompound(model_seed_id) "
                                 "WHERE c.model_seed_id = $model_seed_id "
                                 "RETURN d, d.boimmg_id as id",
                                  model_seed_id=model_seed_id)

            data = result.single()
            if data:
                node_properties = data.get("d")
                node_id = data.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)

                return node

            else: return None

    def get_all_aliases(self, node_id):
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c)-[:is_db_link_of]->(d:Compound) "
                                 "WHERE d.boimmg_id = $node_id "
                                 "RETURN c",
                                 node_id=node_id)

            data = result.data()
            res=[]
            if data:
                for node in data:
                    node_properties = node.get("c")
                    res.append(node_properties)


            return res

    def get_node_id_from_model_seed_id(self, model_seed_id: str) -> int:
        self.login()
        with self.tx.session() as session:

            result = session.run("MATCH (c:ModelSeedCompound)-[:is_db_link_of]->(d:Compound) "
                                 "USING INDEX c:ModelSeedCompound(model_seed_id) "
                                 "WHERE c.model_seed_id = $model_seed_id "
                                 "RETURN d.boimmg_id as id",
                                 model_seed_id=model_seed_id)

            data = result.single()
            if data:
                node_id = data.get("id")

                return node_id

            else:
                return None

    def get_node_from_lipid_maps_id(self, lipid_maps_id: str) -> CompoundNode:
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound) "
                                 "WHERE c.lipidmaps_id = $lipid_maps_id "
                                 "RETURN c, c.boimmg_id as id",
                                  lipid_maps_id=lipid_maps_id)

            data = result.single()
            if data:
                node_properties = data.get("c")
                node_id = data.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)

                return node

            else: return None

    def get_node_by_alias(self, key, alias) -> int:
        self.login()

        boimmg_property = AliasesTransformer.convert_model_seed_format_into_boimmg(key)

        if boimmg_property:
            with self.tx.session() as session:
                # result = session.run("MATCH (c:Compound) "
                #                      "WHERE c." + boimmg_property + "=$alias "
                #                      "RETURN ID(c) as id",
                #                      alias=alias)

                result = session.run("MATCH (c:Compound {" + boimmg_property + ": $alias}) "
                                     "RETURN c.boimmg_id as id",
                                     alias=alias)

                data = result.single()

                if data:
                    # node_properties = data.get("c")
                    node_id = data.get("id")
                    # node = CompoundNode(node_id, node_properties)

                    return node_id

        return

    def get_conjugates(self,ont_id : int) -> list:
        """Get conjugated acids and bases predecessors using as parameter the ontology identifier"""

        res = []
        conjugated_bases = \
            self.get_predecessors_by_ont_id_rel_type(ont_id,"conjugated_base_of")

        conjugated_acids =  \
            self.get_predecessors_by_ont_id_rel_type(ont_id,"conjugated_acid_of")

        res.extend(conjugated_bases)
        res.extend(conjugated_acids)

        conjugated_bases = \
            self.get_successors_by_ont_id_rel_type(ont_id, "conjugated_base_of")

        conjugated_acids = \
            self.get_successors_by_ont_id_rel_type(ont_id, "conjugated_acid_of")

        res.extend(conjugated_bases)
        res.extend(conjugated_acids)

        return res



    def get_predecessors_by_ont_id(self, ont_id: int) -> list:
        """Get predecessors using as parameter the ontology identifier"""
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r]-(p:Compound) "
                                 "WHERE c.boimmg_id = $ont_id "
                                 "RETURN p.boimmg_id as id",
                                 ont_id=ont_id)

            data = result.data()
            for node in data:
                node_id = node.get("id")
                predecessors.append(node_id)

        return predecessors

    def get_predecessors_only_from_lipid_maps_and_ms(self,ont_id: int):
        self.login()
        predecessors = []

        with self.tx.session() as session:
            result = session.run(
                "match (c:Compound)<-[:is_a]-(d:Compound) where c.boimmg_id = $ont_id "
                "and not exists(d.swiss_lipids_id) and (not exists(d.annotated) or d.annotated) "
                "return d.boimmg_id as id",
                ont_id=ont_id)

            data = result.data()
            for node in data:
                node_id = node.get("id")
                predecessors.append(node_id)

        return predecessors


    def get_predecessors_with_same_component(self, ont_id: int, relationship_type : str) -> list:
        """Get predecessors with only one component using as parameter the database identifier and the relationship type"""
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("match (c:Compound)<-[:"+relationship_type+"]-(d:Compound)<-[:component_of]-(g) where c.boimmg_id = $ont_id "
                                "with collect(g) as components, d ,c "
                                 "where size(components) = 1 return d.boimmg_id as id",
                                 ont_id=ont_id,
                                 rel_type =relationship_type)

            data = result.data()
            if data:
                for node in data:
                    node_id = node.get("id")
                    predecessors.append(node_id)

        return predecessors

    def get_compounds_with_only_one_component(self, ont_id: int, components_ids_list: list) -> list:
        """Get predecessors with only one component using as parameter the database identifier and the relationship type"""
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("match (c:Compound)<-[:is_a]-(d:Compound)<-[:component_of]-(g) where c.boimmg_id = $ont_id "
                                "with collect(g) as components, d "
                                "where size(components) = 1 with d as result "
                                "match (result)<-[:component_of]-(t) where t.boimmg_id in $components_ids_list "
                                "return result.boimmg_id as id ",
                                 ont_id=ont_id,
                                 components_ids_list =components_ids_list)

            data = result.data()
            if data:
                for node in data:
                    node_id = node.get("id")
                    predecessors.append(node_id)

        return predecessors

    def get_all_model_seed_ids(self,ont_id):
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r:is_db_link_of]-(p:ModelSeedCompound) "
                                 "WHERE c.boimmg_id = $ont_id "
                                 "RETURN p.model_seed_id as model_seed_id",
                                 ont_id=ont_id)

            data = result.data()

            model_seed_ids = []
            if data:
                for node in data:
                    model_seed_id = node.get("model_seed_id")
                    model_seed_ids.append(model_seed_id)

        return model_seed_ids


    def get_predecessors_by_ont_id_rel_type(self, ont_id: int, relationship_type : str) -> list:
        """Get predecessors using as parameter the database identifier and the relationship type"""
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r]-(p:Compound) "
                                 "WHERE c.boimmg_id = $ont_id AND TYPE(r) = $rel_type "
                                 "RETURN p.boimmg_id as id",
                                 ont_id=ont_id,
                                 rel_type =relationship_type)

            data = result.data()
            if data:
                for node in data:
                    node_id = node.get("id")
                    predecessors.append(node_id)

        return predecessors

    def get_predecessors_by_ont_id_reaction_id(self, ont_id: int, reaction_id: str) -> list:
        """Get predecessors using as parameter the database identifier and the relationship type"""
        self.login()
        predecessors = []
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r:precursor_of]-(p:Compound) "
                                 "WHERE c.boimmg_id = $ont_id AND r.reaction = $reaction_id "
                                 "RETURN p.boimmg_id as id",
                                 ont_id=ont_id,
                                 reaction_id =reaction_id)

            data = result.data()
            if data:
                for node in data:
                    node_id = node.get("id")
                    predecessors.append(node_id)

        return predecessors

    def get_all_successors_by_ont_id_rel_type(self,ont_id: int, relationship_type : str,threshold=50) -> list:
        self.login()

        parent_container = self.get_node_by_ont_id(ont_id)
        l = [ont_id]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            successors = self.get_successors_by_ont_id_rel_type(node, relationship_type)
            if node != parent_container.id:
                res.append(node)
            for elem in successors:

                if elem not in res and elem not in l:
                    l.append(elem)
        return res

    def get_all_predecessors_by_ont_id_rel_type(self, ont_id: int, relationship_type : str) -> list:
        self.login()

        parent_container = self.get_node_by_ont_id(ont_id)
        l = [ont_id]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            predecessors = self.get_predecessors_by_ont_id_rel_type(node, relationship_type)
            if node != parent_container.id:
                res.append(node)
            for elem in predecessors:
                if elem not in res and elem not in l:
                    l.append(elem)
        return res


    def get_compound_by_inchikey(self,inchikey):
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound) "
                                 "WHERE c.inchikey contains $inchikey "
                                 "RETURN c, c.boimmg_id as id "
                                 ,
                                 inchikey=inchikey[:-1])
            data = result.single()
            if data:

                node_properties = data.get("c")
                node_id = data.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)
                return node

            else: return None

    def get_compound_by_smiles_and_abstraction(self,smiles,generic):
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound) "
                                 "WHERE c.smiles = $smiles and c.generic = $generic "
                                 "RETURN c, c.boimmg_id as id "
                                 ,
                                 smiles=smiles,
                                 generic = generic)

            data = result.single()
            if data:

                node_properties = data.get("c")
                node_id = data.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)
                return node

            else:
                return None

    def get_compound_by_smiles(self,smiles):
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound) "
                                 "WHERE c.smiles = $smiles "
                                 "RETURN c, c.boimmg_id as id "
                                 ,
                                 smiles = smiles)
            data = result.single()
            if data:

                node_properties = data.get("c")
                node_id = data.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)
                return node

            else:
                return None


    def create_compound(self,**kwargs) -> CompoundNode:

        name = kwargs.get("name")
        formula = kwargs.get("formula")
        smiles = kwargs.get("smiles")
        generic = bool(kwargs.get("generic"))
        charge = kwargs.get("charge")
        inchikey = kwargs.get("inchikey")
        annotated = kwargs.get("annotated")
        inchi = kwargs.get("inchi")
        aliases = kwargs.get("aliases")
        self.login()
        with self.tx.session() as session:
            result = session.run("match (c:Compound) return MAX(c.boimmg_id) as max_id")
            data = result.single()

            last_id = data.get("max_id")

            new_id = last_id +1

        if aliases:
            return self.__create_compound_with_aliases_boimmg_format(boimmg_id = new_id, name = name, formula = formula,
                                                                     smiles = smiles,
                                                                     generic = generic,
                                                                     charge = charge,
                                                                     inchikey = inchikey, annotated = annotated,
                                                                     inchi = inchi, aliases = aliases)

        else:
            self.login()
            with self.tx.session() as session:
                result = session.run("MERGE (c:Compound {name : $name,"
                            "formula : $formula, smiles: $smiles,"
                            "generic: $generic ,"
                            "charge: $charge, inchikey:$inchikey, "
                            "annotated:$annotated }) "
                            "on create set c.boimmg_id = $id "
                            "RETURN c, c.boimmg_id as id"
                            ,
                             name = name, formula = formula, smiles = smiles,
                            generic = generic, charge = charge, inchikey = inchikey,
                            annotated=annotated,id = new_id)


                data = result.single()
                node_properties = data.get("c")
                node_id = data.get("id")
                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)

                return node

    def get_node_by_ont_id(self, ont_id: int) -> CompoundNode:
        """Get node container by ontology identifier"""
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound) "
                                 "WHERE c.boimmg_id = $ont_id "
                                 "RETURN c, c.boimmg_id as id",
                                 ont_id=ont_id)

            data = result.single()
            node_properties = data.get("c")
            node_id = data.get("id")

            other_aliases = self.get_all_aliases(node_id)
            node = CompoundNode(node_id, node_properties,other_aliases)

        return node

    def get_compounds_with_specific_parent_within_set_of_components(self, parent, components):

        self.login()
        with self.tx.session() as session:
            result = session.run("match (e:Compound)<-[:is_a]-(c:Compound)<-[:component_of]-(d:Compound) "
                                 "with collect(d) as components,e,c "
                                 "where e.boimmg_id = $parent and "
                                 "all (x in components where x.boimmg_id in $components_par) "
                                 "return c.boimmg_id as id ",
                                 components_par=components,
                                 parent = parent
                                 )

            data = result.data()
            res = []
            for node in data:
                node_id = node.get("id")
                res.append(node_id)

        return res

    def get_compounds_with_specific_parent_set_of_components(self, parent, components):

        self.login()
        with self.tx.session() as session:
            result = session.run("match (e:Compound)<-[:is_a]-(c:Compound)<-[:component_of]-(d:Compound) "
                                 "with collect(id(d)) as components,e,c "
                                 "where e.boimmg_id = $parent and "
                                 "components = $components_par "
                                 "return c.boimmg_id as id ",
                                 components_par=components,
                                 parent=parent
                                 )

            data = result.data()
            res = []
            for node in data:
                node_id = node.get("id")
                res.append(node_id)

        return res

        # res=[]
        # for child in children:
        #     child_components = self.get_predecessors_by_ont_id_rel_type(child,"component_of")
        #
        #     inList = True
        #     if child_components:
        #         for component in child_components:
        #             if component not in components:
        #                 inList = False
        #                 break
        #
        #         if inList:
        #             res.append(child)
        # return res



    def get_leaves_from_ont_id(self,ont_id:int) -> list:

        parent_container = self.get_node_by_ont_id(ont_id)
        if parent_container.generic:
            l = [parent_container.id]
            res = []
            while len(l) > 0:
                node = l.pop(0)
                predecessors = self.get_predecessors_by_ont_id_rel_type(node, "is_a")
                if node != parent_container.id:
                    if not predecessors:
                        res.append(node)
                for elem in predecessors:
                    if elem not in res and elem not in l:
                        l.append(elem)
            return res
        else:
            return []

    def get_all_parents(self,leaf):

        l = [leaf]
        res = []
        while len(l) > 0:

            node = l.pop(0)
            successors = self.get_successors_by_ont_id_rel_type(node, "is_a")

            if node != leaf:
                res.append(node)

            for elem in successors:
                if elem not in res and elem not in l:
                    l.append(elem)

        return res

    def check_if_node_exist_by_model_seed_id(self,model_seed_id):
        self.login()
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound) "
                                 "WHERE c.model_seed_id = $model_seed_id "
                                 "RETURN c, c.boimmg_id as id",
                                 model_seed_id=model_seed_id)

            data = result.single()
            node = None
            if data:
                node_properties = data.get("c")
                node_id = data.get("id")

                other_aliases = self.get_all_aliases(node_id)
                node = CompoundNode(node_id, node_properties,other_aliases)

        return node

    def get_successors_by_ont_id_rel_type(self, ont_id:int, relationship_type:str) -> list:
        self.login()
        successors = []
        with self.tx.session() as session:

            result = session.run("MATCH (c:Compound)<-[r]-(p:Compound) "
                                 "WHERE p.boimmg_id = $ont_id AND TYPE(r) = $rel_type "
                                 "RETURN c.boimmg_id as id",
                                 ont_id=ont_id,
                                 rel_type=relationship_type)

            data = result.data()
            if data:
                for node in data:
                    node_id = node.get("id")
                    successors.append(node_id)
        return successors

    def get_successors_and_biosynthetic_relationship_info_by_ont_id_rel_type(self, ont_id:int, relationship_type:str):
        self.login()
        successors = []
        with self.tx.session() as session:
            result = session.run("MATCH (c:Compound)<-[r]-(p:Compound) "
                                 "WHERE p.boimmg_id = $ont_id AND TYPE(r) = $rel_type "
                                 "return r.pathway as pathway, r.reaction as reaction, c.boimmg_id as id ",
                                 ont_id = ont_id, rel_type=relationship_type)

            data = result.data()
            if data:
                for node in data:
                    node_id = node.get("id")
                    pathway = node.get("pathway")
                    reaction = node.get("reaction")
                    successors.append((node_id,pathway,reaction))

        return successors



    def establish_structural_relationship(self,origin,target):
        self.login()

        with self.tx.session() as session:
            session.run("MATCH (c:Compound),(p:Compound) " 
                        "where c.boimmg_id = $target and p.boimmg_id = $origin "
                        "MERGE (c)<-[r:is_a]-(p) "
                         "on create set r.time_stamp = TIMESTAMP()",
                        target = target, origin = origin)

    def __create_compound_with_aliases_boimmg_format(self,**kwargs):

        name = kwargs.get("name")
        formula = kwargs.get("formula")
        smiles = kwargs.get("smiles")
        generic = kwargs.get("generic")
        charge = kwargs.get("charge")
        inchikey = kwargs.get("inchikey")
        annotated = kwargs.get("annotated")
        inchi = kwargs.get("inchi")
        aliases = kwargs.get("aliases")
        boimmg_id = kwargs.get("boimmg_id")


        self.login()
        with self.tx.session() as session:
            result = session.run("MERGE (c:Compound {name : $name, inchikey:$inchikey}) "
                                 "ON create set "
                                 "c.boimmg_id = $boimmg_id, "
                                 "c.formula = $formula, c.smiles = $smiles,"
                                 "c.generic = $generic ,"
                                 "c.charge = $charge,"
                                 "c.annotated = $annotated , c.inchi = $inchi "
                                 "RETURN c, c.boimmg_id as id"
                                 ,
                                 name=name, formula=formula, smiles=smiles,
                                 generic=generic, charge=charge, inchikey=inchikey,
                                 annotated=annotated, inchi = inchi, boimmg_id=boimmg_id)

            data = result.single()
            node_properties = data.get("c")
            node_id = data.get("id")
            for alias in aliases:
                session.run("MATCH (c:Compound) "
                            "where c.boimmg_id = $id "
                            "SET c."+alias+" = $alias_value ",
                            alias_value = aliases[alias],
                            id = node_id)

            other_aliases = self.get_all_aliases(node_id)
            node = CompoundNode(node_id, node_properties,other_aliases)

            return node




