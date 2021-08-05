
from neo4j import GraphDatabase

from boimmgpy.database.interfaces import node


class ReactionsDBAccessor:

    def login(self):
        uri = "bolt://palsson.di.uminho.pt:6094"
        driver = GraphDatabase.driver(uri, auth=("neo4j", "Tuning999"), encrypted=False)
        self.tx = driver

    @property
    def tx(self):
        return self._tx

    @tx.setter
    def tx(self, tx):
        self._tx = tx

    def add_relationship(self, origin: int, target: int, type: str):
        self.login()
        with self.tx.session() as session:
            session.run("MATCH (c:Reaction), (p:Reaction)"
                        "WHERE ID(c) = $target and ID(p) = $origin "
                        "MERGE (c)<-[r:" + type + "]-(p) "
                                                  "ON CREATE SET r.timestamp = timestamp() ",
                        origin=origin, target=target, type=type)

    def get_predecessors_by_ont_id(self, ont_id: int) -> list:
        """Get predecessors using as parameter the database identifier"""
        pass

    def get_predecessors_by_ont_id_rel_type(self, ont_id: int, relationship_type : str) -> list:
        """Get predecessors using as parameter the database identifier and the relationship type"""

        raise NotImplementedError

    def get_node_by_ont_id(self, ont_id: int) -> node:
        """Get predecessors using as parameter the database identifier and the relationship type"""
        raise NotImplementedError


    def add_reaction(self,**kwargs):
        self.login()

        name = kwargs.get("name")

        onto_stoichiometry = kwargs.get("onto_stoichiometry")

        ms_stoich = kwargs.get("ms_stoich")

        direction = kwargs.get("direction")

        equation = self.generate_equation(onto_stoichiometry, ms_stoich, direction)

        ec_numbers = kwargs.get("ec_numbers") #list of ec numbers

        annotated = kwargs.get("annotated")

        deltaG = kwargs.get("deltaG")

        # aliases = kwargs.get("aliases")
        if "model_seed_id" in kwargs:
            model_seed_id = kwargs.get("model_seed_id")
        else:
            model_seed_id = None

        with self.tx.session() as session:

            if model_seed_id:
                result = session.run("merge (r:Reaction {equation: $equation, ec_numbers: $ec_numbers, "
                            "model_seed_id: $model_seed_id, annotated: $annotated}) "
                             "on create set r.name = $name, r.deltaG = $deltaG "
                            "return id(r) as r_id",
                            name = name, equation = equation, ec_numbers = ec_numbers,
                             annotated = annotated, deltaG = deltaG, model_seed_id = model_seed_id)
            else:
                result = session.run("merge (r:Reaction {equation: $equation, ec_numbers: $ec_numbers, "
                                     "annotated: $annotated}) "
                                     "on create set r.name = $name, r.deltaG = $deltaG "
                                     "return id(r) as r_id",
                                     name=name, equation=equation, ec_numbers=ec_numbers,
                                     annotated=annotated, deltaG=deltaG)

            result_reaction = result.single()
            reaction_id = result_reaction.get("r_id")

            for metabolite in onto_stoichiometry:
                coef = onto_stoichiometry[metabolite]

                if coef < 0:
                    role = "reactant"

                else:
                    role = "product"

                session.run("match (c:Compound), (r:Reaction) "
                            "where c.boimmg_id = $compound_id and id(r) = $reaction_id "
                            "merge (c)-[s:is_intervenient_in {coef: $coef, role : $role}]->(r) "
                            "on create set s.timestamp = timestamp()",
                            coef = coef, role = role, compound_id = metabolite,
                            reaction_id = reaction_id)

            for metabolite in ms_stoich:
                coef = ms_stoich[metabolite]

                if coef < 0:
                    role = "reactant"

                else:
                    role = "product"

                session.run("match (r:Reaction) "
                            "where id(r) = $reaction_id "
                            "merge (c:ModelSeedCompound {model_seed_id: $compound_id})-[s:is_intervenient_in {coef: $coef, role: $role}]->(r) "
                            "on create set s.timestamp = timestamp()",
                            coef = coef, role = role, compound_id = metabolite,
                            reaction_id = reaction_id)


            return reaction_id



    def generate_equation(self,onto_stoichiometry, ms_stoich, direction):

        left_side = []
        right_side = []

        for metabolite in onto_stoichiometry:

            coef = onto_stoichiometry[metabolite]

            if coef < 0:

                if coef < -1:

                    left_side.append(str(abs(coef)) + " C_BOIMMG_" + str(metabolite))

                else:

                    left_side.append("C_BOIMMG_" + str(metabolite))

            elif coef == 1:
                right_side.append("C_BOIMMG_" +str(metabolite))

            else:
                right_side.append(str(coef) + " C_BOIMMG_" + str(metabolite))

        for metabolite in ms_stoich:

            coef = ms_stoich[metabolite]

            if coef < 0:

                if coef < -1:

                    left_side.append(str(abs(coef)) + metabolite)

                else:

                    left_side.append(metabolite)

            elif coef == 1:
                right_side.append(metabolite)

            else:
                right_side.append(str(coef) + metabolite)

        left_side_str = "+".join(left_side)
        right_side_str = "+".join(right_side)

        equation = direction.join([left_side_str,right_side_str])

        return equation











