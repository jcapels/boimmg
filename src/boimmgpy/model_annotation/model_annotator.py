import re
from tqdm import tqdm
from collections import defaultdict
from collections import Counter
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager
from joblib import Parallel, delayed
from typing import List
from cobra import Model,Metabolite


class LipidNameAnnotator:
    """
    Class that given a model implements the annotation of the structurally defined lipids in that model.
    """
    def __init__(self) -> None:
        self.backbone = None
        self.converted_lipid_class_dict = {}
        self.check_annotation_dict = {}
        self.results = {}
        self.counter = {}


    def login(self):
        """Method to create the connection to the database

        Returns:
            _type_: _description_
        """
        driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
        session = driver.session()
        return session

    def treat_data(self,dict_list:List[tuple]):
        """Method responsible for handling the dictionary fractions originating from the multiprocessing of the medel_lipids_finder method.
         Here all partitions are assigned to a corresponding final dictionary.

        Args:
            dict_list (List[tuple]): List of tuples containing 4 dictionaries in each tuple
        """
        
        for converted_lipid_class_dict, check_annotation_dict, counter, results in dict_list:
            if isinstance(converted_lipid_class_dict, defaultdict):
                converted_lipid_class_dict = Counter(converted_lipid_class_dict)
            if isinstance(counter, defaultdict):
                counter = Counter(counter)

            self.counter = {i: self.counter.get(i, 0) + counter.get(i, 0)
            for i in set(self.counter).union(counter)}
            self.converted_lipid_class_dict = {i: self.converted_lipid_class_dict.get(i, 0) + converted_lipid_class_dict.get(i, 0)
            for i in set(self.converted_lipid_class_dict).union(converted_lipid_class_dict)}
        
            self.check_annotation_dict.update(check_annotation_dict)
            self.results.update(results)


    def model_lipids_finder(self,model:Model):
        """
        Method responsible for multiprocessing the metabolites of the model, thus accelerating the whole annotation procedure.

        Args:
            model (Model): Cobrapy metabolic model

        Returns:
            _type_: A tuple with four dictionaries first one refers to the extent of the lipid classes present in the model, 
            the second one refers to the pre-presence of annotation on the lipids, 
            third for the extent of the classes that the algorithm was able to annotate 
            and finally a dictionary with the identifiers in the model for the lipids and their identifier in the Lipid Maps and/or Swiss Lipids database
        """
        n_iterations = len(model.metabolites)
        parallel_callback = Parallel(-1)
        resultados = parallel_callback(delayed(self._model_lipids_finder)(model.metabolites[i]) for i in tqdm(range(n_iterations)))
        self.treat_data(resultados)
        return (
            self.converted_lipid_class_dict,
            self.check_annotation_dict,
            self.counter,
            self.results,
        )

    def _model_lipids_finder(self,metabolite):
        """Method that searchs for lipid metabolites in model and finds their synonyms in the BOIMMG database

        Returns:
            _type_: A tuple with three dictionaries, first one with sorted information about lipid class present in the model, second one with the annotation information of lipids in the model (False or True for the presence of annotation) an the last
            with the classes that the algorithm can annotate.
        """

        matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
        metabolite_original_name = metabolite.name
        metabolite_name = metabolite.name
        backbone = None
        found = False
        lipid_class_dict = defaultdict(int)
        check_annotation_dict={}
        counter = {}
        results = {}
        side_chain = []
        for match in matches:
            found = True
            metabolite_name = metabolite_name.replace(
                match.string[match.start() : match.end()], ""
            )
            side_chain.append(match.string[match.start() : match.end()])

        if found:
            check_annotation_dict = self.annotation_checker(lipid=metabolite)
            backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", metabolite_name)
            if len(side_chain) != 1:
                for a in range(len(side_chain) - 1):
                    backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", backbone)
            lipid_class_dict[backbone] += 1
            #count += 1
            self.backbone = backbone
            counter,results = self.search_lipid_synonyms(metabolite,side_chain)

        sorted_lipid_class_dict = sorted(
            lipid_class_dict.items(), key=lambda x: x[1], reverse=True
        )
        converted_lipid_class_dict = dict(sorted_lipid_class_dict)
        return (
            converted_lipid_class_dict,
            check_annotation_dict,
            counter,
            results,
        )

    def annotation_checker(self, lipid):
        """Method that checks if a given lipid is annotated in the model

        Args:
            lipid (_type_): Given lipid metabolite from the model
        """
        check_annotation_dict = {}
        annotation = lipid.annotation.keys()
        annotated = False
        check_annotation_dict[lipid.id] = annotated
        if "slm" in annotation or "lipidmaps" in annotation:
            annotated = True
            check_annotation_dict[lipid.id] = annotated
        return check_annotation_dict

    def search_lipid_synonyms(self, metabolite:Metabolite,side_chain:List[str]):
        """Method that implements the screening in the database for accurate lipid structure.
        This method calls all StaticMethods to do the screening in the database.

        Args:
            metabolite (Metabolite): Given Metabolite of the cobra model
            side_chain (List[str]): List of the side chains referent to the given Metavolite

        Returns:
            _type_: tuple with two dictionaries, first one for the extend of the lipid classes that the algorithm caugth and second one for the identifiers to do the annotation
        """
        session = self.login()
        results={}
        counter = defaultdict(int)
        is_compound = False
        lipid_id = self.get_synonym_id(session,self.backbone, side_chain)
        if lipid_id != None and not lipid_id[2]:
            for backbone in lipid_id[0]:
                lipid_compound = self.get_coumpound(session,backbone, lipid_id[1], is_compound)
                if lipid_compound != None and lipid_compound[0] != 0:
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(
                        session,lipid_compound[0], lipid_compound[1]
                        )
                    )
                    if structurally_defined_lipids:
                        if metabolite.id in results:
                            results[metabolite.id] = (
                                results[metabolite.id]
                                + structurally_defined_lipids
                            )
                        else:
                            results[metabolite.id] = structurally_defined_lipids
                        counter[self.backbone] += 1

                if lipid_compound != None and lipid_compound[0] == 0:
                    backbone_compound = self.get_compound_from_synonym(session,backbone)
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(session,
                            backbone_compound, lipid_compound[1]
                        )
                    )
                    if structurally_defined_lipids:
                        if metabolite.id in results:
                            results[metabolite.id] = (
                                results[metabolite.id]
                                + structurally_defined_lipids
                            )
                        else:
                            results[metabolite.id] = structurally_defined_lipids
                        counter[self.backbone] += 1

        if lipid_id != None and lipid_id[2] == True:
            is_compound = True
            for backbone in lipid_id[0]:
                lipid_compound = self.get_coumpound(session,backbone, lipid_id[1], is_compound)
                if lipid_compound != None:
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(
                            session,lipid_compound[0], lipid_compound[1]
                        )
                    )
                    if structurally_defined_lipids:
                        if metabolite.id in results:
                            results[metabolite.id] = (
                                results[metabolite.id]
                                + structurally_defined_lipids
                            )
                        else:
                            results[metabolite.id] = structurally_defined_lipids
                        counter[self.backbone] += 1
        return counter,results
    
    @staticmethod
    def get_synonym_id(session,backbone: str, side_chain: list):
        """Method that searches for the synonyms of the backbone and side chains for each lipid in the model.
        In cases where the synonym for the backbone doesnt exist it will search for the generic compound with the same name as the backbone:

        Args:
            backbone (str): Backbone portion of the lipid name
            side_chain (list): List of side chains from the lipid name
        Returns:
            _type_: Tuple with two lists and a Boolean. First list refers to the backbone ID, second for the side_chains ID's. Lastly the Bolean is True when the backbone ID is from a compound.
        """
        sidechain_id = []
        backbone_id = []
        backbone = backbone.replace("'", "")
        backbone = backbone.replace(" ", "")
        backbone_node = session.run(
            "match (s:Synonym) where s.synonym='"
            + str(backbone.lower())
            + "' return ID(s) as boimmg_id"
        )
        node = backbone_node.data()
        compound = False

        for value in node:
            backbone_id.append(value.get("boimmg_id"))

        for synonym in side_chain:
            synonym_split = synonym.split(":")
            if synonym_split[1] == "0":
                synonym = "C" + synonym
            result = session.run(
                "match (s:Synonym) where s.synonym='"
                + str(synonym.lower())
                + "' return ID(s) as boimmg_id"
            )
            data = result.data()
            for value in data:
                sidechain_id.append(value.get("boimmg_id"))

        if len(backbone_id) != 0 and len(sidechain_id) != 0:
            return backbone_id, sidechain_id, compound

        if len(backbone_id) == 0 and len(sidechain_id) != 0:
            result = session.run(
                "match (c:Compound) where c.name='"
                + str(backbone)
                + "' return ID(c) as boimmg_id"
            )
            data = result.data()
            for value in data:
                backbone_id.append(value.get("boimmg_id"))
            if backbone != 0:
                compound = True
                return backbone_id, sidechain_id, compound

    @staticmethod
    def get_coumpound(session,backbone: str, sidechains: list, flag: bool):
        """Method that search for the the specific compounds attached for the side chain and backbone synonym's.
        In case that the ID from the backbone is already from a compound it will search only for the side chains compounds.

        Args:
            backbone (str): Backbone's ID
            sidechains (list): List of side chains ID's
            flag (bool): Bolean is True when the backbone ID is from a compound

        Returns:
            _type_: Tuple with backbone and side chains Compounds ID's
        """
        side_chains_coumpound = []
        backbone_coumpound = 0

        if flag == False:
            for v in sidechains:
                result = session.run(
                    "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$v and exists(t.boimmg_id) return id(t)",
                    v=v,
                )
                data = result.data()
                for value in data:
                    if value.get("id(t)") not in side_chains_coumpound:
                        side_chains_coumpound.append(value.get("id(t)"))
            node = session.run(
                "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$backbone and exists(t.boimmg_id) return t.boimmg_id as boimmg_id",
                backbone=backbone,
            )
            for value in node:
                backbone_coumpound = value.get("boimmg_id")
            if len(side_chains_coumpound) != 0:
                return backbone_coumpound, side_chains_coumpound

        if flag == True:
            for v in sidechains:
                result = session.run(
                    "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$v and exists(t.boimmg_id) return id(t)",
                    v=v,
                )
                data = result.data()
                for value in data:
                    side_chains_coumpound.append(value.get("id(t)"))
                return backbone, side_chains_coumpound

    @staticmethod
    def get_compound_from_synonym(session,synonym_ID: int):
        """Method to get the generic compound directly connected to the backbone synonnym.

        Args:
            synonym_ID (): ID of the backbone synonym

        Returns:
            _type_: Generic Compound with specific backbone synonym
        """
        node = session.run(
            "match(s:Synonym)-[:is_synonym_of]->(t:Compound) where id(s)=$backbone and exists(t.boimmg_id) return t.boimmg_id as boimmg_id",
            backbone=synonym_ID,
        )
        backbone_coumpound = None
        for value in node:
            backbone_coumpound = value.get("boimmg_id")
        return backbone_coumpound

    @staticmethod
    def get_compounds_with_specific_parent_set_of_components(
        session,parent: str, components: list
    ):
        """Get all non genericall compounds with the specific set of components

        Args:
            parent (str): Backbone compound id
            components (list): List of side chains compound ID's

        Returns:
            _type_: List with all the componds found to the specific pack of components
        """
        result = session.run(
            "match (e:Compound)<-[:is_a]-(c:Compound)<-[:component_of]-(d:Compound) "
            "with collect(id(d)) as components,e,c "
            "where e.boimmg_id = $parent and "
            "components = $components_par "
            "return c.boimmg_id as id ",
            components_par=components,
            parent=parent,
        )

        data = result.data()
        res = []
        for node in data:
            node_id = node.get("id")
            res.append(node_id)

        if len(res) != 0:
            return res
        

    def set_annotation():
        
        pass

