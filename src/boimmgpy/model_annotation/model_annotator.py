import re
from tqdm import tqdm
from neo4j import GraphDatabase
from collections import defaultdict
from collections import Counter
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager
from boimmgpy.model_annotation._utils import set_metabolite_annotation_in_model
from joblib import Parallel, delayed
from typing import List, Tuple, Dict
from cobra import Model, Metabolite


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

    def login(self) -> GraphDatabase.driver:
        """
        Method to create the connection to the database

        :return: session linkage to remote database
        :rtype: 
        """
        driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
        session = driver.session()
        return session

    def treat_data(self, dict_list: List[tuple]):
        """
        Method responsible for handling the dictionary fractions generated
        from the multiprocessing of the model_lipids_finder method.
        Here all partitions are assigned to a corresponding final dictionary.

        :param dict_list: List of tuples containing 4 dictionaries in each tuple
        :type dict_list: List[tuple]
        """
        for converted_lipid_class_dict, check_annotation_dict, counter, results in dict_list:
            if isinstance(converted_lipid_class_dict, defaultdict):
                converted_lipid_class_dict = Counter(converted_lipid_class_dict)
            if isinstance(counter, defaultdict):
                counter = Counter(counter)

            self.counter = {i: self.counter.get(i, 0) + counter.get(i, 0)
                            for i in set(self.counter).union(counter)}
            self.converted_lipid_class_dict = {
                i: self.converted_lipid_class_dict.get(i, 0) + converted_lipid_class_dict.get(i, 0)
                for i in set(self.converted_lipid_class_dict).union(converted_lipid_class_dict)}

            self.check_annotation_dict.update(check_annotation_dict)
            self.results.update(results)

    def find_model_lipids(self, model: Model, n_jobs: int = 6) -> Tuple[Dict, Dict, Dict, Dict, Model]:
        """
        Method responsible for multiprocessing the metabolites of the model,
        thus accelerating the whole annotation procedure.

        :param model: Cobrapy metabolic model
        :type model: Model
        :param n_jobs: Number of CPU cores to be used for multiprocessing, defaults to -1
        :type n_jobs: int, optional
        :return:  A tuple with four dictionaries first one refers to the extent
            of the lipid classes present in the model,
            the second one refers to the pre-presence of annotation on the lipids, 
            third for the extent of the classes that the algorithm was able to annotate 
            and finally a dictionary with the identifiers in the model for the lipids and their identifier
            in the Lipid Maps and/or Swiss Lipids database
        :rtype: Tuple(Dict)
        """
        n_iterations = len(model.metabolites)
        parallel_callback = Parallel(n_jobs)
        resultados = parallel_callback(
            delayed(self._find_model_lipids)(model.metabolites[i]) for i in tqdm(range(n_iterations)))
        self.treat_data(resultados)
        final_model = set_metabolite_annotation_in_model(self.login(), self.results, model)

        return (
            self.converted_lipid_class_dict,
            self.check_annotation_dict,
            self.counter,
            self.results,
            final_model
        )

    def _find_model_lipids(self, metabolite: Metabolite) -> Tuple[Dict, Dict, Dict, Dict]:
        """
        Method that searches for lipid metabolites in model and finds their synonyms in the BOIMMG database.
        Firstly lipids are separated from another model Metabolites by the regex
        [0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))* that searches for patterns similar to the lipids side chains
        representation distinct from all other metabolites
        The second Regex  *(\([\-a-zA-Z0-9/|, ]*\)) is used to delete the side
        chain part of the name to get only the backbone part of the lipid

        :param metabolite: Lipid metabolite from GSM model
        :type metabolite: Metabolite
        :return: A tuple with four dictionaries, first one with sorted information about lipid
            class present in the model, second one with the annotation information of lipids in the model
            (False or True for the presence of annotation), third one
            with the classes that the algorithm can annotate and the last one with metabolite ID's in the model as Keys and the annotation id as values.
        :rtype: Tuple(Dict)
        """

        matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
        metabolite_name = metabolite.name
        found = False
        lipid_class_dict = defaultdict(int)
        check_annotation_dict = {}
        counter = {}
        results = {}
        side_chain = []
        for match in matches:
            found = True
            metabolite_name = metabolite_name.replace(
                match.string[match.start(): match.end()], ""
            )
            side_chain.append(match.string[match.start(): match.end()])

        if found:
            check_annotation_dict = self.check_annotation(lipid=metabolite)
            backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", metabolite_name)
            if len(side_chain) != 1:
                for a in range(len(side_chain) - 1):
                    backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", backbone)
            lipid_class_dict[backbone] += 1
            self.backbone = backbone
            
            counter, results = self.search_lipid_synonyms(metabolite, side_chain)

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

    @staticmethod
    def check_annotation(lipid: Metabolite) -> Dict:
        """
        Method that checks if a given lipid is annotated in the model

        :param lipid: given lipid from the model
        :type lipid: Metabolite
        :return: dictionary with lipid model ID as key and bool for the annotation as values
        :rtype: Dict
        """
        check_annotation_dict = {}
        annotation = lipid.annotation.keys()
        annotated = False
        check_annotation_dict[lipid.id] = annotated
        if "slm" in annotation or "lipidmaps" in annotation:
            annotated = True
            check_annotation_dict[lipid.id] = annotated
        return check_annotation_dict

    def search_lipid_synonyms(self, metabolite: Metabolite, side_chain: List[str]) -> Tuple[Dict, Dict]:
        """
        Method that implements the screening in the database for accurate lipid structure.
        This method calls all StaticMethods to do the screening in the database.


        :param metabolite: Given Metabolite of the cobra model
        :type metabolite: Metabolite
        :param side_chain: List of the side chains referent to the given Metabolite
        :type side_chain: List[str]
        :return: tuple with two dictionaries, first one for the extend of the lipid classes that the algorithm caught
        and second one for the identifiers to do the annotation
        :rtype: Tuple[Dict, Dict]
        """
        results = {}
        counter = defaultdict(int)
        get_synonym = self.get_synonym_id(self.backbone, side_chain)
        if get_synonym is not None and not get_synonym[2]:
            backbone_id, sidechain_id, compound = get_synonym[0], get_synonym[1], get_synonym[2]
            if False not in sidechain_id:
                for backbone in backbone_id:
                    get_compound = self.get_compound(backbone, sidechain_id, compound)
                    if get_compound is not None and get_compound[0] != 0:
                        backbone_coumpound, side_chains_compound = get_compound[0], get_compound[1]
                        structurally_defined_lipids = (
                            self.get_compounds_with_specific_parent_set_of_components(
                                backbone_coumpound, side_chains_compound
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
                    if get_compound is not None and get_compound[0] == 0:
                        backbone_coumpound, side_chains_compound = get_compound[0], get_compound[1]
                        results, counter = self.get_id_from_compound_with_no_backbone(backbone, side_chains_compound,
                                                                                    metabolite, results, counter)

        if get_synonym is not None and get_synonym[2]:
            backbone_id, sidechain_id, compound = get_synonym[0], get_synonym[1], get_synonym[2]
            if False not in sidechain_id:
                results, counter = self.get_id_from_compound_entity(backbone_id, sidechain_id, compound, metabolite)

        return counter, results

    def get_id_from_compound_with_no_backbone(self, backbone: int, side_chains_compound: List, metabolite: Metabolite,
                                              results: Dict, counter: Dict) -> tuple():
        """Method to handle compounds that have synonyms ID's for sidechains and backbone but only sidechains ID's for compound in database. This Method will implement a search directly for the 
        generic compound with the same synonym as the backbone

        :param backbone: backbone sinonym id
        :type backbone: int
        :param side_chains_compound: sidechain synonym id
        :type side_chain_compound: list
        :param metabolite: Metabolite being analized
        :type metabolite: Metabolite
        :param results: Dictionary of the identifiers to perform the annotation to be updated in this method
        :type results: Dict
        :param counter: Dictionary of the lipid class that the algorithm has caugth to be updated in this method
        :type counter: Dict
        :return: tuple with two dictionaries,  first one for the extend of the lipid classes that the algorithm caugth and second one for the identifiers to do the annotation
        :rtype: Tuple(Dict)
        """
        backbone_compound = self.get_compound_from_synonym(backbone)
        structurally_defined_lipids = (
            self.get_compounds_with_specific_parent_set_of_components(
                backbone_compound, side_chains_compound
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
        return results, counter

    def get_id_from_compound_entity(self, backbone_ids: List, sidechain_id: List, compound: bool,
                                     metabolite: Metabolite) -> Tuple[Dict, Dict]:
        """Method to handle compounds were the backbone id is already a compound id

        :param backbone_id: backbone sinonym id
        :type backbone_id: List
        :param sidechain_id: sidechain synonym id
        :type sidechain_id: list
        :param compound: flag saying that the backbone id is or not a compound ID
        :type compound: bool
        :param metabolite: Metabolite being analized
        :type metabolite: Metabolite
        :return: tuple with two dictionaries,  first one for the extend of the lipid classes that the algorithm caugth and second one for the identifiers to do the annotation
        :rtype: Tuple(Dict)
        """
        results = {}
        counter = defaultdict(int)
        for backbone in backbone_ids:
            get_compound = self.get_compound(backbone, sidechain_id, compound)
            backbone_coumpound, side_chains_coumpound = backbone, get_compound[1]
            if get_compound is not None and backbone_coumpound != 0:
                structurally_defined_lipids = (
                    self.get_compounds_with_specific_parent_set_of_components(
                        backbone_coumpound, side_chains_coumpound
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
        return results, counter

    def get_synonym_id(self, backbone: str, side_chain: list) -> tuple:
        """
        Method that searches for the synonyms of the backbone and side chains for each lipid in the model.
        In cases where the synonym for the backbone doesn't exist it will search for the generic compound with
        the same name as the backbone:

        :param backbone: Backbone portion of the lipid name
        :type backbone: str
        :param side_chain: List of side chains from the lipid name
        :type side_chain: list
        :return: Tuple with two lists and a Boolean. First list refers to the backbone ID, second for the
            side_chains ID's. Lastly the Boolean is True when the backbone ID is from a compound.
        :rtype: Tuple
        """
        session = self.login()
        sidechain_id = []
        backbone_id = []
        backbone = backbone.replace("'", "")
        backbone = backbone.replace(" ", "")
        backbone_node = session.run(
            "match (s:Synonym) where s.synonym=$backbone return ID(s) as boimmg_id", backbone=backbone.lower())
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
            if data:
                for value in data:
                    sidechain_id.append(value.get("boimmg_id"))
            else:
                sidechain_id.append(False)

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
            if len(backbone_id) != 0:
                compound = True
                return backbone_id, sidechain_id, compound

    def get_compound(self, backbone: str, sidechains: list, flag: bool) -> tuple:
        """Method that search for the the specific compounds attached for the side chain and backbone synonym's.
        In case that the ID from the backbone is already from a compound it will search only for the side chains compounds.

        :param backbone: Backbone's ID
        :type backbone: str
        :param sidechains: List of side chains ID's
        :type sidechains: list
        :param flag: Bolean is True when the backbone ID is from a compound
        :type flag: bool
        :return: Tuple with backbone and side chains Compounds ID's
        :rtype: Tuple
        """
        backbone_compound = 0
        session = self.login()

        if not flag:
            side_chains_compound = self.get_side_chain_compounds(sidechains)
            node = session.run(
                "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$backbone and "
                "exists(t.boimmg_id) return t.boimmg_id as boimmg_id",
                backbone=backbone,
            )
            for value in node:
                backbone_compound = value.get("boimmg_id")
            if len(side_chains_compound) != 0:
                return backbone_compound, side_chains_compound

        if flag:
            side_chains_compound = self.get_side_chain_compounds(sidechains)
            if len(side_chains_compound) != 0:
                return backbone_compound, side_chains_compound

    def get_side_chain_compounds(self, side_chains) -> List[int]:
        """
        Method responsible to get database compounds from the side_chains id's

        :return: list of side chains compounds id's
        :rtype: List(int)
        """
        session = self.login()
        side_chains_compound = []
        for v in side_chains:
            result = session.run(
                "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$v and exists("
                "t.boimmg_id) return id(t)",
                v=v,
            )
            data = result.data()
            for value in data:
                if value.get("id(t)") not in side_chains_compound:
                    side_chains_compound.append(value.get("id(t)"))
        return side_chains_compound

    def get_compound_from_synonym(self, synonym_ID: int) -> int:
        """Method to get the generic compound directly connected to the backbone synonym.

        :param synonym_ID: ID of the backbone synonym
        :type synonym_ID: Generic Compound ID with specific backbone synonym
        :return: _description_
        :rtype: int
        """
        session = self.login()
        node = session.run(
            "match(s:Synonym)-[:is_synonym_of]->(t:Compound) where id(s)=$backbone and exists(t.boimmg_id) return "
            "t.boimmg_id as boimmg_id",
            backbone=synonym_ID,
        )
        for value in node:
            backbone_compound = value.get("boimmg_id")
            return backbone_compound

    def get_compounds_with_specific_parent_set_of_components(
            self, parent: int, components: list
    ) -> List[int]:
        """
        Get all not generic compounds with the specific set of components

        :param parent: Backbone compound id
        :type parent: int
        :param components: List of side chains compound ID's
        :type components: int 
        :return: List with all the componds found to the specific pack of components
        :rtype: List(str)
        """
        session = self.login()
        result = session.run("match (e:Compound)<-[:is_a]-(c:Compound)<-[:component_of]-(d:Compound) "
                            "with apoc.coll.sort(collect(id(d))) as components,e,c "
                            "where e.boimmg_id = $parent and "
                            "components = $components_par "
                            "return c.boimmg_id as id",
                             components_par=sorted(components),
                             parent=parent
                             )

        data = result.data()
        res = []
        for node in data:
            node_id = node.get("id")
            res.append(node_id)

        if len(res) != 0:
            return res
