import re
from tqdm import tqdm
from neo4j import GraphDatabase
from collections import defaultdict
from collections import Counter
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager
from boimmgpy.model_annotation._utils import set_metabolite_annotation_in_model
from joblib import Parallel, delayed
from typing import List,Tuple,Dict
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



    def login(self)->GraphDatabase.driver:
        """Method to create the connection to the database

        :return: session linkage to remote database
        :rtype: 
        """
        driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
        session = driver.session()
        return session


    def treat_data(self,dict_list:List[tuple]):
        """Method responsible for handling the dictionary fractions originating from the multiprocessing of the medel_lipids_finder method.
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
            self.converted_lipid_class_dict = {i: self.converted_lipid_class_dict.get(i, 0) + converted_lipid_class_dict.get(i, 0)
            for i in set(self.converted_lipid_class_dict).union(converted_lipid_class_dict)}
        
            self.check_annotation_dict.update(check_annotation_dict)
            self.results.update(results)


      
    def find_model_lipids(self,model:Model, n_jobs:int=-1)->Tuple(Dict):
        """ Method responsible for multiprocessing the metabolites of the model, thus accelerating the whole annotation procedure.

        :param model: Cobrapy metabolic model
        :type model: Model
        :param n_jobs: Number of CPU cores to be used for multiprocessing, defaults to -1
        :type n_jobs: int, optional
        :return:  A tuple with four dictionaries first one refers to the extent of the lipid classes present in the model, 
            the second one refers to the pre-presence of annotation on the lipids, 
            third for the extent of the classes that the algorithm was able to annotate 
            and finally a dictionary with the identifiers in the model for the lipids and their identifier in the Lipid Maps and/or Swiss Lipids database
        :rtype: Tuple(Dict)
        """
        n_iterations = len(model.metabolites)
        parallel_callback = Parallel(n_jobs)
        resultados = parallel_callback(delayed(self._find_model_lipids)(model.metabolites[i]) for i in tqdm(range(n_iterations)))
        self.treat_data(resultados)
        session = self.login()
        final_model = set_metabolite_annotation_in_model(session,self.results,model)

        return (
            self.converted_lipid_class_dict,
            self.check_annotation_dict,
            self.counter,
            self.results,
            final_model
        )


    def _find_model_lipids(self,metabolite:Metabolite)->Tuple(Dict):
        """Method that searchs for lipid metabolites in model and finds their synonyms in the BOIMMG database.
        Firstly Lipids are separeted from another model Metabolites by the regex [0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))* that searches for patterns similar to the lipids side chains representatation distinct from all other metabolites
        The second Regex  *(\([\-a-zA-Z0-9/|, ]*\)) is used to delete the side chain part of the name to get only the backbone part of the lipid

        :param metabolite: Lipid metabolite from GSM model
        :type metabolite: Metabolite
        :return: A tuple with three dictionaries, first one with sorted information about lipid class present in the model, second one with the annotation information of lipids in the model (False or True for the presence of annotation) an the last
            with the classes that the algorithm can annotate.
        :rtype: Tuple(Dict)
        """

        matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
        metabolite_name = metabolite.name
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
            check_annotation_dict = self.check_annotation(lipid=metabolite)
            backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", metabolite_name)
            if len(side_chain) != 1:
                for a in range(len(side_chain) - 1):
                    backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", backbone)
            lipid_class_dict[backbone] += 1
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


    def check_annotation(self, lipid:Metabolite)->Dict:
        """Method that checks if a given lipid is annotated in the model

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

    def search_lipid_synonyms(self, metabolite:Metabolite,side_chain:List[str])->Tuple(Dict):
        """Method that implements the screening in the database for accurate lipid structure.
        This method calls all StaticMethods to do the screening in the database.


        :param metabolite: Given Metabolite of the cobra model
        :type metabolite: Metabolite
        :param side_chain: List of the side chains referent to the given Metavolite
        :type side_chain: List[str]
        :return: tuple with two dictionaries, first one for the extend of the lipid classes that the algorithm caugth and second one for the identifiers to do the annotation
        :rtype: Tuple(Dict)
        """
        session = self.login()
        results={}
        counter = defaultdict(int)
        backbone_id, sidechain_id, compound = self.get_synonym_id(session,self.backbone, side_chain)
        if backbone_id != None and not compound:
            for backbone in backbone_id:
                get_compound = self.get_coumpound(session,backbone, sidechain_id, compound)
                backbone_coumpound, side_chains_compound = get_compound[0],get_compound[1]
                if  get_compound != None and backbone_coumpound != 0:
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(
                        session, backbone_coumpound, side_chains_compound
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
                if get_compound != None and backbone_coumpound == 0:
                    results,counter = self.get_id_from_compound_with_no_backbone(session,backbone,side_chains_compound,metabolite,results,counter)

        if backbone_id != None and compound:
            results,counter = self.get_id_from_compound_entitie(session,backbone_id,sidechain_id,compound,metabolite)
        
        
        return counter,results
    
    def get_id_from_compound_with_no_backbone(self,session:GraphDatabase.driver,backbone:str,side_chains_compound:List,metabolite:Metabolite,results:Dict,counter:Dict)->Tuple(Dict):
        """Method to handle compounds that have synonyms ID's for sidechains and backbone but only sidechains ID's for compound in database. This Method will implement a search directly for the 
        generic compound with the same synonym as the backbone

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
        :param backbone: backbone sinonym id
        :type backbone: str
        :param sidechain_id: sidechain synonym id
        :type sidechain_id: list
        :param metabolite: Metabolite being analized
        :type metabolite: Metabolite
        :return: tuple with two dictionaries,  first one for the extend of the lipid classes that the algorithm caugth and second one for the identifiers to do the annotation
        :rtype: Tuple(Dict)
        """
        backbone_compound = self.get_compound_from_synonym(session,backbone)
        structurally_defined_lipids = (
            self.get_compounds_with_specific_parent_set_of_components(session,
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
        return results,counter 
    
    def get_id_from_compound_entitie(self,session:GraphDatabase.driver,backbone_id:List,sidechain_id:List,compound:bool,metabolite:Metabolite)->Tuple(Dict):
        """Method to handle compounds were the backbone id is already a compound id

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
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
        results={}
        counter = defaultdict(int)
        for backbone in backbone_id:
            get_compound = self.get_coumpound(session,backbone, sidechain_id, compound)
            backbone_coumpound, side_chains_coumpound = get_compound[0],get_compound[1]
            if  get_compound != None and backbone_coumpound != 0:
                structurally_defined_lipids = (
                    self.get_compounds_with_specific_parent_set_of_components(
                        session, backbone_coumpound, side_chains_coumpound
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
        return results,counter 




    @staticmethod
    def get_synonym_id(session:GraphDatabase.driver,backbone: str, side_chain: list)->Tuple:
        """Method that searches for the synonyms of the backbone and side chains for each lipid in the model.
        In cases where the synonym for the backbone doesnt exist it will search for the generic compound with the same name as the backbone:

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
        :param backbone: Backbone portion of the lipid name
        :type backbone: str
        :param side_chain: List of side chains from the lipid name
        :type side_chain: list
        :return: Tuple with two lists and a Boolean. First list refers to the backbone ID, second for the side_chains ID's. Lastly the Bolean is True when the backbone ID is from a compound.
        :rtype: Tuple
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
            if len(backbone_id) != 0:
                compound = True
                return backbone_id, sidechain_id, compound


    
    def get_coumpound(self,session:GraphDatabase.driver,backbone: str, sidechains: list, flag: bool)->Tuple:
        """Method that search for the the specific compounds attached for the side chain and backbone synonym's.
        In case that the ID from the backbone is already from a compound it will search only for the side chains compounds.

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
        :param backbone: Backbone's ID
        :type backbone: str
        :param sidechains: List of side chains ID's
        :type sidechains: list
        :param flag: Bolean is True when the backbone ID is from a compound
        :type flag: bool
        :return: Tuple with backbone and side chains Compounds ID's
        :rtype: Tuple
        """
        backbone_coumpound = None

        if not flag:
            side_chains_coumpound = self.get_side_chain_compounds(session,sidechains)
            node = session.run(
                "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$backbone and exists(t.boimmg_id) return t.boimmg_id as boimmg_id",
                backbone=backbone,
            )
            for value in node:
                backbone_coumpound = value.get("boimmg_id")
            if len(side_chains_coumpound) != 0:
                return backbone_coumpound, side_chains_coumpound

        if flag:
            side_chains_coumpound = self.get_side_chain_compounds(session,sidechains)
            if len(side_chains_coumpound) != 0:
                return backbone_coumpound, side_chains_coumpound
               
    @staticmethod        
    def get_side_chain_compounds(session:GraphDatabase.driver,sidechains)->List(int):
        """Method responsible to get database componds from the side_chains id's

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
        :return: list of sidechains compounds id's 
        :rtype: List(int)
        """
        side_chains_compound = []
        for v in sidechains:
            result = session.run(
                "match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$v and exists(t.boimmg_id) return id(t)",
                v=v,
            )
            data = result.data()
            for value in data:
                if value.get("id(t)") not in side_chains_compound:
                    side_chains_compound.append(value.get("id(t)"))
        return side_chains_compound



    @staticmethod
    def get_compound_from_synonym(session:GraphDatabase.driver,synonym_ID: int)->int:
        """Method to get the generic compound directly connected to the backbone synonnym.

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
        :param synonym_ID: ID of the backbone synonym
        :type synonym_ID: Generic Compound ID with specific backbone synonym
        :return: _description_
        :rtype: int
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
        session:GraphDatabase.driver,parent: str, components: list
    )->List(str):
        """Get all non genericall compounds with the specific set of components

        :param session: driver linkage to database
        :type session: GraphDatabase.driver
        :param parent: Backbone compound id
        :type parent: str
        :param components: List of side chains compound ID's
        :type components: list
        :return: List with all the componds found to the specific pack of components
        :rtype: List(str)
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
        

            
        


