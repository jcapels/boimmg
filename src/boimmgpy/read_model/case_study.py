from cobra.io import read_sbml_model
import re
from tqdm import tqdm
from collections import defaultdict
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager

driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
session = driver.session()


class LipidNameAnnotator:
    def __init__(self, model) -> None:
        self.model = model
        self.lipid_class_dict = defaultdict(int)
        self.check_annotation_dict = {}
        self.side_chain = []
        self.backbone = None
        self.counter = defaultdict(int)
        self.results = {}

    def model_lipids_finder(self):
        """Method that searchs for lipid metabolites in model and finds their synonyms in the BOIMMG database

        Returns:
            _type_: A tuple with three dictionaries, first one with sorted information about lipid class present in the model, second one with the annotation information of lipids in the model (False or True for the presence of annotation) an the last
            with the classes that the algorithm can annotate.
        """

        count = 0
        for metabolite in tqdm(self.model.metabolites):
            matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
            metabolite_original_name = metabolite.name
            metabolite_name = metabolite.name
            backbone = None
            found = False
            self.side_chain = []
            for match in matches:
                found = True
                metabolite_name = metabolite_name.replace(
                    match.string[match.start() : match.end()], ""
                )
                self.side_chain.append(match.string[match.start() : match.end()])

            if found:
                self.annotation_checker(lipid=metabolite)
                backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", metabolite_name)
                if len(self.side_chain) != 1:
                    for a in range(len(self.side_chain) - 1):
                        backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", backbone)
                self.lipid_class_dict[backbone] += 1
                count += 1
                self.backbone = backbone
                self.search_lipid_synonyms(metabolite)

        sorted_lipid_class_dict = sorted(
            self.lipid_class_dict.items(), key=lambda x: x[1], reverse=True
        )
        converted_lipid_class_dict = dict(sorted_lipid_class_dict)
        return (
            converted_lipid_class_dict,
            self.check_annotation_dict,
            self.counter,
            self.results,
        )

    def annotation_checker(self, lipid):
        """Method that checks if a given lipid is annotated in the model

        Args:
            lipid (_type_): Given lipid metabolite from the model
        """
        annotation = lipid.annotation.keys()
        annotated = False
        self.check_annotation_dict[lipid.id] = annotated
        if "slm" in annotation or "lipidmaps" in annotation:
            annotated = True
            self.check_annotation_dict[lipid.id] = annotated

    def search_lipid_synonyms(self, metabolite):
        """Method that implements the screening in the database for accurate lipid structure.
        This method calls all StaticMethods to do the screening in the database.

        Args:
            metabolite (_type_): Lipid metabolit of the model.
        """
        is_compound = False
        lipid_id = self.get_synonym_id(self.backbone, self.side_chain)
        if lipid_id != None and not lipid_id[2]:
            for backbone in lipid_id[0]:
                lipid_compound = self.get_coumpound(backbone, lipid_id[1], is_compound)
                if lipid_compound != None and lipid_compound[0] != 0:
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(
                            lipid_compound[0], lipid_compound[1]
                        )
                    )
                    if structurally_defined_lipids:
                        if metabolite.id in self.results:
                            self.results[metabolite.id] = (
                                self.results[metabolite.id]
                                + structurally_defined_lipids
                            )
                        else:
                            self.results[metabolite.id] = structurally_defined_lipids
                        self.counter[self.backbone] += 1

                if lipid_compound != None and lipid_compound[0] == 0:
                    backbone_compound = self.get_compound_from_synonym(backbone)
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(
                            backbone_compound, lipid_compound[1]
                        )
                    )
                    if structurally_defined_lipids:
                        if metabolite.id in self.results:
                            self.results[metabolite.id] = (
                                self.results[metabolite.id]
                                + structurally_defined_lipids
                            )
                        else:
                            self.results[metabolite.id] = structurally_defined_lipids
                        self.counter[self.backbone] += 1

        if lipid_id != None and lipid_id[2] == True:
            is_compound = True
            for backbone in lipid_id[0]:
                lipid_compound = self.get_coumpound(backbone, lipid_id[1], is_compound)
                if lipid_compound != None:
                    structurally_defined_lipids = (
                        self.get_compounds_with_specific_parent_set_of_components(
                            lipid_compound[0], lipid_compound[1]
                        )
                    )
                    if structurally_defined_lipids:
                        if metabolite.id in self.results:
                            self.results[metabolite.id] = (
                                self.results[metabolite.id]
                                + structurally_defined_lipids
                            )
                        else:
                            self.results[metabolite.id] = structurally_defined_lipids
                        self.counter[self.backbone] += 1

    @staticmethod
    def get_synonym_id(backbone: str, side_chain: list):
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
    def get_coumpound(backbone: str, sidechains: list, flag: bool):
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
    def get_compound_from_synonym(synonym_ID: int):
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
        parent: str, components: list
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

