from neo4j import GraphDatabase
import pandas as pd
from tqdm import tqdm
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.utilities import chemo_utilities
from boimmgpy.utilities.LipidMapsStructureDB import LipidMapsStructureDB
from boimmgpy.utilities.chemo_utilities import neutralise_charges
from boimmgpy.ontology_generators.side_chain_checker import SidechainChecker
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles, MolFromSmarts
from rdkit.Chem import AllChem,MolToInchiKey
from rdkit import Chem


class LipidMapsRelationships:
    def establish_relationships(self,core=None):
        driver = GraphDatabase.driver(uri="bolt://localhost:7687",auth=("neo4j","potassio19"))
        accessor = CompoundsDBAccessor()
        generic_targets = self.get_swiss_generics()
        lipid_maps_db = LipidMapsStructureDB()
        generic_seen = []
        for generic_smile in tqdm(generic_targets):
            if not pd.isna(generic_smile) and generic_smile not in generic_seen:
                parent_smile = generic_smile
                parent_id = accessor.get_node_id_from_smiles(parent_smile)
                filter_compounds = self.get_lm_compound_similar_to_generic(parent_smile,lipid_maps_db)
                for lm_compound in tqdm(filter_compounds):
                    lm_compound_mol = MolFromSmiles(lm_compound.getSmiles())
                    parent_neutralised_smiles,replaced = neutralise_charges(parent_smile)
                    smarts = MolFromSmarts(parent_neutralised_smiles)
                    
                    lm_compound_sidechains = AllChem.DeleteSubstructs(lm_compound_mol,smarts)
                    #check_sidechain = self.process_and_check_side_chain(sidechains=lm_compound_sidechains, compound_mol=lm_compound_mol, core_smart = core)
                    sidechain, sidechains_smiles_list = chemo_utilities.retrieve_fragments(lm_compound_sidechains)
                    compound_chains_counter = len(sidechain)
                    parent_chains_counter = parent_smile.count("*")
                    if parent_chains_counter == 0 or compound_chains_counter==parent_chains_counter:
                        compound_node = accessor.get_node_from_lipid_maps_id(lm_compound.getDbId())
                        if compound_node and 'SWISS_LIPIDS' not in compound_node.aliases.keys():
                            compound_id = compound_node.id
                            #successors = accessor.get_all_successors_by_ont_id_rel_type(lm_compound.getDbId(), "is_a")
                            #if parent_id not in successors:
                            with driver.session() as session:
                                session.run("MATCH (c:Compound),(d:Compound) "
                                                    "WHERE ID(c)=$target and ID(d) = $origin "
                                                    "MERGE (c)<-[r:is_a]-(d) "
                                                    "ON CREATE SET r.time_stamp = TIMESTAMP() "
                                                    ,
                                                    target=parent_id, origin=compound_id)
                            self.establish_components_relationships(sidechains_smiles_list, core, driver, compound_node)
                    generic_seen.append(generic_smile)
            



    def get_lm_compound_similar_to_generic(self,parent_smile,lipid_maps_db):
        lm_smiles_db = lipid_maps_db.getSmilesDatabase()
        neutralized_parent_smiles,replaced = neutralise_charges(parent_smile)
        parent_smarts = MolFromSmarts(neutralized_parent_smiles)
        smiles_keys = list(lm_smiles_db.keys())

        similar_compounds, new_smiles_keys = self.get_similar_compounds(smiles_keys,parent_smarts,lm_smiles_db)

        return similar_compounds
    

    @staticmethod
    def establish_components_relationships(sidechains_smiles_list, core, tx, compound_ont):
        for sidechain in sidechains_smiles_list:
            sidechain += core
            try:
                mol = MolFromSmiles(sidechain)

                inchikey = MolToInchiKey(mol)

                with tx.session() as session:
                    session.run("MATCH (c:Compound),(d:Compound) "
                                "WHERE ID(c)=$target and d.inchikey contains $inchikey "
                                "MERGE (c)<-[r:component_of]-(d) "
                                "ON CREATE SET r.time_stamp = TIMESTAMP() ",
                                inchikey=inchikey[:-1], target=compound_ont.id)

            except:
                print("Warning: side chains not viable")
    
    @staticmethod
    def process_and_check_side_chain(sidechains, compound_mol, core_smart,
                                 side_chain_atoms_rules={},side_chain_molecule_rules={}):
        sidechains.UpdatePropertyCache()
        Chem.GetSymmSSSR(sidechains)

        sidechain_checker = SidechainChecker()
        sidechain_checked = sidechain_checker(core_smart, compound_mol,
                                side_chain_atoms_rules,
                                side_chain_molecule_rules)

        return sidechain_checked


    @staticmethod
    def get_similar_compounds(smiles_keys,smarts,smiles_database):
        new_smiles_list = smiles_keys.copy()
        filtered_compounds = []
        for smile in smiles_keys:
            molecule = MolFromSmiles(smile)
            if molecule:
                match = molecule.GetSubstructMatch(smarts)

                if match:
                    filtered_compounds.append(smiles_database[smile])
                    new_smiles_list.remove(smile)

        smiles_keys = new_smiles_list.copy()

        return filtered_compounds, smiles_keys


    @staticmethod
    def get_swiss_generics():
        sl_nodes = pd.read_csv("entities.csv")
        list_generic_smiles = []
        for i,row in sl_nodes.iterrows():
            smiles = row["smiles"]
            generic =row["generic"]
            if generic:
                list_generic_smiles.append(smiles)
        return list_generic_smiles


if __name__ == "__main__":
    test = LipidMapsRelationships()
    test.establish_relationships(core="C(=O)O")