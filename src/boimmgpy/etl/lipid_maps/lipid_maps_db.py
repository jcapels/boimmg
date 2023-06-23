from typing import List
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from neo4j import GraphDatabase
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager
from boimmgpy.utilities.LipidMapsStructureDB import LipidMapsStructureDB
from boimmgpy.database.databases_babel import BOIMMGDatabases
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles, MolFromSmarts


class Wrapper(object):
    """Wrapper class to enable pickling of functions that cannot be pickled in multiprocessing.

    This class acts as a wrapper for calling functions that may not be picklable when using
    multiprocessing methods. It dynamically imports the module and retrieves the specified
    method to execute it.
    
    :param method_name: The name of the method to be called.
    :type method_name: str
    :param module_name: The name of the module containing the method.
    :type module_name: str
        
    """

    def __init__(self, method_name, module_name):
        
        self.method_name = method_name
        self.module_name = module_name

    def __call__(self, *args, **kwargs):
        try:
            method = __import__(
                self.module_name,
                globals(),
                locals(),
                [
                    self.method_name,
                ],
            )
            _function = getattr(method, self.method_name)
            result = _function(*args, **kwargs)
            return result
        except Exception as e:
            print("Error calling function:", self.method_name)
            print("Module:", self.module_name)
            print("Error message:", str(e))
            return None


class LipidMapsDB:
    def treat_dataframe(self):
        """Iterates through all LipidMaps lipids in a dictionary and inserts them into a Neo4j database using joblib multiprocessing.

        This method retrieves the LipidMaps lipids from a database, and then uses joblib's Parallel function to
        parallelize the insertion process into a Neo4j database. The lipids are processed in parallel, with a
        maximum of 6 concurrent processes.


        :return: None
        """
        lipids_db = LipidMapsStructureDB()
        lipid_maps_db = lipids_db.getDatabase()
        parallel_callback = Parallel(6)
        data_treated = parallel_callback(
            delayed(self.insert_lm_database)(value)
            for key, value in tqdm(lipid_maps_db.items())
        )

    def insert_lm_database(self, lipid_container: List):
        """Inserts lipid information into a Neo4j database by creating a new node for each lipid.

        The method receives a lipid container and accesses its information to create a new node in the Neo4j database.
        Relevant information such as lipidmaps_id, smiles, inchikey, formula, kegg_id, chebi_id, lipid_bank_id,
        pubchem_cid, hmdb_id, and name are used to set the attributes of the new node.

        If a compound with the same inchikey or smiles already exists in the database, the method updates the
        existing compound with the lipidmaps_id and kegg_id if they are available.

        :param lipid_container: Lipid container
        :type lipid_container: List

        """
        driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()

        accessor = CompoundsDBAccessor()
        smiles = lipid_container.getSmiles()
        if pd.isna(smiles):
            canonical_smiles = None
        else:
            try:
                mol_from_smiles = Wrapper("MolFromSmiles", "rdkit.Chem.rdmolfiles")
                mol = mol_from_smiles(smiles)
                mol_to_smiles = Wrapper("MolToSmiles", "rdkit.Chem.rdmolfiles")
                canonical_smiles = mol_to_smiles(mol)
            except:
                canonical_smiles = None
        aliases = lipid_container.getAliases()
        chebi_id = None
        lipid_bank_id = None
        pubchem_id = None
        hmdb_id = None
        kegg_id = None
        name = lipid_container.getName()

        if "KEGG" in aliases.keys():
            kegg_id = aliases.get("KEGG")[0]

        if "HMDB" in aliases.keys():
            hmdb_id = aliases.get("HMDB")[0]

        if "ChEBI" in aliases.keys():
            chebi_id = aliases.get("ChEBI")[0]

        if "LipidBank" in aliases.keys():
            lipid_bank_id = aliases.get("LipidBank")[0]

        if "PubChem" in aliases.keys():
            pubchem_id = aliases.get("PubChem")[0]

        inchikey = lipid_container.getInchiKey()

        inchikey_wlast = None
        if type(inchikey) == str:
            inchikey_wlast = inchikey[:-1]

        if inchikey_wlast:
            res = accessor.get_compound_by_inchikey(inchikey)
            node_smiles = accessor.get_compound_by_smiles(smiles)
            if not res and not node_smiles:

                with driver.session() as session:
                    session.run(
                        "MERGE (c:Compound { lipidmaps_id: $lipid_maps_id}) "
                        "ON CREATE SET c.lipidmaps_id = $lipid_maps_id, "
                        "c.smiles = $smiles, c.generic = False, c.inchikey = $inchikey, "
                        "c.formula = $formula, c.charge = 0, c.kegg_id = $kegg_id, c.inchi = $inchi, "
                        "c.chebi_id = $chebi_id, c.lipid_bank_id = $lipid_bank_id,"
                        "c.name = $name,"
                        "c.pubchem_cid = $pubchem_id, c.hmdb_id = $hmdb_id ",
                        lipid_maps_id=lipid_container.getDbId(),
                        smiles=canonical_smiles,
                        inchikey=lipid_container.getInchiKey(),
                        formula=lipid_container.getFormula(),
                        kegg_id=kegg_id,
                        chebi_id=chebi_id,
                        lipid_bank_id=lipid_bank_id,
                        pubchem_id=pubchem_id,
                        hmdb_id=hmdb_id,
                        inchi=lipid_container.getInchi(),
                        name=name,
                    )

            else:
                if BOIMMGDatabases.LIPID_MAPS.value not in res.aliases.keys():
                    with driver.session() as session:
                        session.run(
                            "MATCH (c:Compound) "
                            "where id(c) = $ont_id "
                            "set c.lipidmaps_id = $lipid_maps_id,"
                            "c.kegg_id = $kegg_id",
                            ont_id=res.id,
                            lipid_maps_id=lipid_container.getDbId(),
                            kegg_id=kegg_id,
                        )

