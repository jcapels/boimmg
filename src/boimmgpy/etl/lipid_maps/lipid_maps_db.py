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
    """Wrapper class to get functions that cannot be picklized

    Args:
        object (_type_): _description_
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
        """Method that iterates through all lipid maps lipids dictionary and gives it to insert_lm_database method be inserted in a neo4j database.
        This process is done using joblib multiprocessing
        """
        lipids_db = LipidMapsStructureDB()
        lipid_maps_db = lipids_db.getDatabase()
        parallel_callback = Parallel(4)
        data_treated = parallel_callback(
            delayed(self.insert_lm_database)(value)
            for key, value in tqdm(lipid_maps_db.items())
        )
        data_treated = pd.concat(data_treated)

    def insert_lm_database(self, lipid_container: List):
        """Method that receives a lipid container, acesses its information and sets a new node for each lipid with all relevant information

        Args:
            lipid_container (_type_): Lipid container
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
                        "MERGE (c:Compound:LipidMapsCompound { lipidmaps_id: $lipid_maps_id}) "
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
                            lipidmaps_id=lipid_container.getDbId(),
                            kegg_id=kegg_id,
                        )


if __name__ == "__main__":
    set = LipidMapsDB()
    set.treat_dataframe()