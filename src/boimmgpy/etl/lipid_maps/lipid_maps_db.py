import pandas as pd
from tqdm import tqdm
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.utilities.LipidMapsStructureDB import LipidMapsStructureDB
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles, MolFromSmarts


class LipidMapsDB:

    def treat_dataframe():
        lipids_db = LipidMapsStructureDB()
        lipid_maps_db = lipids_db.getDatabase()
        accessor = CompoundsDBAccessor()
        for lipid in tqdm(lipid_maps_db):
            lipid_container = lipid_maps_db[lipid]
            smiles = lipid_container.getSmiles()
            if pd.isna(smiles):
                canonical_smiles = None
            else:
                try:
                    canonical_smiles = MolToSmiles(MolFromSmiles(smiles))
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

            if not res:

                with tx.session() as session:
                    session.run("MERGE (c:Compound { lipidmaps_id: $lipid_maps_id}) "
                                "ON CREATE SET c.lipid_maps_id = $lipid_maps_id, "
                                "c.smiles = $smiles, c.generic = False, c.inchikey = $inchikey, "
                                "c.formula = $formula, c.charge = 0, c.kegg_id = $kegg_id, c.inchi = $inchi, "
                                "c.chebi_id = $chebi_id, c.lipid_bank_id = $lipid_bank_id,"
                                "c.name = $name,"
                                "c.pubchem_cid = $pubchem_id, c.hmdb_id = $hmdb_id, c.time_stamp = timestamp() "
                                ,
                                lipid_maps_id=lipid,
                                smiles=canonical_smiles,
                                inchikey=lipid_container.getInchiKey(),
                                formula=lipid_container.getFormula(),
                                kegg_id=kegg_id,
                                chebi_id=chebi_id, lipid_bank_id=lipid_bank_id,
                                pubchem_id=pubchem_id, hmdb_id=hmdb_id,
                                inchi=lipid_container.getInchi(), name=name
                                )

            else:
                to_merge_node = res[0]
                if BOIMMGDatabases.LIPID_MAPS.value not in to_merge_node.aliases.keys():
                    with tx.session() as session:
                        session.run("MATCH (c:Compound) "
                                    "where id(c) = $ont_id "
                                    "set c.lipid_maps_id = $lipid_maps_id,"
                                    "c.kegg_id = $kegg_id",
                                    ont_id=to_merge_node.id,
                                    lipid_maps_id=lipid_container.getDbId(),
                                    kegg_id=kegg_id)



        pass