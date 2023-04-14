from tqdm import tqdm
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager
from boimmgpy.etl.model_seed.compounds.model_seed_compounds_DB import (
    ModelSeedCompoundsDB,
)


def integrate_all_model_seed_compounds():
    driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
    ms_compounds_db = ModelSeedCompoundsDB()

    accessor = CompoundsDBAccessor()

    for smiles in tqdm(ms_compounds_db.isomer_smiles_database):
        if "*" in smiles:
            if len(ms_compounds_db.isomer_smiles_database[smiles]) > 1:

                isomer_list = [
                    iso.getDbId()
                    for iso in ms_compounds_db.isomer_smiles_database[smiles]
                ]

                chosen_boimmg_id = None
                for iso in isomer_list:
                    node = accessor.get_node_from_model_seed_id(iso)

                    if node:
                        chosen_boimmg_id = node.id
                        break

                if chosen_boimmg_id:
                    for iso in isomer_list:
                        with driver.session() as session:
                            session.run(
                                "match (c:Compound) where id(c) = $node_id "
                                "merge (c)<-[:is_db_link_of]-(:ModelSeedCompound {model_seed_id: $iso})",
                                node_id=chosen_boimmg_id,
                                iso=iso,
                            )


integrate_all_model_seed_compounds()
