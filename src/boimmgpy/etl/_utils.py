import pandas as pd

from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager


def insert_in_database_lipid_maps(df: pd.Series):
    """
    This method creates the queries necessary to upload the treated data into the database
    :param df:  Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
    :type df: pd.DataFrame
    :return: List of queries necessary to the upload of the whole dataframe
    :rtype: list
    """

    driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()

    with driver.session() as session:
        for i, row in df.iterrows():
            lipid_maps_id = row["LM_ID"]
            lm_synonym = row["SYNONYMS"]
            session.run('MERGE (s: Synonym {synonym:"%s"})' % str(lm_synonym))
            session.run(
                "match (l:LipidMapsCompound),(s:Synonym) where l.lipidmaps_id=$lipid_maps_id and s.synonym=$synonym "
                "merge (s)-[:is_synonym_of]->(l)",
                synonym=lm_synonym,
                lipid_maps_id=lipid_maps_id,
            )


def insert_in_database_swiss_lipids(df: pd.Series):
    """
    This method creates the queries necessary to upload the treated data into the database
    :param df:  Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
    :type df: pd.DataFrame
    :return: List of queries necessary to the upload of the whole dataframe
    :rtype: list
    """

    driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()

    with driver.session() as session1:
        for i, row in df.iterrows():
            swiss_lipids_id = row["Lipid ID"]
            sl_synonym = row["Synonym"]
            session1.run('MERGE (s: Synonym {synonym:"%s"})' % str(sl_synonym))
            session1.run(
                "match (l:SwissLipidsCompound),(s:Synonym) where l.swiss_lipids_id=$sl_id and "
                "s.synonym=$synonym merge (s)-[:is_synonym_of]->(l)",
                synonym=sl_synonym,
                sl_id=swiss_lipids_id,
            )


"""CALL apoc.periodic.iterate(
  "match (c:Compound) where c.swiss_lipids_id is not null return c",
  "merge (c)<-[l:is_db_link_of]-(lm:SwissLipidsCompound {swiss_lipids_id: c.swiss_lipids_id})",
  {batchSize:10000, parallel:true})
"""
