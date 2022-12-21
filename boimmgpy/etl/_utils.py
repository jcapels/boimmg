import pandas as pd
from neo4j import GraphDatabase
from boimmgpy.database.accessors.database_access_manager import read_config_file


def insert_in_database_lipid_maps(df: pd.Series):
    """
    This method creates the queries necessary to upload the treated data into the database
    :param df:  Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
    :type df: pd.DataFrame
    :return: List of queries necessary to the upload of the whole dataframe
    :rtype: list
    """
    log, user, password = read_config_file()
    data_base_connection = GraphDatabase.driver(uri=log, auth=(user, password))

    with data_base_connection.session() as session:
        for i, row in df.iterrows():
            lipid_maps_id = row["LM_ID"]
            lm_synonym = row["SYNONYMS"]
            session.run('MERGE (s: Synonym {synonym:"%s"})' % str(lm_synonym))
            session.run(
                "match (l:LipidMapsCompound),(s:Synonym) where l.lipidmaps_id=$lipid_maps_id and s.synonym=$synonym "
                "merge (s)-[:is_synonym_of]->(l)",
                synonym=lm_synonym, lipid_maps_id=lipid_maps_id)


def insert_in_database_swiss_lipids(df: pd.Series):
    """
    This method creates the queries necessary to upload the treated data into the database
    :param df:  Treated pandas dataframe with a column for ID and another column for synonym and abbreviation to be load
    :type df: pd.DataFrame
    :return: List of queries necessary to the upload of the whole dataframe
    :rtype: list
    """
    log, user, password = read_config_file()
    data_base_connection = GraphDatabase.driver(uri=log, auth=(user, password))

    with data_base_connection.session() as session1:
        for i, row in df.iterrows():
            swiss_lipids_id = row["Lipid ID"]
            sl_synonym = row["Synonym"]
            session1.run('MERGE (s: Synonym {synonym:"%s"})' % str(sl_synonym))
            session1.run("match (l:SwissLipidsCompound),(s:Synonym) where l.swiss_lipids_id=$sl_id and "
                         "s.synonym=$synonym merge (s)-[:is_synonym_of]->(l)", synonym=sl_synonym,
                         sl_id=swiss_lipids_id)
