from boimmgpy.etl.lipid_maps.lipid_maps_db import LipidMapsDB
from boimmgpy.etl.lipid_maps.lipid_maps_relationships import LipidMapsRelationships


def integrate_lm_db_pipeline():
    """This method implements the ETL pipeline for the integration of 
    Lipid Maps lipids in LipidGEM database"""
    set_lm_db = LipidMapsDB()
    set_lm_db.treat_dataframe()
    set_lm_rel = LipidMapsRelationships()
    set_lm_rel.establish_relationships(core="C(=O)O")