from boimmgpy.etl.lipid_maps.lipid_maps_db import LipidMapsDB
from boimmgpy.etl.lipid_maps.lipid_maps_synonyms import LipidMapsExtractor, LipidMapsTransformer, LipidMapsLoader
from boimmgpy.etl.swiss_lipids.swiss_lipids_synonyms import SwissLipidsExtractor, SwissLipidsTransformer, SwissLipidsLoader



def run_set_lipid_maps_db():
    set = LipidMapsDB()
    set.treat_dataframe()

def run_lm_pipeline():
    """ This method implements the ETL pipeline for Lipid Maps dataset 
    """
    extractor = LipidMapsExtractor()
    scrape_data = extractor.extract()
    transformer = LipidMapsTransformer()
    treated_dataframe = transformer.transform(scrape_data)
    loader = LipidMapsLoader()
    loader.load(treated_dataframe)

def run_sl_pipeline():
    """ This method implements the ETL pipeline for Swiss Lipids dataset 
    """
    extractor = SwissLipidsExtractor()
    scrape_data = extractor.extract()
    transformer = SwissLipidsTransformer()
    treated_dataframe = transformer.transform(scrape_data)
    loader = SwissLipidsLoader()
    loader.load(treated_dataframe)




if __name__ == '__main__':
    run_set_lipid_maps_db()
    run_lm_pipeline()
    run_sl_pipeline()