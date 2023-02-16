from boimmgpy.etl.lipid_maps import LipidMapsExtractor, LipidMapsTransformer, LipidMapsLoader
from boimmgpy.etl.swiss_lipids import SwissLipidsExtractor, SwissLipidsTransformer, SwissLipidsLoader



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
    run_lm_pipeline()
    run_sl_pipeline()