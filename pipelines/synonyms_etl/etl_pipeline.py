from boimmgpy.etl.lipid_maps import LipidMapsExtractor, LipidMapsTransformer, LipidMapsLoader
from boimmgpy.etl.swiss_lipids import SwissLipidsExtractor, SwissLipidsTransformer, SwissLipidsLoader

def run_lm_pipeline():
    extractor = LipidMapsExtractor()
    transformer = LipidMapsTransformer()
    loader = LipidMapsLoader()
    scrape_data = extractor.extract()
    treated_dataframe = transformer.transform(scrape_data)
    loader.load(treated_dataframe)

def run_sl_pipeline():
    extractor = SwissLipidsExtractor()
    transformer = SwissLipidsTransformer()
    loader = SwissLipidsLoader()
    scrape_data = extractor.extract()
    treated_dataframe = transformer.transform(scrape_data)
    loader.load(treated_dataframe)




if __name__ == '__main__':
    run_lm_pipeline()
    run_sl_pipeline()