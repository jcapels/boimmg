import unittest
from boimmgpy.etl.lipid_maps import LipidMapsExtractor,LipidMapsTransformer
from boimmgpy.etl.swiss_lipids import SwissLipidsExtractor

class TestETL(unittest.TestCase):
    def test_lipid_maps_extract(self):
        extractor = LipidMapsExtractor()
        df = extractor.extract()
    
    def test_lipid_maps_extract(self):
        extractor = SwissLipidsExtractor()
        df = extractor.extract()

    def test_swiss_lipids_extract(self):
        extractor = SwissLipidsExtractor()
        df = extractor.extract()

    def test_lipid_maps_transformer(self):
        transformer=LipidMapsTransformer()
        df_treted=transformer.transform()

if __name__ == '__main__':
    unittest.main()