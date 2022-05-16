import unittest
from boimmgpy.etl.lipid_maps import LipidMapsExtractor
from boimmgpy.etl.swiss_lipids import SwissLipidsExtractor

class TestETL(unittest.TestCase):

    def test_lipid_maps_extract(self):
        extractor = LipidMapsExtractor()
        df = extractor.extract()
    
    def test_lipid_maps_extract(self):
        extractor = SwissLipidsExtractor()
        df = extractor.extract()

if __name__ == '__main__':
    unittest.main()