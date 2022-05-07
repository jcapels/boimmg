import unittest
from boimmgpy.etl.lipid_maps import LipidMapsExtractor

class TestETL(unittest.TestCase):

    def test_lipid_maps_extract():
        extractor = LipidMapsExtractor()
        df = extractor.extract()

if __name__ == '__main__':
    unittest.main()