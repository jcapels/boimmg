import unittest
import sys
sys.path.insert(1,'.')
from boimmgpy.etl.lipid_maps import LipidMapsExtractor,LipidMapsTransformer,LipidMapsLoader
from boimmgpy.etl.swiss_lipids import SwissLipidsExtractor,SwissLipidsTransformer,SwissLipidsLoader

class TestLipidMapsETL(unittest.TestCase):

    def setUp(self) -> None: 
        self.extractor = LipidMapsExtractor()
        self.transformer = LipidMapsTransformer()
        self.loader = LipidMapsLoader()
        self.scrape_data = self.extractor.extract()
        self.treated_dataframe = self.transformer.transform(self.scrape_data)

    def test_lipid_maps_extract(self):
        df=self.scrape_data
        self.assertEqual(df.shape[1],23) 

        self.assertEqual(self.treated_dataframe.shape[1],2)
    
        self.loader.load(self.treated_dataframe)


class TestSwissLipidsETL(unittest.TestCase):

    def setUp(self) -> None:  
        self.extractor = SwissLipidsExtractor()
        self.transformer = SwissLipidsTransformer()
        self.loader = SwissLipidsLoader()
        self.scrape_data = self.extractor.extract()
        self.treated_dataframe = self.transformer.transform(self.scrape_data)

    def test_swiss_lipids_extract(self):
        df = self.scrape_data
        self.assertEqual(df.shape[1],29) 
        self.assertEqual(self.treated_dataframe.shape[1],2)
        print("aqui")
        self.loader.load(self.treated_dataframe)


if __name__ == '__main__':
    unittest.main()