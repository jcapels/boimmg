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
        self.assertEqual(self.scrape_data.shape[1],23) 

    def test_lipid_maps_transformer(self):
        self.assertEqual(self.treated_dataframe.shape[1],2)
    
    def test_lipid_maps_loader(self):
        self.loader.load(self.treated_dataframe)


class TestSwissLipidsETL(unittest.TestCase):

    def setUp(self) -> None:  
        self.extractor = SwissLipidsExtractor()
        self.transformer = SwissLipidsTransformer()
        self.loader = SwissLipidsLoader()
        self.scrape_data = self.extractor.extract()
        self.treated_dataframe = self.transformer.transform(self.scrape_data)

    def test_swiss_lipids_extract(self):
        df = self.extractor.extract()
        self.assertEqual(df.shape[1],29) 

    def test_swiss_lipids_transformer(self):
        self.assertEqual(self.treated_dataframe.shape[1],2)
        
        
    def test_lipid_maps_loader(self):
        self.loader.load(self.treated_dataframe)



if __name__ == '__main__':
    unittest.main()