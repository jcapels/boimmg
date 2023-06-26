import unittest
from boimmgpy.etl.lipid_maps.lipid_maps_synonyms import LipidMapsExtractor, LipidMapsTransformer, LipidMapsLoader
from boimmgpy.etl.swiss_lipids.swiss_lipids_synonyms import SwissLipidsExtractor, SwissLipidsTransformer, SwissLipidsLoader


class TestLipidMapsETL(unittest.TestCase):

    def setUp(self) -> None:
        """Sets up the necessary components for data extraction, transformation, and loading.

        This method initializes the `extractor`, `transformer`, and `loader` components for the LipidMaps data pipeline.
        It calls the `extract` method of the `LipidMapsExtractor` to extract data from a source.
        Then, it passes the extracted data to the `transform` method of the `LipidMapsTransformer` to perform data transformation.
        The transformed data is stored in the `treated_dataframe` attribute.

        """
        self.extractor = LipidMapsExtractor()
        self.transformer = LipidMapsTransformer()
        self.loader = LipidMapsLoader()
        self.scrape_data = self.extractor.extract()[:100]
        print(self.scrape_data)
        self.treated_dataframe = self.transformer.transform(self.scrape_data)

    def test_lipid_maps_extract(self):
        """
        Test the LipidMaps data extraction, transformation and loading.

        It checks the shape of the extracted data and the transformed data to ensure they have the expected dimensions.
        It also calls the `load` method of the `LipidMapsLoader` to test the data loading functionality.
        """
        df = self.scrape_data
        self.assertEqual(df.shape[1], 23)

        self.assertEqual(self.treated_dataframe.shape[1], 2)

        self.loader.load(self.treated_dataframe)


class TestSwissLipidsETL(unittest.TestCase):

    def setUp(self) -> None:
        """
        Sets up the necessary components for data extraction, transformation, and loading.

        This method initializes the `extractor`, `transformer`, and `loader` components for the SwissLipids data pipeline.
        It calls the `extract` method of the `SwissLipidsExtractor` to extract data from a source.
        Then, it passes the extracted data to the `transform` method of the `SwissLipidsTransformer` to perform data transformation.
        The transformed data is stored in the `treated_dataframe` attribute.
        """
        self.extractor = SwissLipidsExtractor()
        self.transformer = SwissLipidsTransformer()
        self.loader = SwissLipidsLoader()
        self.scrape_data = self.extractor.extract()[:100]
        self.treated_dataframe = self.transformer.transform(self.scrape_data)

    def test_swiss_lipids_extract(self):
        """
        Test the SwissLipids data extraction, transformation and loading.

        It checks the shape of the extracted data and the transformed data to ensure they have the expected dimensions.
        It also calls the `load` method of the `SwissLipidsLoader` to test the data loading functionality.
        """
        df = self.scrape_data
        self.assertEqual(df.shape[1], 29)
        self.assertEqual(self.treated_dataframe.shape[1], 2)
        self.loader.load(self.treated_dataframe)


if __name__ == '__main__':
    unittest.main()
