from etl.airflow_interfaces import AirflowLoader, AirflowTransformer, AirflowExtractor


class SwissLipidsExtractor(AirflowExtractor):

    def extract(self):
        pass

    def scrape_data(self):
        pass


class SwissLipidsTransformer(AirflowTransformer):

    def transform(self):
        pass


class SwissLipidsLoader(AirflowLoader):

    def load(self):
        pass