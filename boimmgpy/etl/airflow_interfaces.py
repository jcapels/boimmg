from abc import abstractmethod

import pandas as pd

from boimmgpy.etl.abstract_classes import AbstractExtractor, AbstractTransformer, AbstractLoader, AbstractETLPipeline


class AirflowExtractor(AbstractExtractor):

    @abstractmethod
    def extract(self) -> pd.DataFrame:
        """
        Extract method that implements the 'task' decorator of the airflow package
        """
        pass

    @abstractmethod
    def scrape_data(self):
        """
        Main method to scrape data.
        """
        pass


class AirflowTransformer(AbstractTransformer):

    @abstractmethod
    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Transform method that implements the 'task' decorator of the airflow package
        """
        pass


class AirflowLoader(AbstractLoader):

    @abstractmethod
    def load(self, df: pd.DataFrame):
        """
        Load method that implements the 'task' decorator of the airflow package
        """
        pass


class AirflowPipeline(AbstractETLPipeline):

    # @task
    @abstractmethod
    def extract(self, df: pd.DataFrame):
        """
        Abstract method where the extraction method will be added.
        """
        pass

    # @task
    @abstractmethod
    def transform(self, df: pd.DataFrame):
        """
        Abstract method where the transform method will be added.
        """
        pass

    # @task
    @abstractmethod
    def load(self, df: pd.DataFrame):
        """
        Abstract method where the load method will be added.
        """
        pass

    # @dag
    @abstractmethod
    def run(self):
        """
        Abstract method that MUST implement the 'dag' decorator
        :return:
        """
        pass
