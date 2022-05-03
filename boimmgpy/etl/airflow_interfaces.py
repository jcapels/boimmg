from abc import abstractmethod
from typing import List

from airflow.decorators import task
from airflow.models.dag import dag

from etl.abstract_classes import AbstractExtractor, AbstractTransformer, AbstractLoader, AbstractETLPipeline


class AirflowExtractor(AbstractExtractor):

    @task
    @abstractmethod
    def extract(self):
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

    @task
    @abstractmethod
    def transform(self):
        """
        Transform method that implements the 'task' decorator of the airflow package
        """
        pass


class AirflowLoader(AbstractLoader):

    @task
    @abstractmethod
    def load(self):
        """
        Load method that implements the 'task' decorator of the airflow package
        """
        pass


class AirflowPipeline(AbstractETLPipeline):

    def __init__(self, steps: List):
        """
        Constructor of the ETL pipeline
        :param List steps: list of steps to run the pipeline
        """
        self._steps = steps

    @task
    @abstractmethod
    def extract(self):
        """
        Abstract method where the extraction method will be added.
        """
        pass

    @task
    @abstractmethod
    def transform(self):
        """
        Abstract method where the transform method will be added.
        """
        pass

    @task
    @abstractmethod
    def load(self):
        """
        Abstract method where the load method will be added.
        """
        pass

    @dag
    @abstractmethod
    def run(self):
        """
        Abstract method that MUST implement the 'dag' decorator
        :return:
        """
        pass
