from abc import abstractmethod

from airflow.decorators import task
from airflow.models.dag import dag

from etl.abstract_classes import AbstractExtractor, AbstractTransformer, AbstractLoader, AbstractETLPipeline


class AirflowExtractor(AbstractExtractor):

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

    @abstractmethod
    def transform(self):
        """
        Transform method that implements the 'task' decorator of the airflow package
        """
        pass


class AirflowLoader(AbstractLoader):

    @abstractmethod
    def load(self):
        """
        Load method that implements the 'task' decorator of the airflow package
        """
        pass


class AirflowPipeline(AbstractETLPipeline):

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
