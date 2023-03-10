from abc import ABCMeta, abstractmethod

import pandas as pd


class AbstractExtractor(metaclass=ABCMeta):
    """
    Abstract class to extract the desired data
    """

    @abstractmethod
    def extract(self) -> pd.DataFrame:
        """
        Method to extract data.
        """
        pass

    @abstractmethod
    def scrape_data(self):
        """
        Main method to scrape data.
        """
        pass


class AbstractTransformer(metaclass=ABCMeta):
    """
    Abstract class to transform the desired data
    """

    @abstractmethod
    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Method to transform data.
        """
        pass


class AbstractLoader(metaclass=ABCMeta):
    """
    Abstract class to load the desired data
    """

    @abstractmethod
    def load(self, df: pd.DataFrame):
        """
        Method to transform data.
        """
        pass


class AbstractETLPipeline(metaclass=ABCMeta):
    """
    Abstract class to represent a generic pipeline.
    """

    @abstractmethod
    def extract(self, df: pd.DataFrame):
        """
        Abstract method where the extraction method will be added.
        """
        pass

    @abstractmethod
    def transform(self, df: pd.DataFrame):
        """
        Abstract method where the transform method will be added.
        """
        pass

    @abstractmethod
    def load(self, df: pd.DataFrame):
        """
        Abstract method where the load method will be added.
        """
        pass

    @abstractmethod
    def run(self):
        """
        Abstract method to run the pipeline.
        """
        pass
