from abc import ABCMeta, abstractmethod


class AbstractExtractor(metaclass=ABCMeta):
    """
    Abstract class to extract the desired data
    """

    @abstractmethod
    def extract(self):
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
    def transform(self):
        """
        Method to transform data.
        """
        pass


class AbstractLoader(metaclass=ABCMeta):
    """
    Abstract class to load the desired data
    """

    @abstractmethod
    def load(self):
        """
        Method to transform data.
        """
        pass


class AbstractETLPipeline(metaclass=ABCMeta):
    """
    Abstract class to represent a generic pipeline.
    """

    @abstractmethod
    def extract(self):
        """
        Abstract method where the extraction method will be added.
        """
        pass

    @abstractmethod
    def transform(self):
        """
        Abstract method where the transform method will be added.
        """
        pass

    @abstractmethod
    def load(self):
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
