import abc


class RepresentationProblemSolver(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'swap_and_gap_fill') and
                callable(subclass.swap_and_gap_fill) and
                hasattr(subclass, 'generalize_model') and
                callable(subclass.generalize_model) and
                hasattr(subclass, 'gap_fill_model_by_target') and
                callable(subclass.gap_fill_model_by_target) and
                hasattr(subclass, 'swap_from_generic') and
                callable(subclass.swap_from_generic)
                or
                NotImplemented)

    @abc.abstractmethod
    def swap_and_gap_fill(self, target_ontology_id):
        """This method aims at swapping all the biosynthetic precursors, and the conjugated acid and base of a given target.
        It only accepts type 1 changes."""
        raise NotImplementedError

    @abc.abstractmethod
    def generalize_model(self, target_ontology_id):
        """generalize the model"""
        raise NotImplementedError

    @abc.abstractmethod
    def gap_fill_model_by_target(self, target_ontology_id: int, components_ont_ids: list):
        """

        :param target_ontology_id:
        :param components_ont_ids:
        :return:
        """
        raise NotImplementedError


    # @abc.abstractmethod
    # def swap_from_generic(self, target_generic_ontology_id: int, **kwargs):
    #     """
    #
    #     :param target_ontology_id:
    #     :param components:
    #     :return:
    #     """
    #     raise NotImplementedError
