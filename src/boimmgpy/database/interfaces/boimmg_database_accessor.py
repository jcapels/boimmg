import abc

from src.boimmgpy.database.interfaces import node


class BOIMMGDatabaseAccessor(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'login') and
                callable(subclass.login) and
                hasattr(subclass, 'get_predecessors_by_ont_id') and
                callable(subclass.get_predecessors_by_ont_id) and
                hasattr(subclass, 'get_predecessors_by_ont_id_rel_type') and
                callable(subclass.get_predecessors_by_ont_id_rel_type) and
                hasattr(subclass, 'get_node_by_ont_id') and
                callable(subclass.get_node_by_ont_id)
                or
                NotImplemented)

    @abc.abstractmethod
    def get_predecessors_by_ont_id(self, ont_id: int) -> list:
        """Get predecessors using as parameter the database identifier"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_predecessors_by_ont_id_rel_type(self, ont_id: int, relationship_type: str) -> list:
        """Get predecessors using as parameter the database identifier and the relationship type"""

        raise NotImplementedError

    @abc.abstractmethod
    def get_node_by_ont_id(self, ont_id: int) -> node:
        """Get predecessors using as parameter the database identifier and the relationship type"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_node_from_model_seed_id(self, model_seed_id: str) -> node:
        """Get node using as parameter the model seed id"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_node_id_from_model_seed_id(self, model_seed_id: str) -> int:
        """Get database id using as parameter the model seed id"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_conjugates(self, ont_id: int) -> list:
        """Get conjugates using as parameter the database id"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_compounds_with_only_one_component(self, onto_id: int, components: list):
        """Get compounds database identifiers using as parameter a list of components"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_compounds_with_specific_parent_within_set_of_components(self, onto_id: id, components: list, sources: list):
        raise NotImplementedError

    @abc.abstractmethod
    def get_all_predecessors_by_ont_id_rel_type(self, onto_id: id, relationship_type: str):
        raise NotImplementedError

    @abc.abstractmethod
    def get_compounds_with_specific_parent_set_of_components(self, onto_id: id, components: list):
        raise NotImplementedError

    @abc.abstractmethod
    def get_leaves_from_ont_id(self, ont_id: int):
        raise NotImplementedError

    @abc.abstractmethod
    def get_successors_by_ont_id_rel_type(self, onto_id: id, relationship_type: str):
        raise NotImplementedError
