import abc

from boimmgpy.database.interfaces import node


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

    @property
    @abc.abstractmethod
    def tx(self):
        return self._tx

    @tx.setter
    @abc.abstractmethod
    def tx(self, tx):
        self._tx = tx

    @abc.abstractmethod
    def login(self):
        """Access to the database"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_predecessors_by_ont_id(self, ont_id: int) -> list:
        """Get predecessors using as parameter the database identifier"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_predecessors_by_ont_id_rel_type(self, ont_id: int, relationship_type : str) -> list:
        """Get predecessors using as parameter the database identifier and the relationship type"""

        raise NotImplementedError

    @abc.abstractmethod
    def get_node_by_ont_id(self, ont_id: int) -> node:
        """Get predecessors using as parameter the database identifier and the relationship type"""
        raise NotImplementedError