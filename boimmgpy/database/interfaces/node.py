import abc


class OntologyNode(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'extract_node_properties') and
                callable(subclass.extract_node_properties) and
                hasattr(subclass, 'id') and
                hasattr(subclass, 'name') and
                hasattr(subclass, 'inchi') and
                hasattr(subclass, 'inchikey') and
                hasattr(subclass, 'smiles') and
                hasattr(subclass, 'mass') and
                hasattr(subclass, 'generic') and
                hasattr(subclass, 'charge') and
                hasattr(subclass, 'aliases') and
                hasattr(subclass, 'formula')
                or
                NotImplemented)

    @abc.abstractmethod
    def extract_node_properties(self, node_properties: dict, other_aliases: dict):
        raise NotImplementedError
