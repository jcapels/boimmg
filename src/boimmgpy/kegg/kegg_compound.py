from src.boimmgpy.kegg.kegg_entity import KeggEntity
import re


class KeggCompound(KeggEntity):

    def __init__(self, entry):
        super().__init__(entry)
        self.smiles = ""
        self.inchikey = ""
        self.formula = ""
        self.__get_all_information()

    def __get_all_information(self):
        info = self.get_raw_data()

        for text in info.split("\n"):
            new_line = re.sub('\s+', "\t", text)

            if re.search("\ANAME", new_line):
                feature = new_line.split("\t")
                self.name = feature[1]

            elif re.search("\AFORMULA", new_line):
                feature = new_line.split("\t")
                self.formula = feature[1]

    def get_name(self):
        return self.name

    def get_formula(self):
        return self.formula

    def get_kegg_id(self):
        return self.id

    def get_smiles(self):
        return self.smiles

    def get_inchikey(self):
        return self.inchikey
