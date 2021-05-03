import re

from boimmgpy.kegg.kegg_entity import KeggEntity
import Bio.KEGG.REST as kegg_api


class KeggPathway(KeggEntity):

    def __init__(self,entry, byModules=True):
        super().__init__(entry)
        self.byModules = byModules
        self.__get_all_information()


    def __get_all_information(self):
        if self.byModules:
            info = self.get_raw_data().split("\n")
            self.modules = []
            go = False
            for text in info:

                if re.search("\AMODULE", text):
                    new_line = re.sub('\s+', "\t", text)
                    module = new_line.split("\t")[1]
                    self.modules.append(module)
                    go = True
                elif (re.search("\ADISEASE", text) or re.search("\ADBLINKS", text)
                    or re.search("\AREFERENCE", text) or re.search("\AKO_PATHWAY", text)):
                    go = False
                elif go:
                    new_line = re.sub('\s+', "\t", text)
                    module = new_line.split("\t")[1]
                    self.modules.append(module)
        else:

            data = kegg_api.kegg_link("reaction", "path:" + self.id)
            self.reactions = []
            for line in data:
                reaction = line.strip().split("\t")[1].split(":")[1]
                self.reactions.append(reaction)


    def get_modules(self):
        return self.modules

    def get_compounds(self):
        return self.compounds

    def get_reactions(self):
        return self.reactions