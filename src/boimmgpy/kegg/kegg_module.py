import re

from src.boimmgpy.kegg.kegg_entity import KeggEntity


class KeggModule(KeggEntity):

    def __init__(self, entry):
        super().__init__(entry)
        self.__get_all_information()

    def get_id(self):
        return self.id

    def __get_all_information(self):
        info = self.get_raw_data().split("\n")
        self.compounds = []
        self.reactions = []
        reaction_seen = False
        compound_seen = False

        for text in info:

            if re.search("\AREACTION", text):
                reaction_seen = True
                new_line = re.sub('\s+', "\t", text)
                reaction = new_line.split("\t")[1]
                self.reactions.append(reaction)

            elif re.search("\ACOMPOUND", text):
                if reaction_seen:
                    reaction_seen = False
                    compound_seen = True
                new_line = re.sub('\s+', "\t", text)
                compound = new_line.split("\t")[1]
                self.compounds.append(compound)

            elif reaction_seen:
                new_line = re.sub('\s+', "\t", text)
                reaction = new_line.split("\t")[1]
                self.reactions.append(reaction)

            elif compound_seen:

                if re.search("\A///", text):
                    compound_seen = False

                else:
                    new_line = re.sub('\s+', "\t", text)
                    compound = new_line.split("\t")[1]
                    if re.search("\AC", compound):
                        self.compounds.append(compound)
                    else:
                        compound_seen = False

    def get_reactions(self):
        return self.reactions

    def get_compounds(self):
        return self.compounds
