import Bio.KEGG.REST as kegg_api

class KeggEntity:

    def __init__(self, entry):
        self.id = entry
        self.retrieve_raw_data()

    def retrieve_raw_data(self):
        info = kegg_api.kegg_get(self.id)
        self.raw_data = ""

        for text in info.readlines():
           self.raw_data+=text

    def get_raw_data(self):
        return self.raw_data
