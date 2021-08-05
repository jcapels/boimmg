class IDConverter:

    def __init__(self):
        self.modelSeedToDb = {}
        self.dbToModelSeed = {}

    def convert_db_id_to_model_seed_by_db_id(self, db_id):
        i = 0
        keys = list(self.dbToModelSeed.keys())
        latent_id = None
        while i < len(keys):
            db = keys[i]
            ids = self.dbToModelSeed[db].keys()
            if db_id in ids:
                if db_id == "BiGG":
                    latent_id = self.dbToModelSeed[db][db_id]
                else:
                    return self.dbToModelSeed[db][db_id]

            i += 1

        if latent_id:
            return latent_id

        return None

    def convert_modelSeedId_into_other_dbID(self, modelSeedId, database_name):

        if database_name in self.modelSeedToDb.get(modelSeedId):
            return self.modelSeedToDb.get(modelSeedId).get(database_name)

        else:
            return []

    def get_all_aliases_by_modelSeedID(self, modelSeedId):
        if modelSeedId in self.modelSeedToDb.keys():
            res = self.modelSeedToDb[modelSeedId]
            if "BiGG" in res and "BiGG1" in res:
                del res["BiGG"]
                bigg_ids = res["BiGG1"]
                res["BiGG"] = bigg_ids

            res["ModelSEED"] = [modelSeedId]

            res = self.__prune_dictionary(res)

            return res
        else:
            return {}

    @staticmethod
    def __prune_dictionary(dictionary):
        res = dictionary.copy()
        for key in dictionary:
            if "" in dictionary[key]:
                del res[key]
        return res

    def convert_dbID_into_modelSeedId(self, database_name, db_id):
        if db_id in self.dbToModelSeed.get(database_name).keys():
            return self.dbToModelSeed.get(database_name).get(db_id)
        else:
            return None

    def get_modelSeedIdToDb(self):
        return self.modelSeedToDb

    def get_dbIdToModelSeed(self):
        return self.dbToModelSeed
