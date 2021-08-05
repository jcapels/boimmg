import json

from boimmgpy.database.containers.compound_node import CompoundNode
from boimmgpy.definitions import COMPOUNDS_ANNOTATION_CONFIGS_PATH
from boimmgpy.utilities import file_utilities


class ModelMapper:

    def __init__(self, model, compoundsIdConverter, accessor):

        self.model = model
        self.__compoundsAnnotationConfigs = file_utilities.read_conf_file(
            COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        self.__compoundsIdConverter = compoundsIdConverter
        self.__compounds_ontology = accessor

        self.mapped = False

        self.boimmg_db_model_map = {}

        self.boimmg_db_model_map_reverse = {}

        self.compounds_aliases_indexation = {}

        self.compound_inchikey_indexation = {}

        self.__compounds_aliases_indexation_reverse = {}

    @property
    def model(self):
        return self.__model

    @model.setter
    def model(self, value):
        self.__model = value

    @property
    def boimmg_db_model_map(self):
        return self.__boimmg_db_model_map

    @boimmg_db_model_map.setter
    def boimmg_db_model_map(self, value):
        if isinstance(value, dict):
            self.__boimmg_db_model_map = value
        else:
            raise ValueError("introduce a dictionary")

    @property
    def boimmg_db_model_map_reverse(self):
        return self.__boimmg_db_model_map_reverse

    @boimmg_db_model_map_reverse.setter
    def boimmg_db_model_map_reverse(self, value):
        if isinstance(value, dict):
            self.__boimmg_db_model_map_reverse = value
        else:
            raise ValueError("introduce a dictionary")

    @property
    def compounds_aliases_indexation(self):
        return self.__compounds_aliases_indexation

    @compounds_aliases_indexation.setter
    def compounds_aliases_indexation(self, value):
        if isinstance(value, dict):
            self.__compounds_aliases_indexation = value
        else:
            raise ValueError("introduce a dictionary")

    @property
    def compounds_aliases_indexation_reverse(self):
        return self.__compounds_aliases_indexation_reverse

    @compounds_aliases_indexation_reverse.setter
    def compounds_aliases_indexation_reverse(self, value):
        if isinstance(value, dict):
            self.__compounds_aliases_indexation_reverse = value
        else:
            raise ValueError("introduce a dictionary")

    @property
    def compound_inchikey_indexation(self):
        return self.__compound_inchikey_indexation

    @compound_inchikey_indexation.setter
    def compound_inchikey_indexation(self, value):
        if isinstance(value, dict):
            self.__compound_inchikey_indexation = value
        else:
            raise ValueError("introduce a dictionary")

    @property
    def mapped(self):
        return self._mapped

    @mapped.setter
    def mapped(self, value: bool):
        self._mapped = value

    @staticmethod
    def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='|'):
        """
      Call in a loop to create terminal progress bar
      @params:
          iteration   - Required  : current iteration (Int)
          total       - Required  : total iterations (Int)
          prefix      - Optional  : prefix string (Str)
          suffix      - Optional  : suffix string (Str)
          decimals    - Optional  : positive number of decimals in percent complete (Int)
          length      - Optional  : character length of bar (Int)
          fill        - Optional  : bar fill character (Str)
          printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
      """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="", flush=True)
        # Print New Line on Complete
        if iteration == total:
            print()

    def map_model(self, database):

        if database in self.__compoundsAnnotationConfigs.keys():
            database_annotation_key = self.__compoundsAnnotationConfigs[database]
        else:
            raise ValueError("introduce a valid database name")

        j = 0
        model_metabolites = list(self.model.metabolites)
        for compound in self.model.metabolites:

            self.printProgressBar(j, len(model_metabolites))
            j += 1
            compound_annotation = compound.annotation

            annotation_keys = ["seed.compound", database_annotation_key, "boimmg.compound", "inchi_key"]

            compound_id = compound.id

            compound_id = compound_id.split("_")[0]

            boimmg_compound_found = self.map_compound(compound_id, compound)

            for key in annotation_keys:
                if key in compound_annotation.keys():

                    if isinstance(compound_annotation[key], str):
                        annotation_values = [compound_annotation[key]]
                    else:
                        annotation_values = compound_annotation[key]

                    for annotation_value in annotation_values:

                        if not boimmg_compound_found:
                            boimmg_compound_found = self.map_compound(annotation_value, compound)

                        if key not in self.compounds_aliases_indexation:
                            self.compounds_aliases_indexation[key] = {}

                        if annotation_value not in self.compounds_aliases_indexation[key]:
                            self.compounds_aliases_indexation[key][annotation_value] = []

                        if compound.id not in self.compounds_aliases_indexation_reverse:
                            self.compounds_aliases_indexation_reverse[compound.id] = {}

                        if key not in self.compounds_aliases_indexation_reverse[compound.id]:
                            self.compounds_aliases_indexation_reverse[compound.id] = {key: []}

                        if annotation_value not in self.compounds_aliases_indexation_reverse[compound.id][key]:
                            self.compounds_aliases_indexation_reverse[compound.id][key].append(annotation_value)

                        if compound.id not in self.compounds_aliases_indexation_reverse[compound.id][key]:
                            self.compounds_aliases_indexation[key][annotation_value].append(compound.id)

                    if not boimmg_compound_found:

                        if key == "boimmg.compound":
                            construction_sub_string = self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]
                            value = compound_annotation[key]

                            boimmg_id = int(value.replace(construction_sub_string, ""))

                            if boimmg_id not in self.boimmg_db_model_map_reverse.keys():
                                self.boimmg_db_model_map_reverse[boimmg_id] = []

                            self.boimmg_db_model_map_reverse[boimmg_id].append(compound.id)

                            if boimmg_id not in self.boimmg_db_model_map.keys():
                                self.boimmg_db_model_map_reverse[compound.id] = []

                            self.boimmg_db_model_map_reverse[compound.id].append(compound.id)

                        elif key == "inchi_key":
                            value = compound_annotation[key][:-1]
                            if value not in self.compound_inchikey_indexation:
                                self.compound_inchikey_indexation[value] = [compound.id]
                            else:
                                self.compound_inchikey_indexation[value].append(compound.id)

        self.mapped = True

    def map_compound(self, annotation_value, compound):

        boimmg_compound_found = False

        if "cpd" in annotation_value:
            model_seed_id = [annotation_value]
        else:
            model_seed_id = self.__compoundsIdConverter.convert_db_id_to_model_seed_by_db_id(annotation_value)

        if model_seed_id:
            found = False
            i = 0
            while not found and i < len(model_seed_id):
                node = self.__compounds_ontology.get_node_id_from_model_seed_id(model_seed_id[i])
                i += 1
                if node:
                    found = True

                    boimmg_compound_found = True
                    self.boimmg_db_model_map[compound.id] = node

                    if node in self.boimmg_db_model_map_reverse.keys():
                        self.boimmg_db_model_map_reverse[node].append(compound.id)

                    else:
                        self.boimmg_db_model_map_reverse[node] = [compound.id]

        return boimmg_compound_found

    def create_map_dump(self, destination):

        out_file = open(destination + "boimmg_db_model_map.json", "w")
        json.dump(self.boimmg_db_model_map, out_file)
        out_file.close()

        out_file = open(destination + "boimmg_db_model_map_reverse.json", "w")
        json.dump(self.boimmg_db_model_map_reverse, out_file)
        out_file.close()

        out_file = open(destination + "compounds_aliases_indexation.json", "w")
        json.dump(self.compounds_aliases_indexation, out_file)
        out_file.close()

        out_file = open(destination + "compounds_aliases_indexation_reverse.json", "w")
        json.dump(self.compounds_aliases_indexation_reverse, out_file)
        out_file.close()

        out_file = open(destination + "compound_inchikey_indexation.json", "w")
        json.dump(self.compound_inchikey_indexation, out_file)
        out_file.close()

    def upload_maps(self, folder):

        in_file = open(folder + "boimmg_db_model_map.json", "r")
        self.boimmg_db_model_map = json.load(in_file)
        in_file.close()

        in_file = open(folder + "boimmg_db_model_map_reverse.json", "r")
        temp = json.load(in_file)
        in_file.close()
        self.boimmg_db_model_map_reverse = {int(k): v for k, v in temp.items()}

        in_file = open(folder + "compounds_aliases_indexation.json", "r")
        self.compounds_aliases_indexation = json.load(in_file)
        in_file.close()

        in_file = open(folder + "compounds_aliases_indexation_reverse.json", "r")
        self.compounds_aliases_indexation_reverse = json.load(in_file)
        in_file.close()

        in_file = open(folder + "compound_inchikey_indexation.json", "r")
        self.compound_inchikey_indexation = json.load(in_file)
        in_file.close()

        self.mapped = True

    def add_new_reactions_to_model(self, new_reactions):

        for reaction in new_reactions:
            self.add_new_metabolites_to_maps(reaction.metabolites)

    def __add_metabolites_to_boimmg_indexation(self, key, annotation_value, metabolite):

        if key == self.__compoundsAnnotationConfigs.get("ModelSEED"):
            model_seed_id = [annotation_value]

        else:
            model_seed_id = self.__compoundsIdConverter.convert_db_id_to_model_seed_by_db_id(
                annotation_value)

        if model_seed_id:

            found = False
            i = 0
            while not found and i < len(model_seed_id):
                node = self.__compounds_ontology.get_node_id_from_model_seed_id(model_seed_id[i])
                i += 1
                if node:
                    # found = True

                    self.boimmg_db_model_map[metabolite.id] = node

                    if node in self.boimmg_db_model_map_reverse.keys():
                        self.boimmg_db_model_map_reverse[node].append(metabolite.id)

                    else:
                        self.boimmg_db_model_map_reverse[node] = [metabolite.id]

                    return True

        return False

    def __add_metabolites_to_aliases_indexation(self, annotation_pair, annotation_value, metabolite):

        if annotation_pair not in self.compounds_aliases_indexation:
            self.compounds_aliases_indexation[annotation_pair] = {}

        if annotation_value not in self.compounds_aliases_indexation[annotation_pair]:
            self.compounds_aliases_indexation[annotation_pair][annotation_value] = []

        if metabolite.id not in self.compounds_aliases_indexation_reverse:
            self.compounds_aliases_indexation_reverse[metabolite.id] = {}

        if annotation_pair not in self.compounds_aliases_indexation_reverse[metabolite.id]:
            self.compounds_aliases_indexation_reverse[metabolite.id] = {annotation_pair: []}

        if annotation_value not in self.compounds_aliases_indexation_reverse[metabolite.id][annotation_pair]:
            self.compounds_aliases_indexation_reverse[metabolite.id][annotation_pair].append(annotation_value)

        if metabolite not in self.compounds_aliases_indexation[annotation_pair][annotation_value]:
            self.compounds_aliases_indexation[annotation_pair][annotation_value].append(metabolite.id)

    def add_new_metabolites_to_maps(self, new_metabolites):

        for metabolite in new_metabolites:

            compound_annotation = metabolite.annotation

            boimmg_compound_found = False
            for key in compound_annotation.keys():
                if key in self.__compoundsAnnotationConfigs.keys():
                    annotation_pair = self.__compoundsAnnotationConfigs[key]

                    if isinstance(compound_annotation[key], str):
                        annotation_values = [compound_annotation[key]]
                    else:
                        annotation_values = compound_annotation[key]

                    for annotation_value in annotation_values:

                        if not boimmg_compound_found:
                            boimmg_compound_found = self.__add_metabolites_to_boimmg_indexation(key,
                                                                                                annotation_value,
                                                                                                metabolite)

                        self.__add_metabolites_to_aliases_indexation(annotation_pair, annotation_value, metabolite)

                    if not boimmg_compound_found:

                        if key == "boimmg.compound":
                            construction_sub_string = self.__compoundsAnnotationConfigs["BOIMMG_ID_CONSTRUCTION"]
                            value = compound_annotation[key]

                            boimmg_id = int(value.replace(construction_sub_string, ""))

                            self.boimmg_db_model_map[metabolite.id] = boimmg_id

                            if boimmg_id not in self.boimmg_db_model_map_reverse.keys():
                                self.boimmg_db_model_map_reverse[boimmg_id] = [metabolite.id]

                            else:
                                self.boimmg_db_model_map_reverse[boimmg_id].append(metabolite.id)

                        elif key == "inchi_key":
                            value = compound_annotation[key][:-1]
                            if value not in self.compound_inchikey_indexation:
                                self.compound_inchikey_indexation[value] = [metabolite.id]
                            else:
                                self.compound_inchikey_indexation[value].append(metabolite.id)

    def update_maps(self, old_inchikey, new_inchikey, old_id, new_id, new_ontology_id, old_aliases, new_aliases):

        compound_container = self.model.metabolites.get_by_id(new_id)

        if old_id in self.boimmg_db_model_map:
            del self.boimmg_db_model_map[old_id]
            self.boimmg_db_model_map[new_id] = new_ontology_id

        for key, value in self.boimmg_db_model_map_reverse.items():
            if old_id in value:
                del self.boimmg_db_model_map_reverse[key]
                break

        if new_ontology_id not in self.boimmg_db_model_map_reverse:
            self.boimmg_db_model_map_reverse[new_ontology_id] = [compound_container.id]

        elif compound_container.id not in self.boimmg_db_model_map_reverse[new_ontology_id]:
            self.boimmg_db_model_map_reverse[new_ontology_id].append(compound_container.id)

        if old_inchikey:
            if old_inchikey[:-1] in self.compound_inchikey_indexation:
                del self.compound_inchikey_indexation[old_inchikey[:-1]]

        if new_inchikey:
            if new_inchikey[:-1] in self.compound_inchikey_indexation:
                self.compound_inchikey_indexation[new_inchikey[:-1]].append(compound_container.id)
            else:
                self.compound_inchikey_indexation[new_inchikey[:-1]] = [compound_container.id]

        self.__update_aliases_indexation(old_id, new_aliases, compound_container)

    def __update_aliases_indexation(self, old_id, new_aliases, compound_container):

        if old_id in self.compounds_aliases_indexation_reverse:
            old_aliases = self.compounds_aliases_indexation_reverse[old_id]
            for db in self.compounds_aliases_indexation:
                if db in old_aliases:
                    for alias in old_aliases[db]:
                        if alias in self.compounds_aliases_indexation[db]:
                            del self.compounds_aliases_indexation[db][alias]

            del self.compounds_aliases_indexation_reverse[old_id]

        for key in new_aliases:
            if key in self.__compoundsAnnotationConfigs.keys() and \
                    key in self.compounds_aliases_indexation:

                for alias in new_aliases[key]:
                    if alias in self.compounds_aliases_indexation[key]:
                        self.compounds_aliases_indexation[key][alias].append(compound_container.id)
                    else:
                        self.compounds_aliases_indexation[key][alias] = [compound_container.id]

    def check_metabolites_in_model(self, inchikey: str, aliases: dict, boimmg_container: CompoundNode = None,
                                   boimmg_id: int = None) -> list:
        """
        This method checks whether a given metabolite with the :param inchikey and the :param aliases.
        If a BOIMMG node is available the search starts with it.

        :param str inchikey: InchiKey of the metabolite to be searched.
        :param dict aliases: databases links of the metabolite to be searched.
        :param CompoundNode boimmg_container: BOIMMG node
        :param int boimmg_id: boimmg identifier
        :return list: metabolites in model
        """

        if not boimmg_container and not boimmg_id:

            if inchikey:
                inchikey_without_protonation = inchikey[:-1]
                if inchikey_without_protonation in self.compound_inchikey_indexation.keys():
                    compound_ids = self.compound_inchikey_indexation[inchikey_without_protonation]
                    compounds_in_model = [self.model.metabolites.get_by_id(compound_id)
                                          for compound_id in compound_ids]

                    return compounds_in_model

            for key in aliases:
                new_key = None
                if key in self.__compoundsAnnotationConfigs:
                    new_key = self.__compoundsAnnotationConfigs[key]

                if key in self.compounds_aliases_indexation:
                    for alias in aliases[key]:
                        if alias in self.compounds_aliases_indexation.get(key):
                            compound_ids = self.compounds_aliases_indexation[key].get(alias)
                            compounds_in_model = [self.model.metabolites.get_by_id(compound_id)
                                                  for compound_id in compound_ids]

                            return compounds_in_model

                elif new_key and new_key in self.compounds_aliases_indexation:
                    for alias in aliases[key]:
                        if alias in self.compounds_aliases_indexation.get(new_key):
                            compound_ids = self.compounds_aliases_indexation[new_key].get(alias)
                            compounds_in_model = [self.model.metabolites.get_by_id(compound_id)
                                                  for compound_id in compound_ids]

                            return compounds_in_model

        elif boimmg_container:
            ont_id = boimmg_container.id
            res = []

            if ont_id in self.__boimmg_db_model_map_reverse.keys():
                for model_id in self.__boimmg_db_model_map_reverse[ont_id]:
                    res.append(self.model.metabolites.get_by_id(model_id))

            return res

        else:
            res = []

            if boimmg_id in self.__boimmg_db_model_map_reverse.keys():
                for model_id in self.__boimmg_db_model_map_reverse[boimmg_id]:
                    res.append(self.model.metabolites.get_by_id(model_id))

            return res

        return []

    def check_if_boimmg_metabolite_in_model(self, boimmg_id, aliases=None):

        if aliases is None:
            aliases = {}
        if boimmg_id in self.boimmg_db_model_map_reverse.keys():
            return self.boimmg_db_model_map_reverse[boimmg_id]

        elif aliases:
            for db in aliases:
                db_aliases = aliases[db]

                for alias in db_aliases:
                    if db in self.compounds_aliases_indexation.keys():
                        if alias in self.compounds_aliases_indexation[db]:
                            compounds = self.compounds_aliases_indexation[db][alias]
                            return compounds
        else:
            return None

        return None

    def get_boimmg_id_from_model_compound_id(self, model_id):

        if model_id in self.boimmg_db_model_map.keys():
            return self.boimmg_db_model_map[model_id]

        return None
