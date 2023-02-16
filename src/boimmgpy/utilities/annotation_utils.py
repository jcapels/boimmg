from boimmgpy.utilities import file_utilities
from boimmgpy.definitions import COMPOUNDS_ANNOTATION_CONFIGS_PATH, REACTIONS_ANNOTATION_CONFIGS_PATH


class AnnotationUtils(object):

    @staticmethod
    def get_compound_annotation_format_by_aliases(aliases):
        compoundsAnnotationConfigs = file_utilities.read_conf_file(COMPOUNDS_ANNOTATION_CONFIGS_PATH)

        res = {}
        for alias in aliases:
            if alias in compoundsAnnotationConfigs.keys():
                new_alias = compoundsAnnotationConfigs[alias]
                res[new_alias] = aliases[alias]
        return res

    @staticmethod
    def get_reaction_annotation_format_by_aliases(aliases):
        reactionsAnnotationConfigs = file_utilities.read_conf_file(
            REACTIONS_ANNOTATION_CONFIGS_PATH)

        res = {}
        for alias in aliases:
            if alias in reactionsAnnotationConfigs.keys():
                new_alias = reactionsAnnotationConfigs[alias]
                res[new_alias] = aliases[alias]
        return res

    @staticmethod
    def get_annotation_from_cobra_annotation(compound):
        compoundsAnnotationConfigs = file_utilities.read_conf_file(COMPOUNDS_ANNOTATION_CONFIGS_PATH)
        inverse = {value: key for key, value in compoundsAnnotationConfigs.items()}

        res = {}
        for annotation_key in compound.annotation:
            if annotation_key in inverse.keys():
                boimmg_database = inverse[annotation_key]

                aliases = compound.annotation[annotation_key]

                if isinstance(aliases, list):
                    res[boimmg_database] = aliases
                else:
                    res[boimmg_database] = [aliases]

        return res
