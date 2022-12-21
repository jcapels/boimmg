from enum import Enum


class BOIMMGDatabases(Enum):
    HMDB = "hmdb_id"
    LIPID_MAPS = "lipid_maps_id"
    MODEL_SEED = "model_seed_id"
    METACYC = "metacyc_id"
    SWISS_LIPIDS = "swiss_lipids_id"
    KEGG_ID = "kegg_id"
    CHEBI = "chebi_id"
    METANETX_ID = "metanetx_id"
    BIGG = "bigg_id"
    PUBMED = "pubmed_cid"


class ModelSEEDDatabases(Enum):
    MODEL_SEED = "ModelSEED"
    BIGG = "BiGG"
    KEGG_ID = "KEGG"
    METACYC = "MetaCyc"
    METANETX_ID = "metanetx.chemical"
    SWISS_LIPIDS = "SWISS_LIPIDS"
    LIPID_MAPS = "LIPID_MAPS"


class AliasesTransformer:

    @staticmethod
    def transform_dictionary_in_boimmg_databases(aliases):
        boimmg_format_aliases = {}
        for database in ModelSEEDDatabases:
            if database.value in aliases:
                alias_lst = aliases.get(database.value)
                boimmg_alias = BOIMMGDatabases[database.name].value

                boimmg_format_aliases[boimmg_alias] = alias_lst[0]

        return boimmg_format_aliases

    @staticmethod
    def transform_boimmg_aliases_into_model_seed(aliases):
        boimmg_format_aliases = {}
        model_seed_dbs = [database.name for database in ModelSEEDDatabases]

        for database in BOIMMGDatabases:
            if database.value in aliases:
                alias_lst = aliases.get(database.value)
                if database.name in model_seed_dbs:
                    boimmg_alias = ModelSEEDDatabases[database.name].value

                    boimmg_format_aliases[boimmg_alias] = alias_lst

        return boimmg_format_aliases

    @staticmethod
    def convert_model_seed_format_into_boimmg(database_name: str) -> str:
        for database in ModelSEEDDatabases:
            if database.value == database_name:
                return BOIMMGDatabases[database.name].value
