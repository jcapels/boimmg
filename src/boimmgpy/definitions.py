import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
TOOL_CONFIG_PATH = os.path.join(ROOT_DIR, 'configs/tool_configs.conf')
COMPOUNDS_ANNOTATION_CONFIGS_PATH = os.path.join(ROOT_DIR, "configs/compounds_annotation_configs.conf")
REACTIONS_ANNOTATION_CONFIGS_PATH = os.path.join(ROOT_DIR, "configs/reactions_annotation_configs.conf")
DATABASE_CONFIGS = os.path.join(ROOT_DIR, "configs/databases_files_settings.conf")
EXCEPTIONS = os.path.join(ROOT_DIR, "configs/generic_compounds_exceptions.conf")
BOIMMG_DATABASE = os.path.join(ROOT_DIR, "configs/my_database_settings.conf")
REST_ACCESS_DATABASE = os.path.join(ROOT_DIR, "configs/database_settings.conf")
