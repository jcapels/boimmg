import re
from typing import List
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles, MolFromSmarts
from boimmgpy.etl.model_seed.compounds.ModelSeedCompound import ModelSeedCompound
from boimmgpy.etl.model_seed.compounds.model_seed_compounds_DB import ModelSeedCompoundsDB
from boimmgpy.etl.swiss_lipids.swiss_lipids_synonyms import SwissLipidsExtractor
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter




class SwissLipidsDb:


    def treat_dataframe(self):
        self.__modelSeedDB = ModelSeedCompoundsDB()
        self.__idConverter = CompoundsIDConverter()
        data = self.extract_database()
        new_data = []
        rel_data = [] 
        iteration = len(data)
        for i in tqdm(range(iteration)):
            new_line = self._treat_dataframe(data.iloc[[i]])
            if new_line is not None:
                new_data.append(new_line)
            new_relations = self.set__relationships(data.iloc[[i]]) 
            if new_relations is not None:
                rel_data.append(new_relations)       
        
       
        self.create_entities_file(new_data)
        self.create_realationships_file(rel_data)
    
    def create_realationships_file(self,data:List[pd.DataFrame]):
        data_frame = pd.concat(data)
        data_frame.to_csv("rel.csv", sep=",", index=False)
    
        
    def create_entities_file(self,data:List[list]):
        """Method to create Swiss Lipids entities files

        Args:
            data (List[list]): List of lists to create entities file
        """
        header = ["swiss_lipids_id:ID", "name", "smiles", "inchi", "inchikey", "formula", "charge", "mass",
              "hmdb_id", "chebi_id", "lipidmaps_id", "pubchem_cid", "kegg_id", "bigg_id", "metanetx_id",
              "metacyc_id", "generic", "model_seed_id"]
        data_frame = pd.DataFrame(data)
        data_frame.to_csv("entities.csv", sep=",", header=header, index=False)


    def _treat_dataframe(self,data:pd.DataFrame):
        modelSeedDB = self.__modelSeedDB
        idConverter = self.__idConverter
        flag = False
        inchikey,smiles,level = self.get_df_info(data,flag)
        flag = True
        canonical_smiles,swisslipids_id,name,formula,inchi,hmdb_id,chebi_id,lipidmaps_id,pubchem_cid,charge,mass = self.get_df_info(data,flag,smiles)
        if not pd.isna(hmdb_id):
            hmdb_id = hmdb_id.split("|")
            if hmdb_id:
                hmdb_id = hmdb_id[0]
            else:
                hmdb_id = None

        if not pd.isna(chebi_id):
            chebi_id = chebi_id.split("|")
            if chebi_id:
                chebi_id = chebi_id[0]
            else:
                chebi_id = None

        if not pd.isna(lipidmaps_id):
            lipidmaps_id = lipidmaps_id.split("|")
            if lipidmaps_id:
                lipidmaps_id = lipidmaps_id[0]
            else:
                lipidmaps_id = None

        if not pd.isna(pubchem_cid):
            pubchem_cid = pubchem_cid.split("|")
            if pubchem_cid:
                pubchem_cid = pubchem_cid[0]
            else:
                pubchem_cid = None
        
        generic = False

        model_seed_compound = None
        if pd.isna(inchikey):
            generic = True
            if pd.isna(smiles) and canonical_smiles:
                model_seed_compound = modelSeedDB.get_compound_by_canonical_smiles(canonical_smiles)

        else:
            model_seed_compound = modelSeedDB.get_compound_by_inchi_key(inchikey)

        if model_seed_compound:

            kegg_id, bigg_id, metanetx_id, metacyc_id = self.integrate_model_ids(idConverter, model_seed_compound)

            new_line = [swisslipids_id, name, canonical_smiles, inchi, inchikey, formula, charge, mass, hmdb_id,
                        chebi_id, lipidmaps_id, pubchem_cid, kegg_id, bigg_id, metanetx_id, metacyc_id, generic,
                        model_seed_compound.getDbId()]
            return new_line

        else:
            kegg_id = None
            bigg_id = None
            metanetx_id = None
            metacyc_id = None

            new_line = [swisslipids_id, name, canonical_smiles, inchi, inchikey, formula, charge, mass, hmdb_id,
                        chebi_id, lipidmaps_id, pubchem_cid, kegg_id, bigg_id, metanetx_id, metacyc_id, generic,
                        None]
            return new_line

    def set__relationships(self,data:pd.DataFrame)->List:
        rel_list = pd.DataFrame(columns=[":START_ID",":END_ID",":TYPE"])
        counter = 0
        for i,row in data.iterrows():
            compound_id = row["Lipid ID"]
            component = row["Components*"]
            if component is not None and not pd.isnull(component):
                component_splits = component.split("/")
                for split in component_splits:
                    split = re.sub(" *\((.*?)\)", "", split)
                    split = split.replace(" ", "")
                    rel_list.at[counter,":START_ID"] = split
                    rel_list.at[counter,":END_ID"] = compound_id
                    rel_list.at[counter,":TYPE"] = "component_of"
                    counter += 1
            parent = row["Lipid class*"]
            if parent is not None and not pd.isnull(parent):
                parent_splits = parent.split("|")
                for split in parent_splits:
                    split = split.replace(" ", "")
                    rel_list.at[counter,":START_ID"] = compound_id
                    rel_list.at[counter,":END_ID"] = split
                    rel_list.at[counter,":TYPE"] = "is_a"
                    counter += 1
        return rel_list


    @staticmethod
    def get_df_info(data:pd.DataFrame,flag:bool,smile=None):
        smiles = None
        if flag == False:
            for i, row in data.iterrows():
                inchikey = row["InChI key (pH7.3)"]
                if not pd.isna(inchikey):
                    inchikey = inchikey.replace("InChIKey=", "")

                else:
                    inchikey = None

                smiles = row["SMILES (pH7.3)"]
                level = row["Level"]
                if pd.isna(smiles):
                    smiles = ""
                return inchikey,smiles,level
        if flag == True:
            for i, row in data.iterrows():
                swisslipids_id = row["Lipid ID"]
                name = row["Name"]
                formula = row["Formula (pH7.3)"]
                inchi = row["InChI (pH7.3)"]
                hmdb_id = row["HMDB"]
                chebi_id = str(row["CHEBI"])
                lipidmaps_id = row["LIPID MAPS"]
                pubchem_cid = str(row["PMID"])
                charge = row["Charge (pH7.3)"]
                mass = row["Mass (pH7.3)"]
                if not pd.isna(inchi) and inchi.split("=")[1] == "none":
                    inchi = None
                canonical_smiles = None
                if not pd.isna(smile):
                    mol_smiles = MolFromSmiles(smile)
                    if mol_smiles:
                        canonical_smiles = MolToSmiles(mol_smiles)

                return canonical_smiles,swisslipids_id,name,formula,inchi,hmdb_id,chebi_id,lipidmaps_id,pubchem_cid,charge,mass



    @staticmethod
    def integrate_model_ids(idConverter, model_seed_compound:ModelSeedCompound):
        kegg_ids = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "KEGG")
        bigg_ids = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "BiGG")
        metanetxs = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "metanetx.chemical")
        metacyc = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "MetaCyc")

        if kegg_ids:
            kegg_id = kegg_ids[0]
        else:
            kegg_id = None

        if bigg_ids:
            bigg_id = bigg_ids[0]

        else:
            bigg_id = None

        if metanetxs:
            metanetx_id = metanetxs[0]

        else:
            metanetx_id = None

        if metacyc:
            metacyc_id = metacyc[0]

        else:
            metacyc_id = None

        return kegg_id, bigg_id, metanetx_id, metacyc_id
    


    @staticmethod
    def extract_database():
        extractor = SwissLipidsExtractor()
        dataframe = extractor.extract()
        return dataframe


if __name__ =='__main__':   
    db = SwissLipidsDb()
    db.treat_dataframe()