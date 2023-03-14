
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles, MolFromSmarts
from boimmgpy.etl.model_seed.compounds.ModelSeedCompoundsDB import ModelSeedCompoundsDB
from boimmgpy.etl.swiss_lipids.swiss_lipids_synonyms import SwissLipidsExtractor
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter




class SwissLipidsDb:
    
    def treat_dataframe(self):
        data = self.extract_database()   
        iteration = len(data)
        data_treated = []
        for i in tqdm(range(iteration)):
            data_ = self._treat_dataframe(data.iloc[[i]])
            data_treated.append(data_)
        data_treated = pd.concat(data_treated)


    def _treat_dataframe(self,data):
        modelSeedDB = ModelSeedCompoundsDB()
        idConverter = CompoundsIDConverter()
        flag = False
        inchikey,smiles,level = self.get_df_info(data,flag)
        if "*" not in smiles or "Class" in level:
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

                return [swisslipids_id, name, canonical_smiles, inchi, inchikey, formula, charge, mass, hmdb_id,
                            chebi_id, lipidmaps_id, pubchem_cid, kegg_id, bigg_id, metanetx_id, metacyc_id, generic,
                            model_seed_compound.getDbId()]

            else:
                kegg_id = None
                bigg_id = None
                metanetx_id = None
                metacyc_id = None

                return [swisslipids_id, name, canonical_smiles, inchi, inchikey, formula, charge, mass, hmdb_id,
                            chebi_id, lipidmaps_id, pubchem_cid, kegg_id, bigg_id, metanetx_id, metacyc_id, generic,
                            None]



    @staticmethod
    def get_df_info(data,flag,smile=None):
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
    def integrate_model_ids(idConverter, model_seed_compound):
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



db = SwissLipidsDb()
db.treat_dataframe()