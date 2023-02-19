from cobra.io import  read_sbml_model
import re
from tqdm import tqdm
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib as plt
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager


class LipidNameAnnotator:
    def __init__(self,model_path) -> None:
        self.model = read_sbml_model(model_path)
        self.lipid_class_dict = defaultdict(int)
        self.check_annotation_dict = {}
        #driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
        #session = driver.session()
        
    def model_lipids_finder (self):
        count = 0

        for metabolite in self.model.metabolites:
            backbone = None
            side_chain = []
            matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
            metabolite_original_name = metabolite.name
            metabolite_name=metabolite.name
            found = False
            for match in matches:
                found=True
                metabolite_name = metabolite_name.replace(match.string[match.start():match.end()], "")
                side_chain.append(match.string[match.start():match.end()])

            if found:
                self.annotation_checker(lipid=metabolite)
                backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", metabolite_name)
                if len(side_chain)!= 1:
                    for a in range(len(side_chain)-1):
                        backbone = re.sub(" *(\([\-a-zA-Z0-9/|, ]*\))", "", backbone)
                self.lipid_class_dict[backbone] += 1
                count +=1
        
        sorted_lipid_class_dict = sorted(self.lipid_class_dict.items(), key=lambda x: x[1], reverse=True)
        converted_lipid_class_dict = dict(sorted_lipid_class_dict)
        return converted_lipid_class_dict,self.check_annotation_dict
    


    def annotation_checker (self,lipid):
        annotation=lipid.annotation.keys()
        annotated=False
        self.check_annotation_dict[lipid.id]=annotated
        if "slm" in annotation or "lipidmaps" in annotation:
            annotated=True
            self.check_annotation_dict[lipid.id]=annotated        

    

annotator=LipidNameAnnotator(r"src\boimmgpy\read_model\iLB1027_lipid.xml")
dicts = annotator.model_lipids_finder()
lipids_class = dicts[0]
print(lipids_class)
original_annotations = dicts[1]
print(len(original_annotations))
print(sum(map((True).__eq__, original_annotations.values())))






def read_treat_model(self):

    counter=defaultdict(int)
    results={}
    check_annotation={}
    for metabolite in tqdm(self.model.metabolites):
        side_chain=[]
        backbone=None
        matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
        metabolite_name = metabolite.name
        found = False
        for match in matches:
            found = True
            metabolite_name = metabolite_name.replace(match.string[match.start():match.end()], "")
            side_chain.append(match.string[match.start():match.end()])
            annotation=metabolite.annotation.keys()
            annotated=False
            check_annotation[metabolite.id]=annotated
            if "slm" in annotation or "lipidmaps" in annotation:
                annotated=True
                check_annotation[metabolite.id]=annotated

        if found:
            flag=False
            backbone = re.sub(" *(\([\-a-zA-Z0-9/]*\))", "", metabolite_name)
            lipid_id=get_synonym_id(backbone,side_chain)
            if lipid_id!= None and not lipid_id[2]:
                for backbone in lipid_id[0]:
                    lipid_compound=get_coumpound(backbone,lipid_id[1],flag)
                    if lipid_compound != None:
                        structurally_defined_lipids = get_compounds_with_specific_parent_set_of_components(lipid_compound[0],lipid_compound[1])
                        if  structurally_defined_lipids:
                            if metabolite.id in results:
                                results[metabolite.id] = results[metabolite.id] + structurally_defined_lipids
                            else:
                                results[metabolite.id] = structurally_defined_lipids
                            counter[re.sub(" *(\([\-a-zA-Z0-9/]*\))", "", metabolite_name)]+=1
            
            if lipid_id!= None and lipid_id[2] == True:
                flag=True
                for backbone in lipid_id[0]:
                    lipid_compound=get_coumpound(backbone,lipid_id[1],flag)
                    if lipid_compound != None:
                        structurally_defined_lipids = get_compounds_with_specific_parent_set_of_components(lipid_compound[0],lipid_compound[1])
                        if  structurally_defined_lipids:
                            if metabolite.id in results:
                                results[metabolite.id] = results[metabolite.id] + structurally_defined_lipids
                            else:
                                results[metabolite.id] = structurally_defined_lipids
                            counter[re.sub(" *(\([\-a-zA-Z0-9/]*\))", "", metabolite_name)]+=1      
    
    print(results,len(results))        
    return results,counter,check_annotation
    

## QUERYS TO MATCH SYNONYM AND GET DATABSE ID
def get_synonym_id(backbone,side_chain):
    sidechain_id=[]
    backbone_id=[]
    backbone_node=session.run("match (s:Synonym) where s.synonym='"+str(backbone.lower())+"' return ID(s) as boimmg_id")
    node= backbone_node.data()
    compound=False
    
    for value in node:
        backbone_id.append(value.get("boimmg_id"))
            

    for synonym in side_chain:
        synonym_split = synonym.split(":")
        if synonym_split[1] == "0":
            synonym = "C" + synonym
        result =(session.run("match (s:Synonym) where s.synonym='"+str(synonym.lower())+"' return ID(s) as boimmg_id"))
        data = result.data()
        for value in data:
            sidechain_id.append(value.get("boimmg_id"))
        
    
    if len(backbone_id) !=0 and len(sidechain_id)!=0:
        return backbone_id,sidechain_id,compound
    
    if len(backbone_id) ==0 and len(sidechain_id)!=0:
        result =(session.run("match (c:Compound) where c.name='"+str(backbone)+"' return ID(c) as boimmg_id"))
        data = result.data()
        for value in data:
            backbone_id.append(value.get("boimmg_id"))
        if backbone != 0:
            compound=True
            return backbone_id,sidechain_id,compound
    
        
def get_coumpound(backbone,sidechains,flag):
    side_chains_coumpound=[]
    backbone_coumpound=0
    
    if flag==False:
        for v in sidechains:
            result=(session.run("match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$v and exists(t.boimmg_id) return id(t)",v=v))
            data = result.data()
            for value in data:
                side_chains_coumpound.append(value.get("id(t)"))
        node=session.run("match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$backbone and exists(t.boimmg_id) return t.boimmg_id as boimmg_id",backbone=backbone)
        for value in node:
                backbone_coumpound = value.get("boimmg_id")
        if backbone_coumpound !=0 and len(side_chains_coumpound)!=0:
            return backbone_coumpound, side_chains_coumpound
    
    if flag==True:
        for v in sidechains:
            result=(session.run("match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound) where id(s)=$v and exists(t.boimmg_id) return id(t)",v=v))
            data = result.data()
            for value in data:
                side_chains_coumpound.append(value.get("id(t)"))
            return backbone,side_chains_coumpound


## GET COMPOUNDS WITH SPECIFIC SETS OF COMPONENTS
def get_compounds_with_specific_parent_set_of_components( parent, components):
    result = session.run("match (e:Compound)<-[:is_a]-(c:Compound)<-[:component_of]-(d:Compound) "
                            "with collect(id(d)) as components,e,c "
                            "where e.boimmg_id = $parent and "
                            "components = $components_par "
                            "return c.boimmg_id as id ",
                            components_par=components,
                            parent=parent
                            )

    data = result.data()
    res = []
    for node in data:
        node_id = node.get("id")
        res.append(node_id)

    if len(res) != 0:
        return res