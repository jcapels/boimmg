from cobra.io import  read_sbml_model
import re
from neo4j import GraphDatabase
from tqdm import tqdm

data_base_connection = GraphDatabase.driver(uri="bolt://localhost:7687",auth=("neo4j","potassio19"))
session = data_base_connection.session()

def read_treat_model():
    lplantarum_model = read_sbml_model(r"boimmgpy\read_model\iLB1027_lipid.xml")
    lipid_ids=[]
    lipid_coumpounds=[]
    compounds=[]
    for metabolite in lplantarum_model.metabolites:
        side_chain=[]
        backbone=None
        matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
        metabolite_name = metabolite.name
        found = False
        for match in matches:
            found = True
            metabolite_name = metabolite_name.replace(match.string[match.start():match.end()], "")
            side_chain.append(match.string[match.start():match.end()])

        if found:
            backbone = re.sub(" *(\([\-a-zA-Z0-9/]*\))", "", metabolite_name)
            lipid_id=get_synonym_id(side_chain,backbone)
            if lipid_id!= None:
                lipid_ids.append(lipid_id)
    print(lipid_ids)
    for dict in lipid_ids:
        backbone_id = dict.keys()
        side_chain_id = dict.values()
        for key in backbone_id:
            backbone_id=key
        for value in side_chain_id:
            side_chain_id = value
        lipid_coumpound=get_coumpound(backbone_id,side_chain_id)
        #compounds.append(get_compounds_with_specific_parent_set_of_components(backbone_id,side_chain_id))
        lipid_coumpounds.append(lipid_coumpound)
    print(lipid_coumpounds)

    for dict in lipid_coumpounds:
        backbone_compound_id = dict.keys()
        side_chain_compound_id = dict.values()
        for key in backbone_compound_id:
            backbone_compound_id=key
        for value in side_chain_compound_id:
            side_chain_compound_id = value
        compounds.append(get_compounds_with_specific_parent_set_of_components(backbone_compound_id,side_chain_compound_id))
    print(compounds)
    


## QUERYS TO MATCH SYNONYM AND GET DATABSE ID
def get_synonym_id(side_chain,backbone):
    lipid={}
    sidechain_id=[]
    backbone_id=0
    backbone_node=session.run("match (s:Synonym) where s.synonym='"+str(backbone)+"' return ID(s)")
    node= backbone_node.data()
    for synonym in side_chain:
        result =(session.run("match (s:Synonym) where s.synonym='"+str(synonym)+"' return ID(s)"))
        data = result.data()
        for value in data:
            sidechain_id.append(value.get("ID(s)"))
        for value in node:
            backbone_id = value.get("ID(s)")
    if backbone_id !=0 and len(sidechain_id)!=0:
        lipid[backbone_id] = sidechain_id
    if lipid:
        return lipid
      
            
        
def get_coumpound(backbone,sidechains):
    lipid_coumpound={}
    side_chains_coumpound=[]
    backbone_coumpound=0
    for v in sidechains:
        result=(session.run("match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound)where id(s)=$v return id(t)",v=v))
        data = result.data()
        for value in data:
            side_chains_coumpound.append(value.get("id(t)"))
    node=session.run("match(s:Synonym)-[:is_synonym_of]->(l)-[:is_db_link_of]->(t:Compound)where id(s)=$backbone return id(t)",backbone=backbone)
    for value in node:
            backbone_coumpound = value.get("id(t)")
    lipid_coumpound[backbone_coumpound]=side_chains_coumpound
    return lipid_coumpound


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

    return res


read_treat_model()
