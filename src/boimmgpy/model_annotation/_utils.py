
import re


def set_annotation(session,dictionary_results,model):
    for key,value in dictionary_results.items():
        lipid_maps_ids=[]
        swiss_lipids_ids=[]
        for id in value:
            result=session.run("match(c:Compound)where id(c)=$boimmg_id return c.lipidmaps_id,c.swiss_lipids_id as ids", boimmg_id=id)
            data= result.data()
            for node in data:
                node_lipid_maps_id = node.get("c.lipidmaps_id")
                if node_lipid_maps_id != None:
                    lipid_maps_ids.append(node_lipid_maps_id)
                
                node_swiss_lipids_id = node.get("ids")
                if node_swiss_lipids_id!= None:
                    swiss_lipids_ids.append(node_swiss_lipids_id)
        for metabolite in model.metabolites:
            if metabolite.id == key:    
                if len(lipid_maps_ids) != 0:
                    for values in lipid_maps_ids:
                        metabolite.annotation["lipidmaps"] = values
                if len(swiss_lipids_ids) != 0:
                    for values in swiss_lipids_ids:  
                        metabolite.annotation["slm"] = values 

    return model