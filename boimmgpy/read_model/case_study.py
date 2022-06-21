
from pathlib import Path
from cobra.io import  read_sbml_model
import re



lplantarum_model = read_sbml_model(r"C:\Users\ampsi\OneDrive\Ambiente de Trabalho\boimmg\boimmgpy\read_model\iLB1027_lipid.xml")
for metabolite in lplantarum_model.metabolites:
    matches = re.finditer("[0-9]+:[0-9]+(\([a-zA-Z0-9,]*\))*", metabolite.name)
    metabolite_name = metabolite.name
    found = False
    for match in matches:
        #print(metabolite_name)
        found = True
        metabolite_name = metabolite_name.replace(match.string[match.start():match.end()], "")
        print(match.string[match.start():match.end()])

    if found:
        metabolite_name = re.sub(" *(\([\-a-zA-Z0-9/]*\))", "", metabolite_name)
        print(metabolite_name)