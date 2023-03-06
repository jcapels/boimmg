from boimmgpy.read_model.case_study import LipidNameAnnotator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model


def get_info(path):
    model = read_sbml_model(path)
    #metabolite_names = [metabolite.name for metabolite in model.metabolites]

    annotator = LipidNameAnnotator(model)
    info = annotator.model_lipids_finder()
    model_id = model.id
    if model_id == '':
        model_id = "iBD1106"
    path = (f"models\\results\{model_id}.xlsx")
    lipids_class = pd.Series(info[0])
    lipids_class = pd.DataFrame(lipids_class)
    original_annotations = pd.Series(info[1])
    original_annotations = pd.DataFrame(original_annotations)
    class_annotated = pd.Series(info[2])
    class_annotated = pd.DataFrame(class_annotated)
    #annotations = pd.Series(info[3])

    with pd.ExcelWriter(path) as writer:
        # to store the dataframe in specified sheet
        lipids_class.to_excel(writer, sheet_name=str(model_id), index=True, startrow=0 , startcol=0)
        original_annotations.to_excel(writer, sheet_name=str(model_id), index=True, startrow=0 , startcol=8)
        class_annotated.to_excel(writer, sheet_name=str(model_id), index=True,startrow=0, startcol=4)

"""
    print(lipids_class)
    print(sum(lipids_class.values()))
    print(len(original_annotations))
    print(sum(map((True).__eq__, original_annotations.values())))
"""

if __name__ == '__main__':
    get_info(r"models\iBD1106.xml")
    get_info(r"models\chorella\PP2016-00593DR3_Data1Model_Heterotrophy.xml")
    get_info(r"models\12918_2017_441_MOESM3_ESM.xml")
    get_info(r"models\iLB1027_lipid.xml")