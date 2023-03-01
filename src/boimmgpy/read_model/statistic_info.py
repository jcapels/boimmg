from boimmgpy.read_model.case_study import LipidNameAnnotator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_info(path):
    annotator = LipidNameAnnotator(path)
    model_id = annotator.model.id
    if model_id == '':
        model_id = "iBD1106"
    info = annotator.model_lipids_finder()
    lipids_class = pd.Series(info[0])
    lipids_class = pd.DataFrame(lipids_class)
    print(lipids_class)
    original_annotations = pd.Series(info[1])
    original_annotations = pd.DataFrame(original_annotations)
    class_annotated = pd.Series(info[2])
    class_annotated = pd.DataFrame(class_annotated)
    #annotations = pd.Series(info[3])

    with pd.ExcelWriter(r"models\results\annotation_results.xlsx") as writer:
    
        # use to_excel function and specify the sheet_name and index
        # to store the dataframe in specified sheet
        print(model_id)
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
