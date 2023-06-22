from boimmgpy.model_annotation.model_annotator import LipidNameAnnotator
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
from boimmgpy.model_annotation._utils import transform_boimmg_id_in_annotation_id


def get_info(path):
    """Function that gets all statisticall information from LipidNameAnnotator class relative to lipids class caugth and number of lipids annotated.
    This data is stored in a spreadsheet specific for each model analised.

    :param path: Path to the model to be analised
    :type path: _type_
    """
    model = read_sbml_model(path)
    annotator = LipidNameAnnotator()
    info = annotator.find_model_lipids(model)
    model_id = model.id
    if model_id == "":
        model_id = "iBD1106"
    lipids_class = pd.Series(info[0])
    lipids_class = pd.DataFrame(lipids_class)
    original_annotations = pd.Series(info[1])
    original_annotations = pd.DataFrame(original_annotations)
    class_annotated = pd.Series(info[2])
    class_annotated = pd.DataFrame(class_annotated)
    sugested_annotations = info[3]

    ########## Set annotations #############
    path = f"models/model_case_study/{model_id}.xml"
    model_final = info[4]
    write_sbml_model(model_final, path)

    ######### Set statistical information ##########
    path = f"models/results/{model_id}.xlsx"

    with pd.ExcelWriter(path) as writer:
        lipids_class.to_excel(
            writer, sheet_name=str(model_id), index=True, startrow=0, startcol=0
        )
        original_annotations.to_excel(
            writer, sheet_name=str(model_id), index=True, startrow=0, startcol=8
        )
        class_annotated.to_excel(
            writer, sheet_name=str(model_id), index=True, startrow=0, startcol=4
        )

    ######## Annotations to be curated ########
    sugested_annotations = transform_boimmg_id_in_annotation_id(sugested_annotations)
    path_txt = f"models/results/{model_id}.txt"
    output = open(path_txt, "w")
    for k, v in sugested_annotations.items():
        output.writelines(f"{k} {v}\n")


if __name__ == "__main__":
    get_info(r"models/iBD1106.xml")
    get_info(r"models/chorella/PP2016-00593DR3_Data1Model_Heterotrophy.xml")
    get_info(r"models/12918_2017_441_MOESM3_ESM.xml")
    get_info(r"models/iLB1027_lipid.xml")
