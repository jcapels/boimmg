from boimmgpy.model_annotation.model_annotator import LipidNameAnnotator
from cobra.io import read_sbml_model,write_sbml_model


def test_annotation(path):

    model = read_sbml_model(path)
    annotator = LipidNameAnnotator()
    info = annotator.find_model_lipids(model)
    print(info[3])
    model_final = info[4]   
    for metabolite in model_final.metabolites:
        print (metabolite.annotation)
    write_sbml_model(model_final, "tests/unit_tests/model_annotation_tests/test_fbc23.xml")
        


if __name__ == '__main__':
    
    test_annotation(r"tests/unit_tests/model_annotation_tests/test_fbc2.xml")