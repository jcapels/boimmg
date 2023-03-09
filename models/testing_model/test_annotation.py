from boimmgpy.model_annotation.model_annotator import LipidNameAnnotator
from cobra.io import read_sbml_model


def test_annotation(path):

    model = read_sbml_model(path)
    annotator = LipidNameAnnotator()
    info = annotator.model_lipids_finder(model)
    print (info[3])


if __name__ == '__main__':
    test_annotation(r"models/testing_model/test_fbc2.xml")