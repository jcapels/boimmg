import unittest 
from boimmgpy.model_annotation.model_annotator import LipidNameAnnotator
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Model, Reaction, Metabolite
from pathlib import Path


def create_test_model(metabolites_dict,path):
    model = Model('test_model')
    for key,value in metabolites_dict.items():

        metabolite = Metabolite(key, name=value,
                                        compartment='c')
        model.add_metabolites(metabolite)
    write_sbml_model(model, path)



class TestModelAnnotation(unittest.TestCase):

    def setUp(self) -> None:
        self.model_path = Path()/"tests/unit_tests/model_annotation_tests/testing_model.xml"
        self.metabolites = {'M_pail3p1601619Z_c':"1-Phosphatidyl-1D-myo-inositol 3-phosphate(16:0/16:1(9Z))","M_12dgr140160_h":"1,2-Diacyl-sn-glycerol(14:0/16:0)","M_pc205n3203n6_c":"Phosphatidylcholine(20:5(5Z,8Z,11Z,14Z,17Z)/20:3(8Z,11Z,14Z))","M_12dgr182n6183n6_c":"1,2-Diacyl-sn-glycerol(18:2(9Z,12Z)/18:3(6Z,9Z,12Z))",
                       "M_clpn182_m":"Triacylglycerol (20:4/20:3/21:0)","M_pe140180_c":"Phosphatidylethanolamine(14:0/18:0)"}
        create_test_model(self.metabolites,self.model_path)
        self.final_model_path = Path()/"tests/unit_tests/model_annotation_tests/testing_model_anotated.xml"
        self.model = read_sbml_model("tests/unit_tests/model_annotation_tests/testing_model.xml")
        
        return super().setUp()

    def test_testing_model(self):
        self.assertTrue(self.model_path.is_file())
        self.assertEqual(len(self.model.metabolites),len(self.metabolites))
        annotations = []
        for metabolite in self.model.metabolites:
            if metabolite.annotation:
                annotations.append(metabolite.annotation)
        self.assertEqual(len(annotations),0)

    
    def test_annotation_class(self):
        annotator = LipidNameAnnotator()    
        info = annotator.find_model_lipids(self.model)
        model_final = info[4]
        self.assertEqual(sum(info[0].values()),len(self.metabolites))
        self.assertEqual(len(info[1]),len(self.metabolites))
        self.assertEqual(info[0],info[2])
        self.assertSetEqual(set(map(type, info[1].values())), {bool})
        self.assertIs(type(model_final),Model)
        write_sbml_model(model_final, "tests/unit_tests/model_annotation_tests/testing_model_anotated.xml")
        self.assertTrue(self.final_model_path.is_file())
        self.assertTrue(type(info[3]) is dict)
        self.assertSetEqual(set(map(type, info[3].values())), {list})
        
if __name__== "__main__":
    unittest.main()

