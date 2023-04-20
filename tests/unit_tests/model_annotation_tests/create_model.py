from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model

model = Model('test_model')


M_pail3p1601619Z_c = Metabolite('M_pail3p1601619Z_c' ,name="1-Phosphatidyl-1D-myo-inositol 3-phosphate(16:0/16:1(9Z))",compartment='c')
M_12dgr140160_h = Metabolite("M_12dgr140160_h",name="1,2-Diacyl-sn-glycerol(14:0/16:0)",compartment='c')
M_pc205n3203n6_c = Metabolite("M_pc205n3203n6_c",name="Phosphatidylcholine(20:5(5Z,8Z,11Z,14Z,17Z)/20:3(8Z,11Z,14Z))",compartment='c')
M_pe140180_c = Metabolite("M_pe140180_c",name="Phosphatidylethanolamine(14:0/18:0)",compartment='c')
M_clpn182_m = Metabolite("M_clpn182_m",name="Triacylglycerol (20:4/20:3/21:0)",compartment='c')

M_12dgr182n6183n6_c = Metabolite("M_12dgr182n6183n6_c",name="1,2-Diacyl-sn-glycerol(18:2(9Z,12Z)/18:3(6Z,9Z,12Z))",compartment='c')

metabolites = [M_12dgr140160_h,M_12dgr182n6183n6_c,M_pail3p1601619Z_c,M_pc205n3203n6_c,M_pe140180_c,M_clpn182_m]
model.add_metabolites(metabolite_list=metabolites)


write_sbml_model(model, "tests/unit_tests/model_annotation_tests/test_fbc2.xml")



