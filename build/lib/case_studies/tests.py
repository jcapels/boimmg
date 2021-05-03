import time

import cobra
from cobra import Reaction

from boimmgpy import definitions
from boimmgpy.boimmg.RepresentationRevisors.RepresentationRedundantCaseSolver import RepresentationRedundantCaseSolver
from boimmgpy.boimmg.RepresentationRevisors.RepresentationSimpleCaseSolver import SimpleCaseSolver

# ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
from boimmgpy.case_studies.redundant_representation_case.map_ecoli_model import map_iJR904


def simulate_ecoli_generalize():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iML1515.xml")
    rep_case_solver = SimpleCaseSolver(model, "BiGG")
    rep_case_solver.swap_compound("cpd15560", "cpd11669")
    rep_case_solver.swap_compound("cpd15500", "cpd11451")
    cobra.io.write_sbml_model(rep_case_solver.model, "simple_representation_case/iML1515_generalized.xml")

def simulate_ecoli_granulate():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/case_studies/iML1515_generalized.xml")
    gapfiller = SimpleCaseSolver(model, "BiGG")
    gapfiller.swap_compound("cpd11669","cpd15290")
    gapfiller.swap_compound("cpd11451", "cpd15500")
    cobra.io.write_sbml_model(gapfiller.model, "simple_representation_case/iML1515_granulated.xml")

def simulate_ecoli_granulate_boimmg():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/case_studies/iML1515_generalized.xml")
    gapfiller = SimpleCaseSolver(model, "BiGG")
    gapfiller.swap_compound("cpd11669", "C_BOIMMG_749196")
    # gapfiller.swap_compound("cpd11451", "cpd15500")
    cobra.io.write_sbml_model(gapfiller.model, "simple_representation_case/iML1515_granulated_ubi5.xml")

def simulate_sacc_generalize():
    model = cobra.io.read_sbml_model("C:/Users/Biosystems/Desktop/thesis_jcapela/BOIMMGpy/models/iMM904.xml")
    gapfiller = SimpleCaseSolver(model, "BiGG")
    gapfiller.swap_compound("cpd15290", "cpd11669")
    cobra.io.write_sbml_model(gapfiller.model, "simple_representation_case/iMM904_generalized.xml")

def simulate_sacc_granulate():
    model = cobra.io.read_sbml_model("/BOIMMGpy/case_studies/simple_representation_case/iMM904_generalized.xml")
    gapfiller = SimpleCaseSolver(model, "BiGG")
    gapfiller.swap_compound("cpd11669" , "cpd15290")
    cobra.io.write_sbml_model(gapfiller.model, "simple_representation_case/iMM904_granulated.xml")

def convert_into_json():

    model = cobra.io.read_sbml_model(
        definitions.ROOT_DIR + "/case_studies/iMM904_generalized.xml")

    cobra.io.save_json_model(model,"iMM904_generalized.json")
    # model2 = cobra.io.read_sbml_model(
    #     ROOT_DIR + "/BOIMMGpy/case_studies/iML1515_generalized.xml")
    #
    # model3 = cobra.io.read_sbml_model(
    #     ROOT_DIR + "/BOIMMGpy/case_studies/iML1515_generalized.xml")

def redundant_granulator_ecoli():
    model = map_iJR904()

    components = ["cpd00214", "cpd01080", "cpd03847", "cpd05274"]

    solver = RepresentationRedundantCaseSolver(model, "BiGG")

    solver.swap_from_generic(["cpd22513", "cpd15649"], components, True,True)


    cobra.io.write_sbml_model(solver.model, "enhanced_model_ecoli_latest_version.xml")
    cobra.io.save_json_model(solver.model, "enhanced_model_ecoli.json")

def redundant_granulator_ecoli_without_components():
    model = map_iJR904()

    components = ["cpd00214", "cpd03847", "cpd05274","cpd25615","cpd05237"]

    solver = RepresentationRedundantCaseSolver(model, "BiGG")

    solver.swap_from_generic(["cpd22513", "cpd15649"], components, True)
    solver.generateISAreactions()

    cobra.io.write_sbml_model(solver.model, "redundant_representation_case/enhanced_model_ecoli_without_components.xml")
    cobra.io.save_json_model(solver.model, "enhanced_model_ecoli.json")

def redundant_granulator_ecoli_without_components_isa_reactions_for_all():
    model = map_iJR904()

    components = ["cpd00214", "cpd03847", "cpd05274","cpd25615","cpd05237"]

    solver = RepresentationRedundantCaseSolver(model, "BiGG")

    solver.swap_from_generic(["cpd22513", "cpd15649"], components, True)
    solver.generateTestISAreactions()

    cobra.io.write_sbml_model(solver.model, "enhanced_model_ecoli_without_components_isa_reactions_for_all.xml")
    cobra.io.save_json_model(solver.model, "enhanced_model_ecoli.json")

def gap_fill_ecoli_model():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/case_studies/enhanced_model_ecoli_without_components.xml")

    components = ["palmACP_c", "myrsACP_c", "hdeACP_c", "octeACP_c","tdeACP_c"]

    PASYN_EC(model)
    LPLIPA2(model)
    LPLIPA1(model)
    LPLIPA3(model)
    PLIPA1(model)


    # model.solver = "cplex"
    # gapfiller = GapFiller(model,exchange_reactions=True, demand_reactions=True,integer_threshold=1e-300)
    # reactions = gapfiller.fill(5)
    # print(reactions)

    # model.add_reactions(reactions[0])
    cobra.io.write_sbml_model(model,
        definitions.ROOT_DIR + "/case_studies/enhanced_model_ecoli_without_components_gap_filled.xml")

def PLIPA1(model):
    s1 = {model.metabolites.get_by_id("C_BOIMMG_427_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1,
          model.metabolites.get_by_id("C_BOIMMG_314_c"): 1}

    s2 = {model.metabolites.get_by_id("C_BOIMMG_7454_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1,
          model.metabolites.get_by_id("C_BOIMMG_16363_c"): 1}

    s3 = {model.metabolites.get_by_id("C_BOIMMG_7408_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1,
          model.metabolites.get_by_id("C_BOIMMG_16317_c"): 1}

    s4 = {model.metabolites.get_by_id("C_BOIMMG_7474_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1,
          model.metabolites.get_by_id("C_BOIMMG_452_c"): 1}

    s5 = {model.metabolites.get_by_id("C_BOIMMG_7456_c"): -1,
          model.metabolites.get_by_id("ttdcea_c"): 1,
          model.metabolites.get_by_id("C_BOIMMG_16364_c"): 1}

    r1 = Reaction("PLIPA1_1")
    r2 = Reaction("PLIPA1_2")
    r3 = Reaction("PLIPA1_3")
    r4 = Reaction("PLIPA1_4")
    r5 = Reaction("PLIPA1_5")

    s = [s1, s2, s3, s4, s5]
    r = [r1, r2, r3, r4, r5]

    for i in range(len(r)):
        si = s[i]

        ri = r[i]

        si[model.metabolites.get_by_id("h2o_c")] = -1
        si[model.metabolites.get_by_id("h_c")] = 1

        ri.add_metabolites(si)

    model.add_reactions(r)

def LPLIPA1(model):
    s1 = {model.metabolites.get_by_id("C_BOIMMG_314_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1}

    s2 = {model.metabolites.get_by_id("C_BOIMMG_16363_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1}

    s3 = {model.metabolites.get_by_id("C_BOIMMG_16317_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1}

    s4 = {model.metabolites.get_by_id("C_BOIMMG_452_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1}

    s5 = {model.metabolites.get_by_id("C_BOIMMG_16364_c"): -1,
          model.metabolites.get_by_id("ttdcea_c"): 1}

    r1 = Reaction("LPLIPA1_1")
    r2 = Reaction("LPLIPA1_2")
    r3 = Reaction("LPLIPA1_3")
    r4 = Reaction("LPLIPA1_4")
    r5 = Reaction("LPLIPA1_5")

    s = [s1, s2, s3, s4, s5]
    r = [r1, r2, r3, r4, r5]

    for i in range(len(r)):
        si = s[i]

        ri = r[i]

        si[model.metabolites.get_by_id("h2o_c")] = -1
        si[model.metabolites.get_by_id("h_c")] = 1
        si[model.metabolites.get_by_id("g3pg_c")] = 1

        ri.add_metabolites(si)

    model.add_reactions(r)

def LPLIPA3(model):
    s1 = {model.metabolites.get_by_id("C_BOIMMG_223_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1}

    s2 = {model.metabolites.get_by_id("C_BOIMMG_6942_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1}

    s3 = {model.metabolites.get_by_id("C_BOIMMG_6896_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1}

    s4 = {model.metabolites.get_by_id("C_BOIMMG_319_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1}

    s5 = {model.metabolites.get_by_id("C_BOIMMG_6943_c"): -1,
          model.metabolites.get_by_id("ttdcea_c"): 1}

    r1 = Reaction("LPLIPA3_1")
    r2 = Reaction("LPLIPA3_2")
    r3 = Reaction("LPLIPA3_3")
    r4 = Reaction("LPLIPA3_4")
    r5 = Reaction("LPLIPA3_5")

    s = [s1, s2, s3, s4, s5]
    r = [r1, r2, r3, r4, r5]

    for i in range(len(r)):
        si = s[i]

        ri = r[i]

        si[model.metabolites.get_by_id("h2o_c")] = -1
        si[model.metabolites.get_by_id("h_c")] = 1
        si[model.metabolites.get_by_id("g3pc_c")] = 1

        ri.add_metabolites(si)

    model.add_reactions(r)

def LPLIPA2(model):
    s1 = {model.metabolites.get_by_id("C_BOIMMG_291_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1}

    s2 = {model.metabolites.get_by_id("C_BOIMMG_27923_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1}

    s3 = {model.metabolites.get_by_id("C_BOIMMG_27877_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1}

    s4 = {model.metabolites.get_by_id("C_BOIMMG_27943_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1}

    s5 = {model.metabolites.get_by_id("C_BOIMMG_27925_c"): -1,
          model.metabolites.get_by_id("ttdcea_c"): 1}

    r1 = Reaction("LPLIPA2_1")
    r2 = Reaction("LPLIPA2_2")
    r3 = Reaction("LPLIPA2_3")
    r4 = Reaction("LPLIPA2_4")
    r5 = Reaction("LPLIPA2_5")

    s = [s1, s2, s3, s4, s5]
    r = [r1, r2, r3, r4, r5]

    for i in range(len(r)):
        si = s[i]

        ri = r[i]

        si[model.metabolites.get_by_id("h2o_c")] = -1
        si[model.metabolites.get_by_id("h_c")] = 1
        si[model.metabolites.get_by_id("g3pe_c")] = 1

        ri.add_metabolites(si)

    model.add_reactions(r)

def PASYN_EC(model):
    s1 = {model.metabolites.get_by_id("palmACP_c"): -1,
          model.metabolites.get_by_id("C_BOIMMG_423_c"): 1}

    s2 = {model.metabolites.get_by_id("myrsACP_c"): -1,
          model.metabolites.get_by_id("C_BOIMMG_38435_c"): 1}

    s3 = {model.metabolites.get_by_id("hdeACP_c"): -1,
          model.metabolites.get_by_id("C_BOIMMG_38416_c"): 1}

    s4 = {model.metabolites.get_by_id("octeACP_c"): -1,
          model.metabolites.get_by_id("C_BOIMMG_38371_c"): 1}

    s5 = {model.metabolites.get_by_id("tdeACP_c"): -1,
          model.metabolites.get_by_id("C_BOIMMG_38417_c"): 1}

    r1 = Reaction("reaction_1")
    r2 = Reaction("reaction_2")
    r3 = Reaction("reaction_3")
    r4 = Reaction("reaction_4")
    r5 = Reaction("reaction_5")

    s = [s1, s2, s3, s4, s5]
    r = [r1, r2, r3, r4, r5]

    for i in range(len(r)):
        si = s[i]

        ri = r[i]

        si[model.metabolites.get_by_id("glyc3p_c")] = -1
        si[model.metabolites.get_by_id("ACP_c")] = 2

        ri.add_metabolites(si)

    model.add_reactions(r)

def redundant_granulator_ecoli_all():
    start = time.time()
    model = map_iJR904()

    components = ["cpd00214", "cpd01080", "cpd03847", "cpd05274"]

    solver = RepresentationRedundantCaseSolver(model, "BiGG")

    solver.swap_from_generic(["cpd22513", "cpd15649"], components, False,True)

    cobra.io.write_sbml_model(solver.model, "enhanced_model_ecoli_all.xml")
    cobra.io.save_json_model(solver.model, "enhanced_model_ecoli_all.json")

    end=time.time()

    print("time")
    print(end-start)

if __name__ == "__main__":
    # simulate_ecoli_generalize()
    # simulate_ecoli_granulate()
    # redundant_granulator_ecoli()
    # redundant_granulator_ecoli_without_components()

    # redundant_granulator_ecoli_without_components_isa_reactions_for_all()
    # redundant_granulator_ecoli_without_components()
    gap_fill_ecoli_model()
    # redundant_granulator_ecoli_all()
    # simulate_ecoli_granulate_boimmg()
    # simulate_sacc_generalize()
    # simulate_sacc_granulate()
    # convert_into_json()