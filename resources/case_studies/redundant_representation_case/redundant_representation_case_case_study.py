import cobra
from cobra import Reaction

from boimmgpy.representation_changers import LipidGranulator


def redundant_granulator_yeast():
    model = cobra.io.read_sbml_model("mapped_yeast.xml")
    components = ["C08362", "C00712", "C00249", "C01530"]

    solver = LipidGranulator(model, "KEGG")
    solver.map_model()
    # solver.load_maps("../../dumps/")

    solver.swap_from_generic(["cpd22513", "C00422"], components, False, sources=["LIPID MAPS", "SwissLipids"])
    solver.generate_isa_reactions()

    solver.write_reports("new_report.txt")
    cobra.io.write_sbml_model(solver.model, "granulated_yeast.xml")


def redundant_granulator_ecoli():
    # model = map_iJR904()
    model = cobra.io.read_sbml_model("iJR904_mapped.xml")
    components = ["cpd00214", "cpd03847", "cpd05274", "cpd25615", "cpd05237"]

    solver = LipidGranulator(model, "BiGG")
    solver.map_model()

    solver.swap_from_generic(["cpd22513", "cpd15649"], components, False, sources=["LIPID MAPS"])
    solver.generate_isa_reactions()

    solver.write_reports("new_report.txt")
    cobra.io.write_sbml_model(solver.model, "granulated_iJR904.xml")

    PASYN_EC(model)
    LPLIPA2(model)
    LPLIPA1(model)
    LPLIPA3(model)
    PLIPA1(model)

    cobra.io.write_sbml_model(solver.model, "granulated_gap_filled_iJR904.xml")


def new_test_server():
    model = cobra.io.read_sbml_model("new_model.xml")
    PASYN_EC(model)
    LPLIPA2(model)
    LPLIPA1(model)
    LPLIPA3(model)
    PLIPA1(model)
    cobra.io.write_sbml_model(model, "new_model2.xml")


def PLIPA3(model):
    s1 = {model.metabolites.get_by_id("BMGC427_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1,
          model.metabolites.get_by_id("BMGC314_c"): 1}

    s2 = {model.metabolites.get_by_id("BMGC7454_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1,
          model.metabolites.get_by_id("BMGC16363_c"): 1}

    s3 = {model.metabolites.get_by_id("BMGC7408_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1,
          model.metabolites.get_by_id("BMGC16317_c"): 1}

    s4 = {model.metabolites.get_by_id("BMGC7474_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1,
          model.metabolites.get_by_id("BMGC452_c"): 1}

    s5 = {model.metabolites.get_by_id("BMGC7456_c"): -1,
          model.metabolites.get_by_id("ttdcea_c"): 1,
          model.metabolites.get_by_id("BMGC16364_c"): 1}

    r1 = Reaction("PLIPA3_1")
    r2 = Reaction("PLIPA3_2")
    r3 = Reaction("PLIPA3_3")
    r4 = Reaction("PLIPA3_4")
    r5 = Reaction("PLIPA3_5")

    s = [s1, s2, s3, s4, s5]
    r = [r1, r2, r3, r4, r5]

    for i in range(len(r)):
        si = s[i]

        ri = r[i]

        si[model.metabolites.get_by_id("h2o_c")] = -1
        si[model.metabolites.get_by_id("h_c")] = 1

        ri.add_metabolites(si)

    model.add_reactions(r)


def PLIPA1(model):
    s1 = {model.metabolites.get_by_id("BMGC427_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1,
          model.metabolites.get_by_id("BMGC314_c"): 1}

    s2 = {model.metabolites.get_by_id("BMGC7454_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1,
          model.metabolites.get_by_id("BMGC16363_c"): 1}

    s3 = {model.metabolites.get_by_id("BMGC7408_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1,
          model.metabolites.get_by_id("BMGC16317_c"): 1}

    s4 = {model.metabolites.get_by_id("BMGC7474_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1,
          model.metabolites.get_by_id("BMGC452_c"): 1}

    s5 = {model.metabolites.get_by_id("BMGC7456_c"): -1,
          model.metabolites.get_by_id("ttdcea_c"): 1,
          model.metabolites.get_by_id("BMGC16364_c"): 1}

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
    s1 = {model.metabolites.get_by_id("BMGC314_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1}

    s2 = {model.metabolites.get_by_id("BMGC16363_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1}

    s3 = {model.metabolites.get_by_id("BMGC16317_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1}

    s4 = {model.metabolites.get_by_id("BMGC452_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1}

    s5 = {model.metabolites.get_by_id("BMGC16364_c"): -1,
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
    s1 = {model.metabolites.get_by_id("BMGC223_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1}

    s2 = {model.metabolites.get_by_id("BMGC6942_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1}

    s3 = {model.metabolites.get_by_id("BMGC6896_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1}

    s4 = {model.metabolites.get_by_id("BMGC319_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1}

    s5 = {model.metabolites.get_by_id("BMGC6943_c"): -1,
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
    s1 = {model.metabolites.get_by_id("BMGC291_c"): -1,
          model.metabolites.get_by_id("hdca_c"): 1}

    s2 = {model.metabolites.get_by_id("BMGC27923_c"): -1,
          model.metabolites.get_by_id("hdcea_c"): 1}

    s3 = {model.metabolites.get_by_id("BMGC27877_c"): -1,
          model.metabolites.get_by_id("ocdcea_c"): 1}

    s4 = {model.metabolites.get_by_id("BMGC27943_c"): -1,
          model.metabolites.get_by_id("ttdca_c"): 1}

    s5 = {model.metabolites.get_by_id("BMGC27925_c"): -1,
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
          model.metabolites.get_by_id("BMGC423_c"): 1}

    s2 = {model.metabolites.get_by_id("myrsACP_c"): -1,
          model.metabolites.get_by_id("BMGC38435_c"): 1}

    s3 = {model.metabolites.get_by_id("hdeACP_c"): -1,
          model.metabolites.get_by_id("BMGC38416_c"): 1}

    s4 = {model.metabolites.get_by_id("octeACP_c"): -1,
          model.metabolites.get_by_id("BMGC38371_c"): 1}

    s5 = {model.metabolites.get_by_id("tdeACP_c"): -1,
          model.metabolites.get_by_id("BMGC38417_c"): 1}

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


if __name__ == "__main__":
    redundant_granulator_yeast()
    # redundant_granulator_ecoli()
    # model = cobra.io.read_sbml_model("granulated_gap_filled_iJR904.xml")
    #
    # reactions = []
    # reaction = model.reactions.get_by_id("PASYN_EC")
    # reactions.append(reaction)
    # reaction = model.reactions.get_by_id("LPLIPA2")
    # reactions.append(reaction)
    #
    # reaction = model.reactions.get_by_id("LPLIPA3")
    # reactions.append(reaction)
    #
    # reaction = model.reactions.get_by_id("LPLIPA1")
    # reactions.append(reaction)
    #
    # model.remove_reactions(reactions)

    # reaction = model.reactions.get_by_id("PLIPA1")
    # reactions.append(reaction)

    # cobra.io.save_json_model(model,"new_model_2.json")

    # model.solver = 'cplex'

    # gp = GapFiller(model, exchange_reactions=True, demand_reactions=True, integer_threshold=1e-1000)
    # solution = gp.fill(1)
    #
    # print(solution)
    #
    # model.add_reactions(solution[0])

    # print(model.optimize().objective_value)
    # cobra.io.write_sbml_model(model,"new_model2.xml")

    # new_test_server()
