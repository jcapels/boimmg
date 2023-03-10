import cobra

from boimmgpy import definitions
from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.representation_changers import CofactorSwapper


def sacc_to_9():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iMM904.xml")

    rep_case_solver = CofactorSwapper(model, "ModelSEED")
    rep_case_solver.map_model()
    # rep_case_solver.dump_maps(ROOT_DIR + "/dumps/")
    # rep_case_solver.load_maps(ROOT_DIR + "/dumps/")
    rep_case_solver.swap_compound("cpd15290", "679986")
    rep_case_solver.write_report("test_report.txt")
    cobra.io.write_sbml_model(rep_case_solver.model, "iMM904_9.xml")


def sacc_to_rest():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iMM904.xml")

    rest_accessor = CompoundsDBAccessor()
    rep_case_solver = CofactorSwapper(model, "ModelSEED", rest_accessor)
    rep_case_solver.map_model()
    # rep_case_solver.dump_maps(ROOT_DIR + "/dumps/")
    # rep_case_solver.load_maps(ROOT_DIR + "/dumps/")
    rep_case_solver.swap_compound("cpd15290", "679986")
    rep_case_solver.write_report("test_report.txt")
    cobra.io.write_sbml_model(rep_case_solver.model, "iMM904_9.xml")


def simulate_ecoli_generalize():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iML1515.xml")
    rep_case_solver = CofactorSwapper(model, "BiGG")
    rep_case_solver.swap_compound("cpd15560", "cpd11669")
    rep_case_solver.swap_compound("cpd15500", "cpd11451")
    cobra.io.write_sbml_model(rep_case_solver.model, "iML1515_generalized.xml")


def simulate_ecoli_granulate():
    model = cobra.io.read_sbml_model("iML1515_generalized.xml")
    rep_case_solver = CofactorSwapper(model, "BiGG")
    rep_case_solver.swap_compound("cpd11669", "cpd15290")
    rep_case_solver.swap_compound("cpd11451", "cpd15500")
    cobra.io.write_sbml_model(rep_case_solver.model, "iML1515_granulated.xml")


def simulate_ecoli_granulate_boimmg():
    model = cobra.io.read_sbml_model("iML1515_generalized.xml")
    rep_case_solver = CofactorSwapper(model, "BiGG")
    rep_case_solver.swap_compound("cpd11669", "C_BOIMMG_749196")
    cobra.io.write_sbml_model(rep_case_solver.model, "iML1515_granulated_ubi5.xml")


def simulate_sacc_generalize():
    accessor = CompoundsDBAccessor()

    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iMM904.xml")
    rep_case_solver = CofactorSwapper(model, "BiGG", accessor)
    rep_case_solver.swap_compound("cpd15290", "cpd11669")
    cobra.io.write_sbml_model(rep_case_solver.model, "iMM904_generalized.xml")


def simulate_sacc_granulate():
    model = cobra.io.read_sbml_model("iMM904_generalized.xml")
    rep_case_solver = CofactorSwapper(model, "BiGG")
    rep_case_solver.swap_compound("cpd11669", "cpd15290")
    cobra.io.write_sbml_model(rep_case_solver.model, "iMM904_granulated.xml")


if __name__ == "__main__":
    # sacc_to_rest()
    # simulate_ecoli_generalize()
    # simulate_ecoli_granulate()
    # simulate_ecoli_granulate_boimmg()
    simulate_sacc_generalize()
    # simulate_sacc_granulate()
