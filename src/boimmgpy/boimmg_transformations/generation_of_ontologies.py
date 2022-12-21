from src.boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from src.boimmgpy.ontology_generators.ontology_generator import OntologyGeneratorDA, TransformationsHandler


def run_ontology_generation_for_lmsd_ms_compounds():

    coiso = OntologyGeneratorDA()

    core = "C(=O)O"
    core_fatty_acids = \
        "C(=O)O"
    core_alcohol = "*CO"

    print("Cardiolipins")
    coiso("Cardiolipin-Biosynthesis",[core],"transformations_cardiolipin",["PWY-7509"],all=False)

    print("Triacylglycerides")
    coiso("TRIGLSYN-PWY", [core_fatty_acids,core_alcohol], "transformations_TRIGLSYN-PWY",
      ["PWY-7277","PWY-7835"],["RXN-12383","RXN-1641"],all=False)


    print("Phosphatidylserine")
    coiso("Phosphatidylserine-Biosynthesis",[core],"transformations_Phosphatidylserine")

    print("Phosphatidylethanolamine")
    coiso("PhosphatidylethanolamineBiosynthesis", [core], "transformations_PhosphatidylethanolamineBiosynthesis",
      ["PWY-7509", "PWY-7409"],all=False)

    print("Phosphatidylcholine")
    coiso("PhosphatidylcholineBiosynthesis",[core],"transformations_PhosphatidylcholineBiosynthesis",
          ["PWY-6826","PWY3O-450"],all=False)

    print("Phosphatidylglycerol")
    coiso("PhosphatidylglycerolBiosynthesis", [core], "transformations_Phosphatidylglycerol",
      ["PWY4FS-7"],all=False)

    print("plasmalogen")
    coiso("PWY-7782", [core_fatty_acids, core_alcohol], "transformations_PWY-7782", all=False)

    print("3-phosphoinositide")
    coiso("PWY-6352", [core_fatty_acids, core_alcohol], "transformations_PWY-6352",
          all=False)

    print("CDP-diacylglycerol")
    coiso("CDP-diacylglycerol-Biosynthesis", [core_fatty_acids, core_alcohol], "transformations_CDP-diacylglycerol", all=False)

def generate_transformation():
    transformations_cardiolipin = TransformationsHandler()
    transformations_cardiolipin.choose_transformations_from_pathway("Phosphatidylserine-Biosynthesis")
    transformations_cardiolipin.save_transformations("transformations_Phosphatidylserine")

def integrate_quinones():
    accessor = CompoundsDBAccessor()
    ele_quinone = accessor.get_predecessors_by_ont_id_rel_type(749300,"is_a")
    ele_quinol = accessor.get_predecessors_by_ont_id_rel_type(749301,"is_a")

    for ele_quinon in ele_quinone:

        ele_quinone_precursors = accessor.get_all_predecessors_by_ont_id_rel_type(ele_quinon,"precursor_of")

        for pre in ele_quinone_precursors:
            accessor.add_relationship(pre,321415,"is_a")

    for ele_quino in ele_quinol:

        ele_quinone_precursors = accessor.get_all_predecessors_by_ont_id_rel_type(ele_quino,"precursor_of")

        for pre in ele_quinone_precursors:
            accessor.add_relationship(pre,321415,"is_a")

def run_ontology_generation_sphingolipids():
    coiso = OntologyGeneratorDA()

    core = "C(=O)O"
    core_fatty_acids = \
        "C(=O)O"
    core_alcohol = "*CO"

    coiso("Sphingolipid-Biosynthesis", [core_fatty_acids, core_alcohol], "transformations_Sphingolipid-Biosynthesis")

if __name__ == "__main__":
    run_ontology_generation_sphingolipids()
    # run_ontology_generation_for_lmsd_ms_compounds()
    # integrate_quinones()