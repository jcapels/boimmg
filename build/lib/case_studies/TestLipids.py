import time
from unittest import TestCase, TestLoader, TextTestRunner

import cobra
from cobra import Model
from cobra.flux_analysis.gapfilling import GapFiller

from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.boimmg.RepresentationRevisors.RepresentationSimpleCaseSolver import SimpleCaseSolver
from boimmgpy.boimmg.RepresentationRevisors.RepresentationRedundantCaseSolver import RepresentationRedundantCaseSolver
from boimmgpy.case_studies.redundant_representation_case.map_ecoli_model import map_iJR904


class TestQuinones(TestCase):


    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print('%s: %.9f' % (self.id(), t))

    def test_simulate_ecoli_ubiquinone10(self):

        model = cobra.io.read_sbml_model("C:/Users/Biosystems/Desktop/thesis_jcapela/BOIMMGpy/models/iML1515.xml")
        gapfiller = SimpleCaseSolver(model, "BiGG")
        gapfiller.swap_and_gap_fill("cpd08232")
        cobra.io.write_sbml_model(gapfiller.model,"iML1515_ubq10.xml")

    def test_simulate_ecoli_generalize(self):
        model = cobra.io.read_sbml_model("C:/Users/Biosystems/Desktop/thesis_jcapela/BOIMMGpy/models/iML1515.xml")
        gapfiller = SimpleCaseSolver(model, "BiGG")
        gapfiller.swap_compound("cpd15560","cpd11669")
        cobra.io.write_sbml_model(gapfiller.model, "simple_representation_case/iML1515_generalized.xml")

    def test_database_access(self):

        accessor = CompoundsDBAccessor()
        accessor.get_successors_by_ont_id_rel_type(3,"is_a")

class TestLipids(TestCase):


    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print('%s: %.9f' % (self.id(), t))

    def get_targets_to_replace_by_type(self, target: int):
        compounds_ontology = CompoundsDBAccessor()

        res={}
        precursors = compounds_ontology.get_all_predecessors_by_ont_id_rel_type(target,"precursor_of")

        for precursor in precursors:
            parent = compounds_ontology.get_successors_by_ont_id_rel_type(precursor,"is_a")

            if parent[0] in res:
                res.get(parent[0]).append(precursor)
            else:
                res[parent[0]] = [precursor]
        return res

    def test_precursors(self):
        compounds_ontology = CompoundsDBAccessor()
        all_precursors = compounds_ontology.get_all_predecessors_by_ont_id_rel_type(322790,"precursor_of")
        print(all_precursors)


    def test_simulate_ecoli_ubiquinone10(self):

        model = cobra.io.read_sbml_model("C:/Users/Biosystems/Desktop/thesis_jcapela/BOIMMGpy/models/iML1515.xml")
        gapfiller = SimpleCaseSolver(model, "BiGG")
        gapfiller.swap_and_gap_fill(690028)

    def test_components_scraping(self):
        db_acessor = CompoundsDBAccessor()
        components = [200,437,436,129]

        all_new_met = []
        compounds = db_acessor.get_compounds_with_specific_parent_within_set_of_components(12, components)
        for compound in compounds:
            precursors = db_acessor.get_all_predecessors_by_ont_id_rel_type(compound,"precursor_of")
            all_new_met.extend(precursors)

        # res = []
        # for new_met in all_new_met:
        #     if new_met not in res:
        #         res.append(new_met)

        print(len(all_new_met))
        # print(res)

    def test_identification_of_biosynthesis_pathway_ecoli(self):
        model = map_iJR904()

        solver = RepresentationRedundantCaseSolver(model,"BiGG")

        reactions,other = solver.identify_biosynthesis_pathway(12)

        assert len(reactions) == 10
        # solver.swap_from_generic(12,[200,437,436,129])

    def test_dict_with_predecessors(self):
        model = map_iJR904()

        solver = RepresentationRedundantCaseSolver(model, "BiGG")
        res = solver.get_targets_to_replace_by_type(300374)

        print()


    def test_phospholipid_granulation_ecoli(self):

        model = map_iJR904()

        components = ["cpd00214", "cpd01080", "cpd03847", "cpd05274"]

        solver = RepresentationRedundantCaseSolver(model, "BiGG")

        solver.swap_from_generic(["cpd22513","cpd15649"], components, True)

        cobra.io.write_sbml_model(solver.model, "../models/enhanced_model_ecoli.xml")
        cobra.io.save_json_model(model, "enhanced_model_ecoli.json")

    def test_components_reactions_solver(self):
        model = map_iJR904()

        solver = RepresentationRedundantCaseSolver(model, "BiGG")

        solver.solve_components_reactions(True)


    def test_optimization(self):

        model = cobra.io.read_sbml_model("../models/enhanced_model_ecoli.xml")

        gapfiller = GapFiller(model, Model(),
                              demand_reactions=True,
                              exchange_reactions=True,
                              integer_threshold=1e-300)

        solution = gapfiller.fill(iterations=1)
        print(solution)
        model.add_reactions(solution[0])

        print(model.optimize().objective_value)

        cobra.io.write_sbml_model(model,"enhanced_model_ecoli2.xml")
        cobra.io.save_json_model(model,"enhanced_model_ecoli.json")

    def test_read_universal_model(self):

        model = cobra.io.read_sbml_model("C:/Users/BioSystems/Desktop/universal_model.xml")

    def test_algae_model(self):
        model = cobra.io.read_sbml_model("../models/cvulgaris.xml")

        model.objective = model.reactions.get_by_id("e_Biomass__cytop")

        r_acyl = model.reactions.get_by_id("e_Fatty_acid__cytop")
        reactants = r_acyl.reactants

        components_list = []
        for reactant in reactants:

            if "kegg.compound" in reactant.annotation:
                compound = reactant.annotation.get("kegg.compound")
                components_list.append(compound)


        ### confirmar se há alguma reação com fatty acid

        solver = RepresentationRedundantCaseSolver(model, "KEGG")

        # components = ["lgnc","arach","C06427","clpnd"]

        solver.swap_from_generic(["C00422","C06037"],components_list,True)

        cobra.io.write_sbml_model(model, "very_nice_algae_model.xml")

if __name__ == '__main__':
    suite = TestLoader().loadTestsFromTestCase()
    TextTestRunner(verbosity=0).run(suite)