from src.boimmgpy import definitions
from src.boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from src.boimmgpy.service.network_handlers.pathway_handler import printProgressBar
from src.boimmgpy.utilities import file_utilities
from src.boimmgpy import COMPOUNDS_ANNOTATION_CONFIGS_PATH


def check_if_reaction_exists(reaction_compounds, compounds_in_model_db_ids):
    """
    This method checks whether a given reaction is equal to another in the model. This is performed by removing and
    adding hydrogens to the one of the reactions.

    :param cobrapy.Reaction reaction: model reaction
    :param ModelSeedReaction reaction_container:
    :param list reaction_compounds: modelseed ids of the reactants and products
    :param list compounds_in_model_db_ids: model seed ids from the compounds in the model reaction
    :return boolean:
    """

    reaction_without_hydrogen = reaction_compounds.copy()
    reaction_with_hydrogen = reaction_compounds.copy()
    modelseed_hydrogen = "cpd00067"

    if modelseed_hydrogen in reaction_without_hydrogen:
        reaction_without_hydrogen.remove(modelseed_hydrogen)

    elif modelseed_hydrogen not in reaction_with_hydrogen:
        reaction_with_hydrogen.append(modelseed_hydrogen)

    if sorted(reaction_without_hydrogen) == sorted(compounds_in_model_db_ids) \
            or sorted(reaction_with_hydrogen) == sorted(compounds_in_model_db_ids) \
            or sorted(reaction_compounds) == sorted(compounds_in_model_db_ids):

        return True

    else:
        return False


def compare_models_reactions(model1, model2, model_database, compoundsAnnotationConfigs, compoundsIDConverter,flag=True):

    rxns1 = model1.reactions
    rxns2 = list(model2.reactions)
    reactions = rxns2.copy()
    diff_reactions_list = []
    equal_reactions_list2 = []
    equal_reactions_list1 = []
    equal_reactions = 0
    t=0
    for r1 in rxns1:
        t+=1
        printProgressBar(t,len(rxns1))
        modelseedids1 = None
        biggids1 = None
        if "seed.reaction" in r1.annotation.keys():

            modelseedids1 = r1.annotation.get("seed.reaction")
            if type(modelseedids1) == str:
                modelseedids1 = [modelseedids1]

        if "bigg.reaction" in r1.annotation.keys():
            biggids1 = r1.annotation.get("bigg.reaction")
            if type(biggids1) == str:
                biggids1 = [biggids1]

        found = False
        i = 0
        len_1 = len(r1.metabolites)

        while not found and i < len(reactions):

            r2 = reactions[i]

            len_2 = len(r2.metabolites)

            if len_1-len_2 > 1 or len_1-len_2 < -1:
                pass

            else:

                if r1.id == r2.id:
                    found = True

                if flag:
                    if not found and "seed.reaction" in r2.annotation.keys() and modelseedids1:
                        found = check_by_model_seed_reaction(r2, modelseedids1)

                    elif not found and "bigg.reaction" in r2.annotation.keys() and biggids1:
                        found = check_by_bigg_reactions(r2, biggids1)

                    if not found and r1.annotation.get("sbo") != "SBO:0000629" and r2.annotation.get("sbo") != "SBO:0000629":
                        equal = \
                            check_reaction_by_compounds(r1, r2, model_database, compoundsAnnotationConfigs,
                                                        compoundsIDConverter)

                        if equal:

                            found = True

                    if found:
                        reactions.pop(i)

            i += 1

        if found:
            equal_reactions_list2.append(r2.id)
            equal_reactions_list1.append(r1.id)
            equal_reactions += 1

        else:
            diff_reactions_list.append(r1.id)



    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    # First way to call the 2 group Venn diagram:
    venn = venn2(subsets=(len(rxns1) - equal_reactions, len(rxns2) - equal_reactions, equal_reactions),
          set_labels=(model1.id, model2.id))

    for text in venn.set_labels:
        text.set_fontsize(14)

    for text in venn.subset_labels:
        text.set_fontsize(16)
    plt.show()

    return equal_reactions_list1,equal_reactions_list2


def check_by_model_seed_reaction(r2, modelseedids1):
    modelseedids2 = r2.annotation.get("seed.reaction")
    if type(modelseedids2) == str:
        modelseedids2 = [modelseedids2]
    j = 0
    while j < len(modelseedids1):
        l = 0
        while l < len(modelseedids2):
            if modelseedids2[l] == modelseedids1[j]:
                return True
            else:
                l += 1
        j += 1

    return False


def check_reaction_by_compounds(r1, r2, model_database, compoundsAnnotationConfigs, compoundsIDConverter):
    ids1 = get_database_ids_compounds_in_model(r1.metabolites, model_database, compoundsAnnotationConfigs,
                                               compoundsIDConverter)

    ids2 = get_database_ids_compounds_in_model(r2.metabolites, model_database, compoundsAnnotationConfigs,
                                               compoundsIDConverter)

    for id1 in ids1:
        for id2 in ids2:


            equal = check_if_reaction_exists(id2, id1)

            if equal:
                return True

    return False


def get_database_ids_compounds_in_model(compounds_in_model, model_database, compoundsAnnotationConfigs,
                                        compoundsIDConverter):
    """
    This method searches in a given database the ids of the compounds in the model

    :param dict compounds_in_model: model compounds assigned to a model reaction

    :returns list: list with database ids of the metabolites assigned to a model reaction

    """

    res = []
    for compound in compounds_in_model.keys():
        modelseed_id_annotation = compound.annotation.get("seed.compound")

        modelseed_id = []

        if modelseed_id_annotation and type(modelseed_id_annotation) == str:
            modelseed_id.append(modelseed_id_annotation)

        elif modelseed_id_annotation:
            modelseed_id.extend(modelseed_id_annotation)

        boimmg = False
        if "boimmg.compound" in compound.annotation.keys():
            boimmg_id = compound.annotation["boimmg.compound"]
            boimmg = True

        elif model_database != "ModelSEED":
            db_id = compound.annotation.get(compoundsAnnotationConfigs[model_database])
            modelseed_ids = []
            if type(db_id) == list:

                i = 0
                found = False

                while not found and i < len(db_id):
                    modelseed_ids = compoundsIDConverter.convert_dbID_into_modelSeedId(model_database,
                                                                                      db_id[i])

                    if modelseed_ids:
                        found = True
            else:
                modelseed_ids = compoundsIDConverter.convert_dbID_into_modelSeedId(model_database,
                                                                                  db_id)

            if modelseed_ids:
                if modelseed_ids[0] not in modelseed_id:
                    modelseed_id.append(modelseed_ids[0])



        if not modelseed_id and not boimmg:
            return []

        else:
            if not res:
                if boimmg:
                    res.append([boimmg_id])

                elif modelseed_id:
                    for id in modelseed_id:
                        res.append([id])

            else:

                if boimmg:
                    if modelseed_id:
                        modelseed_id.append(boimmg_id)

                    else:
                        previous = res.copy()
                        res = []
                        for reaction in previous:
                            new_reaction = reaction.copy()
                            new_reaction.append(boimmg_id)
                            res.append(new_reaction)

                if modelseed_id:
                    if len(modelseed_id) == 1:
                        for reaction in res:
                            reaction.append(modelseed_id[0])
                    else:
                        previous = res.copy()
                        res = []
                        for id in modelseed_id:
                            for reaction in previous:
                                new_reaction = reaction.copy()
                                new_reaction.append(id)
                                res.append(new_reaction)


        # elif res:
        #     for reaction in res:
        #         reaction.append(modelseed_id)
        # else:
        #     res.append([db_id])

    final_res = []
    for reaction in res:
        final_res.append(sorted(reaction))
    return final_res


def check_by_bigg_reactions(r2, biggids1):
    biggids2 = r2.annotation.get("bigg.reaction")
    if type(biggids2) == str:
        biggids2 = [biggids2]
    j = 0
    while j < len(biggids1):
        l = 0
        while l < len(biggids2):
            if biggids1[j] == biggids2[l]:
                return True
            else:
                l += 1
        j += 1

    return False


def check_by_model_seed_compounds(r2, modelseedids1):
    modelseedids2 = r2.annotation.get("seed.compound")
    if type(modelseedids2) == str:
        modelseedids2 = [modelseedids2]
    j = 0
    while j < len(modelseedids1):
        l = 0
        while l < len(modelseedids2):
            if modelseedids2[l] == modelseedids1[j]:
                return True,modelseedids2[l]
            else:
                l += 1
        j += 1

    return False,None


def check_by_bigg_compounds(r2, biggids1):
    biggids2 = r2.annotation.get("bigg.metabolite")
    if type(biggids2) == str:
        biggids2 = [biggids2]
    j = 0
    while j < len(biggids1):
        l = 0
        while l < len(biggids2):
            if biggids1[j] == biggids2[l]:
                return True,biggids1[j]
            else:
                l += 1
        j += 1

    return False,None


def compare_models_metabolites(model1, model2, met_print =False):
    mtx1 = model1.metabolites
    mtx2 = list(model2.metabolites)

    equal_reactions = 0
    model_seed_seen = []
    bigg_seen = []
    boimmg_seen = []
    ids_seen = []
    found_metabolites = []
    for r1 in mtx1:
        modelseedids1 = None
        biggids1 = None

        boimmg = False
        if "boimmg.compound" in r1.annotation.keys():
            boimmg_id = r1.annotation.get("boimmg.compound")
            boimmg = True


        if "seed.compound" in r1.annotation.keys():

            modelseedids1 = r1.annotation.get("seed.compound")
            if type(modelseedids1) == str:
                modelseedids1 = [modelseedids1]

        if "bigg.metabolite" in r1.annotation.keys():
            biggids1 = r1.annotation.get("bigg.metabolite")
            if type(biggids1) == str:
                biggids1 = [biggids1]

        found = False
        i = 0
        while not found and i < len(mtx2):
            r2 = mtx2[i]
            if r1.id == r2.id and r2.id not in ids_seen:
                ids_seen.append(r1.id)
                found = True


            # if not found and "seed.compound" in r2.annotation.keys() and modelseedids1 and r1.compartment == r2.compartment:
            #     found = check_by_model_seed_compounds(r2, modelseedids1)
            #
            # if not found and "bigg.metabolite" in r2.annotation.keys() and biggids1 and r1.compartment == r2.compartment:
            #     found = check_by_bigg_compounds(r2, biggids1)
            #
            # if not found and boimmg and "boimmg.compound" in r2.annotation.keys():
            #     boimmg_compound2 = r2.annotation.get("boimmg.compound")
            #     if boimmg_id == boimmg_compound2 and r1.compartment == r2.compartment:
            #         found = True

            if "seed.compound" in r2.annotation.keys() and modelseedids1:
                found, ms_id = check_by_model_seed_compounds(r2, modelseedids1)
                if found:
                    if ms_id not in model_seed_seen:
                        model_seed_seen.append(ms_id)
                    else:
                        found = False

            if "bigg.metabolite" in r2.annotation.keys() and biggids1:
                found,bigg_id = check_by_bigg_compounds(r2, biggids1)
                if found:
                    if bigg_id not in bigg_seen:
                        bigg_seen.append(bigg_id)
                    else:
                        found = False

            if boimmg and "boimmg.compound" in r2.annotation.keys():
                boimmg_compound2 = r2.annotation.get("boimmg.compound")
                if boimmg_id == boimmg_compound2 and boimmg_id not in boimmg_seen:
                    boimmg_seen.append(boimmg_id)
                    found = True

            i += 1

        if found:
            found_metabolites.append(r2.id)
            equal_reactions += 1

        elif met_print:

            print(r1)




    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2


    # First way to call the 2 group Venn diagram:
    venn = venn2(subsets=(len(mtx1) - equal_reactions, len(mtx2) - equal_reactions, equal_reactions),
          set_labels=(model1.id, model2.id), set_colors=("b","g"),alpha=0.4)

    for text in venn.set_labels:
        text.set_fontsize(14)

    for text in venn.subset_labels:
        text.set_fontsize(16)
    plt.show()

    return found_metabolites




def compare_models(model1, model2, model3, model_database, compoundsAnnotationConfigs,
                                                 compoundsIDConverter):
    model3 = cobra.io.read_sbml_model(model3)
    model2 = cobra.io.read_sbml_model(
        model2)
    model1 = cobra.io.read_sbml_model(model1)

    found_reactions_1_1, found_reactions_1_2 = compare_models_reactions(model3, model1, model_database,
                                                                        compoundsAnnotationConfigs,
                                                                        compoundsIDConverter)
    found_1 = compare_models_metabolites(model3, model1)

    found_reactions_2_2, found_reactions_2_2 = compare_models_reactions(model2, model1, model_database,
                                                                        compoundsAnnotationConfigs,
                                                                        compoundsIDConverter)
    found_2 = compare_models_metabolites(model2, model1)

    print("--------------------- Metabolites ---------------------------")
    new_list = []
    for met2 in found_2:
        if met2 not in found_1:
            new_list.append(met2)

    for met2 in sorted(new_list):
        print(met2)

    print("--------------------- Reactions ---------------------------")
    for reaction2 in found_reactions_2_2:
        if reaction2 not in found_reactions_1_2:
            print(reaction2)

def compare_lipids_ecoli_granulated(model_database, compoundsAnnotationConfigs,
                                                 compoundsIDConverter):

    model3 = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iJR904.xml")
    model2 = cobra.io.read_sbml_model(
        definitions.ROOT_DIR + "/case_studies/enhanced_model_ecoli_without_components_gap_filled.xml")
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iAF1260b_mapped.xml")




    found_reactions_1_1,found_reactions_1_2 = compare_models_reactions(model3, model, model_database, compoundsAnnotationConfigs,
                                                 compoundsIDConverter)
    found_1 = compare_models_metabolites(model3, model)

    found_reactions_2_2,found_reactions_2_2 = compare_models_reactions(model2, model, model_database, compoundsAnnotationConfigs,
                                                 compoundsIDConverter)
    found_2 = compare_models_metabolites(model2, model)

    print("--------------------- Metabolites ---------------------------")
    new_list = []
    for met2 in found_2:
        if met2 not in found_1:
            new_list.append(met2)

    for met2 in sorted(new_list):
        print(met2)

    print("--------------------- Reactions ---------------------------")
    for reaction2 in found_reactions_2_2:
        if reaction2 not in found_reactions_1_2:
            print(reaction2)



def granulated_ecoli_models(model_database, compoundsAnnotationConfigs, compoundsIDConverter):
    import cobra
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iMM904.xml")
    model2 = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/case_studies/iML1515_granulated.xml")
    model3 = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iML1515.xml")

    # compare_models_metabolites(model2, model3, met_print = True)

    # print(len(model2.reactions))
    # print(len(model3.reactions))



    found_reactions_1 = compare_models_reactions(model, model3, model_database, compoundsAnnotationConfigs, compoundsIDConverter)
    found_1 =  compare_models_metabolites(model, model3)

    found_reactions_2 = compare_models_reactions(model,model2 , model_database, compoundsAnnotationConfigs, compoundsIDConverter)
    found_2 = compare_models_metabolites(model, model2)

    print(len(found_1))
    print(len(found_2))
    for met2 in found_2:
        if met2 not in found_1:
            print(met2)

    for reaction2 in found_reactions_2:
        if reaction2 not in found_reactions_1:
            print(reaction2)



if __name__ == "__main__":
    import cobra

    compoundsIDConverter = CompoundsIDConverter()
    compoundsAnnotationConfigs = file_utilities.read_conf_file(
        COMPOUNDS_ANNOTATION_CONFIGS_PATH)
    model_database = "BiGG"

    model1 = definitions.ROOT_DIR + "/case_studies/redundant_representation_case/iML1515_mapped.xml"
    model3 = definitions.ROOT_DIR + "/case_studies/redundant_representation_case/iJR904_mapped.xml"
    model2 = definitions.ROOT_DIR + "/case_studies/redundant_representation_case/granulated_gap_filled_iJR904.xml"

    # compare_lipids_ecoli_granulated(model_database, compoundsAnnotationConfigs, compoundsIDConverter)
    compare_models(model1,model2,model3,model_database, compoundsAnnotationConfigs, compoundsIDConverter)
