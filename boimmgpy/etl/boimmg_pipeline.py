import math
import re
from collections import defaultdict

from biocyc import biocyc
from neo4j import GraphDatabase
from rdkit import Chem
from rdkit.Chem import MolToInchiKey, AllChem
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToSmiles, MolFromSmarts
from sympy import solve

import pandas as pd

from boimmgpy.database.accessors.compounds_database_accessor import CompoundsDBAccessor
from boimmgpy.database.accessors.reactions_database_accessor import ReactionsDBAccessor
from boimmgpy.database.databases_babel import BOIMMGDatabases
from boimmgpy.etl.model_seed.compounds.ModelSeedCompoundsDB import ModelSeedCompoundsDB
from boimmgpy.id_converters.compounds_id_converter import CompoundsIDConverter
from boimmgpy.id_converters.reactions_id_converter import ReactionsIDConverter
from boimmgpy.model_seed.model_seed_compound import ModelSeedCompound
from boimmgpy.model_seed.model_seed_reactions_database import ModelSeedReactionsDB
from boimmgpy.ontology_generators.ontology_generator import TransformationsHandler, OntologyGeneratorDA
from boimmgpy.ontology_generators.side_chain_checker import SidechainChecker
from boimmgpy.service.network_handlers.pathway_handler import PathwayHandler
from boimmgpy.utilities import ChemoUtilities, chemo_utilities
from boimmgpy.utilities.LipidMapsStructureDB import LipidMapsStructureDB
from boimmgpy.utilities.chemo_utilities import neutralise_charges


def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='|', printEnd="\r"):
    """
  Call in a loop to create terminal progress bar
  @params:
      iteration   - Required  : current iteration (Int)
      total       - Required  : total iterations (Int)
      prefix      - Optional  : prefix string (Str)
      suffix      - Optional  : suffix string (Str)
      decimals    - Optional  : positive number of decimals in percent complete (Int)
      length      - Optional  : character length of bar (Int)
      fill        - Optional  : bar fill character (Str)
      printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
  """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end="", flush=True)
    # Print New Line on Complete
    if iteration == total:
        print()


def read_swiss_lipids_database():
    import os
    data = pd.read_csv(os.getcwd() + "/lipids.tsv", header=0, delimiter="\t", encoding="ISO-8859-1")
    modelSeedDB = ModelSeedCompoundsDB()

    idConverter = CompoundsIDConverter()

    header = ["swiss_lipids_id:ID", "name", "smiles", "inchi", "inchikey", "formula", "charge:int", "mass:float",
              "hmdb_id", "chebi_id", "lipidmaps_id", "pubchem_cid", "kegg_id", "bigg_id", "metanetx_id",
              "metacyc_id", "generic", "model_seed_id", ":LABEL"]

    new_data = []

    for i in range(data.shape[0]):
        printProgressBar(i, data.shape[0])

        inchikey = data.loc[:, "InChI key (pH7.3)"].iloc[i]
        if not pd.isna(inchikey):
            inchikey = inchikey.replace("InChIKey=", "")

        else:
            inchikey = None

        smiles = data.loc[:, "SMILES (pH7.3)"].iloc[i]
        level = data.loc[:, "Level"].iloc[i]
        if pd.isna(smiles):
            smiles = ""

        if "*" not in smiles or "Class" in level:

            canonical_smiles = None
            if not pd.isna(smiles):
                mol_smiles = MolFromSmiles(smiles)
                if mol_smiles:
                    canonical_smiles = MolToSmiles(mol_smiles)

            swisslipids_id = data.loc[:, "Lipid ID"].iloc[i]
            name = data.loc[:, "Name"].iloc[i]
            formula = data.loc[:, "Formula (pH7.3)"].iloc[i]
            inchi = data.loc[:, "InChI (pH7.3)"].iloc[i]
            if not pd.isna(inchi) and inchi.split("=")[1] == "none":
                inchi = None
            hmdb_id = data.loc[:, "HMDB"].iloc[i]
            chebi_id = str(data.loc[:, "CHEBI"].iloc[i])
            lipidmaps_id = data.loc[:, "LIPID MAPS"].iloc[i]
            pubchem_cid = str(data.loc[:, "PMID"].iloc[i])
            charge = data.loc[:, "Charge (pH7.3)"].iloc[i]
            mass = data.loc[:, "Mass (pH7.3)"].iloc[i]

            if not pd.isna(hmdb_id):
                hmdb_id = hmdb_id.split("|")
                if hmdb_id:
                    hmdb_id = hmdb_id[0]
                else:
                    hmdb_id = None

            if not pd.isna(chebi_id):
                chebi_id = chebi_id.split("|")
                if chebi_id:
                    chebi_id = chebi_id[0]
                else:
                    chebi_id = None

            if not pd.isna(lipidmaps_id):
                lipidmaps_id = lipidmaps_id.split("|")
                if lipidmaps_id:
                    lipidmaps_id = lipidmaps_id[0]
                else:
                    lipidmaps_id = None

            if not pd.isna(pubchem_cid):
                pubchem_cid = pubchem_cid.split("|")
                if pubchem_cid:
                    pubchem_cid = pubchem_cid[0]
                else:
                    pubchem_cid = None

            generic = False

            model_seed_compound = None
            if pd.isna(inchikey):
                generic = True
                if pd.isna(smiles) and canonical_smiles:
                    model_seed_compound = modelSeedDB.get_compound_by_canonical_smiles(canonical_smiles)

            else:
                model_seed_compound = modelSeedDB.get_compound_by_inchi_key(inchikey)

            if model_seed_compound:

                kegg_id, bigg_id, metanetx_id, metacyc_id = integrate_model_ids(idConverter, model_seed_compound)

                new_line = [swisslipids_id, name, canonical_smiles, inchi, inchikey, formula, charge, mass, hmdb_id,
                            chebi_id, lipidmaps_id, pubchem_cid, kegg_id, bigg_id, metanetx_id, metacyc_id, generic,
                            model_seed_compound.getDbId()]

            else:
                kegg_id = None
                bigg_id = None
                metanetx_id = None
                metacyc_id = None

                new_line = [swisslipids_id, name, canonical_smiles, inchi, inchikey, formula, charge, mass, hmdb_id,
                            chebi_id, lipidmaps_id, pubchem_cid, kegg_id, bigg_id, metanetx_id, metacyc_id, generic,
                            None]

            new_data.append(new_line)

    data_frame = pd.DataFrame(new_data)
    data_frame.to_csv("entities.csv", sep=",", header=header, index=False)


def integrate_model_ids(idConverter, model_seed_compound):
    kegg_ids = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "KEGG")
    bigg_ids = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "BiGG")
    metanetxs = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "metanetx.chemical")
    metacyc = idConverter.convert_modelSeedId_into_other_dbID(model_seed_compound.getDbId(), "MetaCyc")

    if kegg_ids:
        kegg_id = kegg_ids[0]
    else:
        kegg_id = None

    if bigg_ids:
        bigg_id = bigg_ids[0]

    else:
        bigg_id = None

    if metanetxs:
        metanetx_id = metanetxs[0]

    else:
        metanetx_id = None

    if metacyc:
        metacyc_id = metacyc[0]

    else:
        metacyc_id = None

    return kegg_id, bigg_id, metanetx_id, metacyc_id


def create_file_for_relationships():
    header = [":START_ID", ":END_ID", ":TYPE"]
    data_list = []
    with open("components.tsv", "r") as file:
        lines = file.readlines()
        for line in lines:
            line_list = line.split("\t")
            end = line_list[0]
            start = line_list[1].strip()
            type = "component_of"
            new_line = [start, end, type]
            data_list.append(new_line)

    with open("../parents.tsv", "r") as file:
        lines = file.readlines()
        for line in lines:
            line_list = line.split("\t")
            start = line_list[0]
            end = line_list[1].strip()
            type = "is_a"
            new_line = [start, end, type]
            data_list.append(new_line)

    data_frame = pd.DataFrame(data_list)
    data_frame.to_csv("rel.csv", sep=",", header=header, index=False)


############################################### LIPID MAPS database integration #############################################

def scrap_lipid_maps_db(tx):
    lipids_db = LipidMapsStructureDB()
    db = lipids_db.getDatabase()
    i = 0
    accessor = CompoundsDBAccessor()
    for lipid in db:
        printProgressBar(i, len(db))
        lipid_container = db[lipid]

        smiles = lipid_container.getSmiles()
        if pd.isna(smiles):
            canonical_smiles = None
        else:
            try:
                canonical_smiles = MolToSmiles(MolFromSmiles(smiles))
            except:
                canonical_smiles = None
        aliases = lipid_container.getAliases()

        chebi_id = None
        lipid_bank_id = None
        pubchem_id = None
        hmdb_id = None
        kegg_id = None
        name = lipid_container.getName()

        if "KEGG" in aliases.keys():
            kegg_id = aliases.get("KEGG")[0]

        if "HMDB" in aliases.keys():
            hmdb_id = aliases.get("HMDB")[0]

        if "ChEBI" in aliases.keys():
            chebi_id = aliases.get("ChEBI")[0]

        if "LipidBank" in aliases.keys():
            lipid_bank_id = aliases.get("LipidBank")[0]

        if "PubChem" in aliases.keys():
            pubchem_id = aliases.get("PubChem")[0]

        inchikey = lipid_container.getInchiKey()

        inchikey_wlast = None
        if type(inchikey) == str:
            inchikey_wlast = inchikey[:-1]

        if inchikey_wlast:
            res = accessor.get_compound_by_inchikey(inchikey)

            if not res:

                with tx.session() as session:
                    session.run("MERGE (c:Compound { lipidmaps_id: $lipid_maps_id}) "
                                "ON CREATE SET c.lipid_maps_id = $lipid_maps_id, "
                                "c.smiles = $smiles, c.generic = False, c.inchikey = $inchikey, "
                                "c.formula = $formula, c.charge = 0, c.kegg_id = $kegg_id, c.inchi = $inchi, "
                                "c.chebi_id = $chebi_id, c.lipid_bank_id = $lipid_bank_id,"
                                "c.name = $name,"
                                "c.pubchem_cid = $pubchem_id, c.hmdb_id = $hmdb_id, c.time_stamp = timestamp() "
                                ,
                                lipid_maps_id=lipid,
                                smiles=canonical_smiles,
                                inchikey=lipid_container.getInchiKey(),
                                formula=lipid_container.getFormula(),
                                kegg_id=kegg_id,
                                chebi_id=chebi_id, lipid_bank_id=lipid_bank_id,
                                pubchem_id=pubchem_id, hmdb_id=hmdb_id,
                                inchi=lipid_container.getInchi(), name=name
                                )

            else:
                to_merge_node = res[0]
                if BOIMMGDatabases.LIPID_MAPS.value not in to_merge_node.aliases.keys():
                    with tx.session() as session:
                        session.run("MATCH (c:Compound) "
                                    "where id(c) = $ont_id "
                                    "set c.lipid_maps_id = $lipid_maps_id,"
                                    "c.kegg_id = $kegg_id",
                                    ont_id=to_merge_node.id,
                                    lipid_maps_id=lipid_container.getDbId(),
                                    kegg_id=kegg_id)

        i += 1


def get_filtered_compounds(smiles_keys, smarts, smiles_database):
    new_smiles_list = smiles_keys.copy()
    filtered_compounds = []
    for smile in smiles_keys:
        molecule = MolFromSmiles(smile)
        if molecule:
            match = molecule.GetSubstructMatch(smarts)

            if match:
                filtered_compounds.append(smiles_database[smile])
                new_smiles_list.remove(smile)

    smiles_keys = new_smiles_list.copy()

    return (filtered_compounds, smiles_keys)


def filter_lipid_maps_compounds_by_model_seed_generic_compound(generic_compound: ModelSeedCompound,
                                                               lipid_maps_database) -> list:
    smiles_database = lipid_maps_database.getSmilesDatabase()
    smiles = generic_compound.getSmiles()

    neutralized_smiles, boo = neutralise_charges(smiles)
    smarts = MolFromSmarts(neutralized_smiles)
    smiles_keys = list(smiles_database.keys())

    filtered_compounds, _ = get_filtered_compounds(smiles_keys, smarts, smiles_database)
    return filtered_compounds


def process_and_check_side_chain(sidechains, compound_mol, core_smart,
                                 side_chain_atoms_rules,
                                 side_chain_molecule_rules):
    sidechains.UpdatePropertyCache()
    Chem.GetSymmSSSR(sidechains)

    finalChecker = SidechainChecker()
    finalCheck = finalChecker(core_smart, compound_mol,
                              side_chain_atoms_rules,
                              side_chain_molecule_rules)

    return finalCheck


def establish_relationship_between_generic_and_lipid_maps_compounds(tx, targets, core, side_chain_atoms_rules={},
                                                                    side_chain_molecule_rules={}):
    accessor = CompoundsDBAccessor()
    model_seed_db = ModelSeedCompoundsDB()

    # id_converter = CompoundsIDConverter()

    seen = []

    lipid_maps_database = LipidMapsStructureDB()
    # lipid = lipid_maps_database.get_compound_by_id("LMFA01000001")
    # aliases = lipid.getAliases()
    # targets = []
    #
    # if "KEGG" in aliases.keys():
    #     kegg_id = aliases.get("KEGG")[0]
    #     model_seed_ids = id_converter.convert_db_id_to_model_seed_by_db_id(kegg_id)
    #     targets.append(model_seed_ids[0])

    for target in targets:
        print("Handling %s" % target)

        if target not in seen:
            modelseed_container = model_seed_db.get_compound_by_id(target)

            if "R" in modelseed_container.getFormula():
                # node = accessor.check_if_node_exist_by_model_seed_id(modelseed_id[0])
                # canonical_smiles = MolToSmiles(MolFromSmiles(modelseed_container.getSmiles()))
                parent_ont_container = accessor.get_node_from_model_seed_id(target)

                filtered_compounds = filter_lipid_maps_compounds_by_model_seed_generic_compound(
                    modelseed_container, lipid_maps_database)

                print("establishing reaction with %s" % parent_ont_container.name)
                i = 0

                for compound in filtered_compounds:
                    i += 1
                    printProgressBar(i, len(filtered_compounds))

                    compound_mol = MolFromSmiles(compound.getSmiles())

                    neutralized_smiles, boo = chemo_utilities.neutralise_charges(modelseed_container.getSmiles())
                    smarts = MolFromSmarts(neutralized_smiles)

                    sidechain_compound = AllChem.DeleteSubstructs(compound_mol, smarts)
                    finalCheck = process_and_check_side_chain(sidechain_compound, compound_mol, core,
                                                              side_chain_atoms_rules,
                                                              side_chain_molecule_rules)

                    if finalCheck:
                        sidechain, sidechains_smiles_list = chemo_utilities.retrieve_fragments(sidechain_compound)

                        side_chain_counter = len(sidechain)

                        n_side_chains = modelseed_container.getSmiles().count("*")

                        if n_side_chains == 0 or side_chain_counter == n_side_chains:

                            compound_ont = accessor.get_node_from_lipid_maps_id(compound.getDbId())

                            if compound_ont:
                                successors = accessor.get_all_successors_by_ont_id_rel_type(compound_ont.id, "is_a")

                                if parent_ont_container.id not in successors:
                                    with tx.session() as session:
                                        session.run("MATCH (c:Compound),(d:Compound) "
                                                    "WHERE ID(c)=$target and ID(d) = $origin "
                                                    "MERGE (c)<-[r:is_a]-(d) "
                                                    "ON CREATE SET r.time_stamp = TIMESTAMP() "
                                                    ,
                                                    target=parent_ont_container.id, origin=compound_ont.id)

                                    establish_components_relationships(sidechains_smiles_list, core, tx, compound_ont)

                if not filtered_compounds:
                    print("not found children for %s" % target)

            seen.append(target)


def establish_components_relationships(sidechains_smiles_list, core, tx, compound_ont):
    for sidechain in sidechains_smiles_list:
        sidechain += core
        try:
            mol = MolFromSmiles(sidechain)

            inchikey = MolToInchiKey(mol)

            with tx.session() as session:
                session.run("MATCH (c:Compound),(d:Compound) "
                            "WHERE ID(c)=$target and d.inchikey contains $inchikey "
                            "MERGE (c)<-[r:component_of]-(d) "
                            "ON CREATE SET r.time_stamp = TIMESTAMP() ",
                            inchikey=inchikey[:-1], target=compound_ont.id)

        except:
            print("Warning: side chains not viable")


def get_targets(pathway):
    ptw = extract_metacyc_pathways(pathway)
    ptw.save_in_file(pathway)
    # read_and_filter(pathway)
    targets = []
    with open("generic_targets", "r") as file:
        lines = file.readlines()
        for line in lines:
            target = line.strip()
            if target not in targets:
                targets.append(target)

    return targets


def extract_metacyc_pathways(parent):
    ptw = PathwayHandler()
    biocyc.set_organism('meta')
    ptw_class = biocyc.get(parent)
    l = [ptw_class]
    res = []
    while len(l) > 0:
        try:
            node = l.pop(0)
            if node != ptw_class:

                if not node.subclasses and not node.instances:
                    res.append(node)
                    print(node.name)

            following = node.subclasses
            following.extend(node.instances)
            for elem in following:
                if elem not in res and elem not in l:
                    l.insert(0, elem)

        except:
            pass

    i = 0
    for instance in res:
        try:
            printProgressBar(i, len(res))
            ptw.add_metacyc_pathway(instance.id)
            i += 1

        except:
            pass

    ptw.save_in_file(parent + "_pathway.txt")
    return ptw


####################################################### reactions ontology ##############################################

def reactions_ontology_generation(to_ignore=[]):
    model_seed_reactions_db = ModelSeedReactionsDB()
    model_seed_compounds_db = ModelSeedCompoundsDB()
    reactions_id_converter = ReactionsIDConverter()
    accessor = CompoundsDBAccessor()

    uri = "bolt://palsson.di.uminho.pt:6094"
    driver = GraphDatabase.driver(uri, auth=("neo4j", "Tuning999"), encrypted=False)

    file = open("statistics", "w")
    file.write("reaction\tbad_relations\ttotal\tratio\n")

    file2 = open("bad_relations_products", "w")
    file2.write("reaction\tcompound\n")

    with driver.session() as session:
        result = session.run("match ()<-[r:precursor_of]-() return distinct r.reaction as reaction")

        reactions = result.data()

        for reaction in reactions:

            reaction = reaction.get("reaction")

            if reaction not in to_ignore:

                model_seed_ids = reactions_id_converter.convert_db_id_to_model_seed_by_db_id(reaction)

                model_seed_id = model_seed_ids[0]

                model_seed_reaction = model_seed_reactions_db.getReaction(model_seed_id)

                # model_seed_reaction = model_seed_reactions_db.getReaction("rxn19256")

                products = model_seed_reaction.getProducts()
                reactants = model_seed_reaction.getReactants()
                ms_stoichiometry = model_seed_reaction.getStoichiometry()

                products_in_ontology = []

                final_ms_stoichiometry = {}
                onto_stoichiometry = {}

                for product in products:

                    compound_in_ontology = accessor.get_node_from_model_seed_id(product)

                    if not compound_in_ontology:
                        session.run("merge (c:ModelSeedCompound {model_seed_id: $compound_ms_id})"
                                    , compound_ms_id=product)

                        final_ms_stoichiometry[product] = ms_stoichiometry[product]

                    else:
                        products_in_ontology.append(compound_in_ontology)
                        onto_stoichiometry[compound_in_ontology.id] = ms_stoichiometry[product]

                for reactant in reactants:
                    compound_in_ontology = accessor.get_node_from_model_seed_id(reactant)

                    if not compound_in_ontology:
                        session.run("merge (c:ModelSeedCompound {model_seed_id: $compound_ms_id})"
                                    , compound_ms_id=reactant)

                        final_ms_stoichiometry[reactant] = ms_stoichiometry[reactant]

                    else:
                        products_in_ontology.append(compound_in_ontology)
                        onto_stoichiometry[compound_in_ontology.id] = ms_stoichiometry[reactant]

                add_generic_reaction_and_generate_other(reaction, onto_stoichiometry,
                                                        final_ms_stoichiometry,
                                                        model_seed_reaction,
                                                        model_seed_compounds_db, model_seed_reactions_db,
                                                        file, file2)

                # add_generic_reaction_and_generate_other("CARDIOLIPSYN-RXN", onto_stoichiometry,
                #                                         final_ms_stoichiometry,
                #                                         model_seed_reaction,
                #                                         model_seed_compounds_db, model_seed_reactions_db,
                #                                         file, file2)

                file = open("statistics", "a+")
                file2 = open("bad_relations_products", "a+")

        file.close()
        file2.close()


def add_generic_reaction_and_generate_other(reaction_id, onto_stoichiometry,
                                            ms_stoichiometry,
                                            model_seed_reaction,
                                            model_seed_compound_db,
                                            model_seed_reactions_db,
                                            file, file2):
    reactions_accessor = ReactionsDBAccessor()
    compounds_accessor = CompoundsDBAccessor()

    parent_reaction = reactions_accessor.add_reaction(name=model_seed_reaction.getName(),
                                                      onto_stoichiometry=onto_stoichiometry,
                                                      ms_stoich=ms_stoichiometry,
                                                      direction=model_seed_reaction.getDirection(),
                                                      ec_numbers=model_seed_reaction.getEcNumbers(),
                                                      annotated=True,
                                                      deltaG=model_seed_reaction.getDeltaG(),
                                                      model_seed_id=model_seed_reaction.getDbId())

    original_ms_stoichiometry = ms_stoichiometry.copy()

    print(reaction_id)

    products = []
    reactants = []

    for compound in onto_stoichiometry:

        if onto_stoichiometry[compound] > 0:
            products.append(compound)
        else:
            reactants.append(compound)

    if len(products) == 1:

        generic_product = products[0]

        product_children = compounds_accessor.get_predecessors_by_ont_id_rel_type(generic_product, "is_a")

        j = 1
        i = 0

        bad_relationships = 0
        for product in product_children:
            printProgressBar(i, len(product_children))
            i += 1
            parent_coef = onto_stoichiometry[generic_product]

            new_product_onto_stoichiometry = onto_stoichiometry.copy()

            ms_stoichiometry = original_ms_stoichiometry.copy()

            del new_product_onto_stoichiometry[generic_product]
            new_product_onto_stoichiometry[product] = parent_coef

            precursors = compounds_accessor.get_predecessors_by_ont_id_reaction_id(product, reaction_id)
            precursor_and_correspondent_parent = {}

            go = True
            seen_parents = [generic_product]
            for precursor in precursors:

                parent = compounds_accessor.get_successors_by_ont_id_rel_type(precursor, "is_a")
                if parent:

                    precursor_and_correspondent_parent[precursor] = parent[0]

                    if parent[0] in onto_stoichiometry.keys():

                        coef = onto_stoichiometry[parent[0]]

                        if parent[0] in new_product_onto_stoichiometry:
                            del new_product_onto_stoichiometry[parent[0]]
                            seen_parents.append(parent[0])

                        new_product_onto_stoichiometry[precursor] = coef

                    else:
                        go = False
                        break
                else:
                    go = False
                    break

            if go and sorted(list(onto_stoichiometry.keys())) == sorted(seen_parents):
                equation = get_equation(new_product_onto_stoichiometry, ms_stoichiometry, model_seed_compound_db,
                                        compounds_accessor)

                go2 = True
                if "R" not in equation:
                    new_onto_stoich, new_ms_stoich = \
                        balance_reaction(equation, new_product_onto_stoichiometry, ms_stoichiometry, compounds_accessor,
                                         model_seed_compound_db)

                    if new_onto_stoich and new_ms_stoich:
                        reactions = check_if_reaction_in_model_seed(new_onto_stoich, new_ms_stoich, compounds_accessor,
                                                                    model_seed_reactions_db)

                    else:
                        bad_relationships += 1
                        file2.write(reaction_id + "\t" + str(product) + "\n")
                        go2 = False


                else:
                    new_ms_stoich = ms_stoichiometry.copy()
                    reactions = check_if_reaction_in_model_seed(new_product_onto_stoichiometry, ms_stoichiometry,
                                                                compounds_accessor,
                                                                model_seed_reactions_db)

                    new_onto_stoich = new_product_onto_stoichiometry.copy()

                if go2:
                    if reactions:
                        ms_reaction_id = reactions[0]
                        new_model_seed_reaction = model_seed_reactions_db.getReaction(ms_reaction_id)

                        ont_reaction_id = reactions_accessor.add_reaction(name=new_model_seed_reaction.getName(),
                                                                          onto_stoichiometry=new_onto_stoich,
                                                                          ms_stoich=new_ms_stoich,
                                                                          direction=new_model_seed_reaction.getDirection(),
                                                                          ec_numbers=new_model_seed_reaction.getEcNumbers(),
                                                                          annotated=True,
                                                                          deltaG=new_model_seed_reaction.getDeltaG(),
                                                                          model_seed_id=new_model_seed_reaction.getDbId())

                    else:

                        ont_reaction_id = reactions_accessor.add_reaction(
                            name=model_seed_reaction.getName() + "_granulated_" + str(j),
                            onto_stoichiometry=new_onto_stoich,
                            ms_stoich=new_ms_stoich,
                            direction=model_seed_reaction.getDirection(),
                            ec_numbers=model_seed_reaction.getEcNumbers(),
                            annotated=False,
                            deltaG=None)

                        j += 1

                    reactions_accessor.add_relationship(ont_reaction_id, parent_reaction, "is_a")


                else:
                    bad_relationships += 1
                    file2.write(reaction_id + "\t" + str(product))

        if len(product_children) != 0:
            file.write(reaction_id + "\t" + str(bad_relationships) + "\t" + str(len(product_children)) + "\t" +
                       str((bad_relationships / len(product_children))) + "\n")

            print("bad_relationships: %d" % bad_relationships)
            print("ratio: %d" % (bad_relationships / len(product_children)))

        file.close()
        file2.close()


def check_if_reaction_in_model_seed(new_onto_stoich, new_ms_stoich, accessor: CompoundsDBAccessor,
                                    reactions_db: ModelSeedReactionsDB):
    reactants = []
    products = []
    for onto_compound in new_onto_stoich:
        node = accessor.get_node_by_ont_id(onto_compound)
        aliases = node.aliases

        if "ModelSEED" in aliases:
            alias = aliases["ModelSEED"][0]

            if new_onto_stoich[onto_compound] < 0:
                reactants.append(alias)

            else:
                products.append(alias)

        else:
            return None

    for ms_compound in new_ms_stoich:

        if new_ms_stoich[ms_compound] < 0:
            reactants.append(ms_compound)

        else:
            products.append(ms_compound)

    reactions = reactions_db.getReactionsByCompounds(reactants, products)
    return reactions


def balance_reaction(eq, onto_stoichiometry, ms_stoichiometry, accessor, ms_db):
    k, Ys = lp_to_balance_reaction(eq)

    hydrogen = "cpd00067"

    isHydrogenInReaction = False
    if hydrogen in ms_stoichiometry:
        isHydrogenInReaction = True

    if len(k) == len(Ys):

        try:
            new_onto_stoich, new_ms_stoich = \
                process_result_of_new_balanced_reaction(k, Ys, eq, onto_stoichiometry, ms_stoichiometry, accessor,
                                                        ms_db)
            return new_onto_stoich, new_ms_stoich

        except:
            return None, None
        # print('->'.join('+'.join(pM(N.pop(0)) + str(t) for t in p.split('+')) for p in eq.split('->')))

    elif not isHydrogenInReaction:

        reactants_products = eq.split("->")
        reactants_products[0] += "+H"

        new_eq = "->".join(reactants_products)
        k, Ys = lp_to_balance_reaction(new_eq)

        if k and len(k) == len(Ys):

            ms_stoichiometry[hydrogen] = -1

            try:
                new_onto_stoich, new_ms_stoich = \
                    process_result_of_new_balanced_reaction(k, Ys, new_eq, onto_stoichiometry, ms_stoichiometry,
                                                            accessor,
                                                            ms_db)
                return new_onto_stoich, new_ms_stoich

            except:
                return None, None
        else:

            reactants_products = eq.split("->")
            reactants_products[1] += "+H"

            new_eq = "->".join(reactants_products)
            k, Ys = lp_to_balance_reaction(new_eq)

            if k and len(k) == len(Ys):

                ms_stoichiometry[hydrogen] = 1

                try:
                    new_onto_stoich, new_ms_stoich = \
                        process_result_of_new_balanced_reaction(k, Ys, new_eq, onto_stoichiometry, ms_stoichiometry,
                                                                accessor,
                                                                ms_db)
                    return new_onto_stoich, new_ms_stoich

                except:
                    return None, None

    else:

        reactants_products = eq.split("->")

        reactants = reactants_products[0].split("+")
        products = reactants_products[1].split("+")

        if "H" in reactants:
            reactants.remove("H")

            new_reactants = "+".join(reactants)
            new_eq = "->".join([new_reactants, reactants_products[1]])

            k, Ys = lp_to_balance_reaction(new_eq)

            if k and len(k) == len(Ys):

                ms_stoichiometry[hydrogen] = -1
                try:
                    new_onto_stoich, new_ms_stoich = \
                        process_result_of_new_balanced_reaction(k, Ys, new_eq, onto_stoichiometry, ms_stoichiometry,
                                                                accessor,
                                                                ms_db)
                    return new_onto_stoich, new_ms_stoich

                except:

                    return None, None

        else:
            products.remove("H")

            new_reactants = "+".join(products)
            new_eq = "->".join([new_reactants, reactants_products[1]])

            k, Ys = lp_to_balance_reaction(new_eq)

            if k and len(k) == len(Ys):

                ms_stoichiometry[hydrogen] = 1

                try:
                    new_onto_stoich, new_ms_stoich = \
                        process_result_of_new_balanced_reaction(k, Ys, new_eq, onto_stoichiometry, ms_stoichiometry,
                                                                accessor,
                                                                ms_db)
                    return new_onto_stoich, new_ms_stoich
                except:
                    return None, None

    return None, None


def lp_to_balance_reaction(eq):
    Ls = list('abcdefghijklmnopqrstuvwxyz')

    Ss, Os, Es, a, i = defaultdict(list), Ls[:], [], 1, 1

    for p in eq.split('->'):
        for k in p.split('+'):
            c = [Ls.pop(0), 1]
            for e, m in re.findall('([A-Z][a-z]?)([0-9]*)', k):
                m = 1 if m == '' else int(m)
                a *= m
                d = [c[0], c[1] * m * i]
                Ss[e][:0], Es[:0] = [d], [[e, d]]
        i = -1

    Ys = dict((s, eval('Symbol("' + s + '")')) for s in Os if s not in Ls)

    Qs = [eval('+'.join('%d*%s' % (c[1], c[0]) for c in Ss[s]), {}, Ys) for s in Ss] + [Ys['a'] - a]

    k = solve(Qs, *Ys)
    return k, Ys


def process_result_of_new_balanced_reaction(k, Ys, eq, onto_stoichiometry, ms_stoich, accessor,
                                            ms_db: ModelSeedCompoundsDB):
    N = []
    for s in sorted(Ys):
        value = k[Ys[s]]
        N.append(value)

    g = N[0]

    for a1, a2 in zip(N[0::2], N[1::2]):
        g = math.gcd(g, a2)

    N = [i / g for i in N]

    pM = lambda c: str(c) if c != 1 else ''

    res_onto_stoichiometry = {}
    res_ms_stoich = {}

    equation_reactants_products = eq.split('->')
    reactants_in_equation = equation_reactants_products[0]
    products_in_equation = equation_reactants_products[1]

    for reactant in reactants_in_equation.split('+'):

        for reactant2 in onto_stoichiometry:
            if onto_stoichiometry[reactant2] < 0:

                node = accessor.get_node_by_ont_id(reactant2)
                if node.formula == reactant:
                    stoich = pM(N.pop(0))
                    if stoich:
                        coef = abs(int(stoich))
                    else:
                        coef = 1

                    if coef != 0:
                        res_onto_stoichiometry[reactant2] = coef * -1
                        break

        for reactant3 in ms_stoich:
            if ms_stoich[reactant3] < 0:

                ms_compound = ms_db.get_compound_by_id(reactant3)
                if ms_compound.getFormula() == reactant:
                    stoich = pM(N.pop(0))
                    if stoich:
                        coef = abs(int(stoich))
                    else:
                        coef = 1

                    if coef != 0:
                        res_ms_stoich[reactant3] = coef * -1
                        break

    for product in products_in_equation.split('+'):

        for product2 in onto_stoichiometry:

            if onto_stoichiometry[product2] > 0:
                node = accessor.get_node_by_ont_id(product2)

                if node.formula == product:
                    stoich = pM(N.pop(0))
                    if stoich:
                        coef = abs(int(stoich))
                    else:
                        coef = 1

                    if coef != 0:
                        res_onto_stoichiometry[product2] = coef
                        break

        for product3 in ms_stoich:
            if ms_stoich[product3] > 0:

                ms_compound = ms_db.get_compound_by_id(product3)
                if ms_compound.getFormula() == product:
                    stoich = pM(N.pop(0))
                    if stoich:
                        coef = abs(int(stoich))
                    else:
                        coef = 1

                    if coef != 0:
                        res_ms_stoich[product3] = coef
                        break

    return res_onto_stoichiometry, res_ms_stoich


def get_equation(onto_stoichiometry: dict, ms_stoichiometry: dict,
                 ms_comp_database: ModelSeedCompoundsDB, accessor: CompoundsDBAccessor):
    reactants_formulas = []
    products_formulas = []
    for compound in onto_stoichiometry:
        node = accessor.get_node_by_ont_id(compound)

        if onto_stoichiometry[compound] < 0:
            reactants_formulas.append(node.formula)
        else:
            products_formulas.append(node.formula)

    for compound in ms_stoichiometry:
        ms_compound = ms_comp_database.get_compound_by_id(compound)

        if ms_stoichiometry[compound] < 0:
            reactants_formulas.append(ms_compound.getFormula())
        else:
            products_formulas.append(ms_compound.getFormula())

    eq = ""

    i = 0
    while i < len(reactants_formulas) - 1:
        reactant_formula = reactants_formulas[i]
        eq += reactant_formula + "+"

        i += 1
    eq += reactants_formulas[i] + "->"

    j = 0
    while j < len(products_formulas) - 1:
        product_formula = products_formulas[j]
        eq += product_formula + "+"

        j += 1

    eq += products_formulas[j]
    return eq


############################################################# pipeline #######################################################
def pipeline():


    read_swiss_lipids_database()

    scrap_lipid_maps_db(driver)

    core = "C(=O)O"

    targets = get_targets("Phospholipid-Biosynthesis")

    establish_relationship_between_generic_and_lipid_maps_compounds(driver, targets=targets, core=core)
    #
    # generator = OntologyGeneratorDA()
    #
    # #### CARDIOLIPIN BIOSYNTHESIS TRANSFORMATION AND RELATIONSHIPS #####################
    # transformations = TransformationsHandler()
    # transformations.choose_transformations_from_pathway("Cardiolipin-Biosynthesis")
    # transformations.save_transformations("transformations_cardiolipin")
    #
    # generator("Cardiolipin-Biosynthesis", core, "transformations_cardiolipin")


def integrate_all_model_seed_compounds(driver):
    ms_compounds_db = ModelSeedCompoundsDB()

    accessor = CompoundsDBAccessor()

    j = 0
    i = 0
    for smiles in ms_compounds_db.isomer_smiles_database:

        printProgressBar(i, len(ms_compounds_db.isomer_smiles_database))
        i += 1
        if "*" in smiles:
            if len(ms_compounds_db.isomer_smiles_database[smiles]) > 1:

                isomer_list = [iso.getDbId() for iso in ms_compounds_db.isomer_smiles_database[smiles]]

                chosen_boimmg_id = None
                for iso in isomer_list:
                    node = accessor.get_node_from_model_seed_id(iso)

                    if node:
                        chosen_boimmg_id = node.id
                        break

                if chosen_boimmg_id:
                    for iso in isomer_list:
                        j += 1
                        with driver.session() as session:
                            session.run("match (c:Compound) where id(c) = $node_id "
                                        "merge (c)<-[:is_db_link_of]-(:ModelSeedCompound {model_seed_id: $iso})",
                                        node_id=chosen_boimmg_id, iso=iso)
    print(j)


if __name__ == "__main__":
    pipeline()

    # targets = get_targets("Sphingolipid-Biosynthesis")
    # establish_relationship_between_generic_and_lipid_maps_compounds(driver, targets=targets, core=core)

    # reactions_db = ModelSeedReactionsDB()
    # reaction = reactions_db.getReaction("rxn23713")
    # add_generic_reaction_and_generate_other()

    reactions_ontology_generation(
        ["RXN-10462", "RXN-1623", "RXN-9591", "1-ACYLGLYCEROL-3-P-ACYLTRANSFER-RXN", "RXN-9590",
         "CDPDIGLYSYN-RXN", "RXN0-7012", "RXN-8141", "CARDIOLIPSYN-RXN", "CARDIOLIPSYN-RXN", "PHOSPHASERDECARB-RXN",
         "ETHANOLAMINEPHOSPHOTRANSFERASE-RXN", "PHOSPHAGLYPSYN-RXN", "PGPPHOSPHA-RXN", "2.3.1.23-RXN",
         "RXN-15380", "RXN4FS-1", "2.1.1.17-RXN", "2.1.1.71-RXN", "RXN4FS-4", "RXN4FS-2",
         "RXN-5781", "2.7.8.24-RXN"])

    # integrate_all_model_seed_compounds(driver)
