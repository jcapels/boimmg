import copy
import itertools
import re

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MolToSmiles, MolFromSmiles, MolFromSmarts


def check_n_to_m_atoms(molecule,n,m,atom_number):
    atoms = molecule.GetAtoms()
    n_atoms = 0
    for atom in atoms:
        if atom.GetAtomicNum() == atom_number:
            n_atoms+=1

    if n_atoms>=n and n_atoms<=m:
        return True

    else:
        return False

def check_at_least_n_atoms(molecule,n,atom_number):
    atoms = molecule.GetAtoms()
    n_atoms = 0
    for atom in atoms:
        if atom.GetAtomicNum() == atom_number:
            n_atoms += 1

    if n_atoms >= n:
        return True

    else:
        return False

def check_if_odd_number_of_atom(molecule,atom_number,backbone_n_atoms=0):
    atoms = molecule.GetAtoms()
    n_atoms = backbone_n_atoms
    for atom in atoms:
        if atom.GetAtomicNum() == atom_number:
            n_atoms += 1

    if n_atoms%2 != 0:
        # print(MolToSmiles(molecule))
        return True

    else:
        return False

def get_side_chains(parent_smiles,child_smiles):

    r_groups_number = parent_smiles.count("*")

    neutralized_smiles, _ = NeutraliseCharges(parent_smiles)

    transformation1 = convert_model_seed_dummy_atoms(parent_smiles) + ">>"
    component_smiles = "*C(=O)(O)"

    component_smiles_list = []
    for i in range(1, r_groups_number + 1):
        component_smiles_list.append(component_smiles.replace("*", "[*:" + str(i) + "]"))

    transformation2 = convert_model_seed_dummy_atoms(neutralized_smiles) + ">>"

    transformation1 += ".".join(component_smiles_list)
    transformation2 += ".".join(component_smiles_list)

    transformation_smarts1 = AllChem.ReactionFromSmarts(transformation1)

    transformation_smarts2 = AllChem.ReactionFromSmarts(transformation2)

    child_mol = MolFromSmiles(child_smiles)

    try:
        new_molecules = [Chem.MolToSmiles(x, 1) for x in transformation_smarts1.RunReactants(
            tuple([child_mol]))[0]]

        return new_molecules

    except:
        try:
            new_molecules = [Chem.MolToSmiles(x, 1) for x in transformation_smarts2.RunReactants(
                tuple([child_mol]))[0]]

            return new_molecules

        except:
            return None

    return None




def check_if_only_cis_stereochemistry(molecule):

    smiles = MolToSmiles(molecule)

    if "=" in smiles:
        regex_cis = re.compile(r"/[A-Z]{1}[a-z]{0,1}=[A-Z]{1}[a-z]{0,1}\\")
        regex_trans = re.compile(r"/[A-Z]{1}[a-z]{0,1}=[A-Z]{1}[a-z]{0,1}/|\\[A-Z]{1}[a-z]{0,1}=[A-Z]{1}[a-z]{0,1}\\")
        cis_stereochemistry = regex_cis.findall(smiles)
        trans_stereochemistry = regex_trans.findall(smiles)
        double_bonds = smiles.count("=")

    else:
        return True

    if cis_stereochemistry and not trans_stereochemistry and len(cis_stereochemistry) == double_bonds:
        return True
    else:
        return  False
    

def get_atoms_idx_from_molecule(molecule):
    res = []
    atoms = molecule.GetAtoms()
    for atom in atoms:
        res.append(atom.GetIdx())
    return res

def checkIfCycle(molecule):
    info = molecule.GetRingInfo()
    numRings = info.NumRings()
    if numRings == 0:
        return True
    else:
        return False

def __check_if_sidechains_are_equal(lst1,lst2):
    lst_1 = copy.deepcopy(lst1)
    lst_2 = copy.deepcopy(lst2)
    for mol1 in lst_1:
        i = 0
        found = False
        while i < len(lst_2):
            fgp1 = AllChem.GetMorganFingerprint(mol1,1)
            fgp2 = AllChem.GetMorganFingerprint(lst2[i],1)
            similarity = DataStructs.TanimotoSimilarity(fgp1, fgp2)
            if similarity == 1:
                lst_2.remove(lst_2[i])
                found = True
            else:
                i += 1

        if not found:
            return False

    if not lst_2:
        return True

    else:
        return False

def check_if_sidechains_are_equal(lst1,lst2):

    if len(lst1) <= len(lst2):
        return __check_if_sidechains_are_equal(lst2,lst1)

    else:
        return __check_if_sidechains_are_equal(lst1, lst2)

def check_if_molecules_are_equal(mol1,mol2):
    fg1 = AllChem.GetMorganFingerprint(mol1,1,useFeatures=True,useChirality=True)
    fg2 = AllChem.GetMorganFingerprint(mol2,1,useFeatures=True,useChirality=True)
    similarity = DataStructs.TanimotoSimilarity(fg1, fg2)
    if similarity==1:
        return True
    else:
        return False

def retrieve_fragments(mol):
    substructures_smiles = MolToSmiles(mol)
    sidechains = []
    sidechains_smiles_list = substructures_smiles.split(".")
    for sidechain_smiles in sidechains_smiles_list:
        sidechain = MolFromSmiles(sidechain_smiles)
        if sidechain_smiles != "":
            sidechains.append(sidechain)
    return sidechains,sidechains_smiles_list

def convert_model_seed_dummy_atoms(modelseed_smiles):
    i=1
    new_smart=""
    match1 = re.finditer("\*",modelseed_smiles)

    previous = None
    for m in match1:
        start = m.start()
        end = m.end()
        if not previous:
            new_smart += modelseed_smiles[:start] + "[*:" + str(i) + "]"
        else:
            new_smart += modelseed_smiles[previous:start] + "[*:" + str(i) + "]"

        i+=1
        previous = end

    if previous!=len(modelseed_smiles)-1:
        new_smart+=modelseed_smiles[previous:]
    return new_smart

def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None

def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return (Chem.MolToSmiles(mol,True), True)
    else:
        return (smiles, False)

def check_similarity_between_generic_and_complete_representation(generic_smiles,complete_smiles):
    complete_smiles,_ = NeutraliseCharges(complete_smiles)
    generic_smiles,_ = NeutraliseCharges(generic_smiles)
    complete_mol = MolFromSmiles(complete_smiles)
    generic_mol = MolFromSmarts(generic_smiles)
    match = complete_mol.GetSubstructMatch(generic_mol)
    if match:
        return True
    return False


def __check_if_dummy_is_undefined(start,strange_metabolite):
    for sm in strange_metabolite:
        if start == sm + 1:
            return True
    return False


def get_combinations_of_dummy_atoms(modelseed_smiles):

    r_groups = modelseed_smiles.count("*")
    to_subtract = modelseed_smiles.count("!*")
    r_groups = r_groups-to_subtract

    lst = [i for i in range(1,r_groups+1)]
    combinations = itertools.permutations(lst)
    combinations = list(combinations)

    res = []
    for comb in combinations:
        lst= list(comb)

        new_smart = ""
        match1 = re.finditer("\*", modelseed_smiles)
        strange_metabolite = re.finditer("!", modelseed_smiles)
        strange_metabolite_positions = [sm.start() for sm in strange_metabolite]

        i=0
        previous = None
        for m in match1:
            start = m.start()
            end = m.end()
            undefined_met = False
            if strange_metabolite:
                undefined_met =  __check_if_dummy_is_undefined(start,strange_metabolite_positions)

            if not undefined_met:
                if not previous:
                    new_smart += modelseed_smiles[:start] + "[*:" + str(lst[i]) + "]"
                else:
                    new_smart += modelseed_smiles[previous:start] + "[*:" + str(lst[i]) + "]"

                i += 1
            else:
                new_smart += modelseed_smiles[previous:start+1]
                new_smart = new_smart.replace("!","")


            previous = end



        if previous != len(modelseed_smiles) - 1:
            new_smart += modelseed_smiles[previous:]

        res.append(new_smart)

    return res

def check_if_molecule_has_double_bonds(molecule):
    atoms = molecule.GetAtoms()
    for atom in atoms:
        bonds = atom.GetBonds()
        for bond in bonds:
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                return True
    return False

def check_if_atom_has_double_bonds(atom):
    bonds = atom.GetBonds()
    for bond in bonds:
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            return True
    return False

def check_if_atom_has_triple_bonds(atom):
    bonds = atom.GetBonds()
    for bond in bonds:
        if bond.GetBondType() == Chem.BondType.TRIPLE:
            return True
    return False

def check_if_two_molecules_are_equal_from_smiles(smiles1,smiles2):
    mol1 = MolFromSmiles(smiles1)
    mol2 = MolFromSmiles(smiles2)

    fgp1 = AllChem.GetMorganFingerprint(mol1, 1,useFeatures =True,useChirality=True)
    fgp2 = AllChem.GetMorganFingerprint(mol2, 1,useFeatures =True,useChirality=True)
    similarity = DataStructs.TanimotoSimilarity(fgp1, fgp2)

    if similarity == 1:
        return True

    else:
        return False



