from rdkit.Chem import AllChem


def matchers_functions(f):
    def wraps(*args, **kwargs):
        return f(*args, **kwargs)

    return wraps


class SidechainChecker:

    def __call__(self, core_smart, mol, side_chain_atoms_rules={}, molecule_atoms_rules={}):

        atoms_checked = True
        molecule_checked = True

        # core_smart = MolFromSmarts(core_smart)
        sidechain_compound = AllChem.DeleteSubstructs(mol, core_smart)

        if side_chain_atoms_rules:
            idx_lst = mol.GetSubstructMatches(core_smart)

            if not idx_lst:
                return False

            tam = mol.GetNumAtoms()
            substructure_idx = []
            for idx in idx_lst:
                for i in range(0, tam):
                    if i not in idx:
                        substructure_idx.append(i)

            atoms_checked = self.__check_sidechain_atoms__(substructure_idx, mol, side_chain_atoms_rules)

        if molecule_atoms_rules:
            molecule_checked = self.__check_molecule_atoms__(molecule_atoms_rules, sidechain_compound)

        if molecule_checked and atoms_checked:
            return True
        else:
            return False

    def __check_sidechain_atoms__(self, vect, mol, matchers_atoms):
        # loop over the atoms we care about:
        boolean_list = []
        for idx in vect:
            atom = mol.GetAtomWithIdx(idx)
            # print(atom.GetAtomicNum())
            for matcher in matchers_atoms:
                if callable(matchers_atoms[matcher]):
                    check = matchers_functions(matchers_atoms[matcher])(atom)

                    boolean_list.append(check)

        if all(boolean_list):
            return True

        else:
            return False

    def __check_molecule_atoms__(self, matchers_molecules, sidechain):
        boolean_list = []
        for matcher in matchers_molecules:
            if callable(matchers_molecules[matcher]):
                check = matchers_functions(matchers_molecules[matcher])(sidechain)
                boolean_list.append(check)

        if all(boolean_list):
            return True
        else:
            return False
