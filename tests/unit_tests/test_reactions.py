from unittest import TestCase

from rdkit import Chem
from rdkit.Chem import MolToSmiles
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


class TestReactionRule(TestCase):

    def test_reaction_rule(self):
        rule = ReactionFromSmarts('[*:1]C(=O)OC[C@H](COP(=O)(O)OCCN)OC([*:2])=O>>'
                                  '[*:1]C(=O)OC[C@H](COP(=O)(O)OC[C@@H](N)C(=O)O)OC([*:2])=O')
        reacts = (
            Chem.MolFromSmiles('CCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(CCCCCCCC)=O'),
        )

        products = rule.RunReactants(reacts)

        for product_tuple in products:
            print(MolToSmiles(product_tuple[0]))
