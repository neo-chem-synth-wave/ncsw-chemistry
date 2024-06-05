""" The ``ncsw_chemistry.compound.bond`` package ``property`` module. """

from typing import Dict, Union

from rdkit.Chem.rdchem import Bond


class CompoundBondPropertyUtilities:
    """ The chemical compound bond property utilities class. """

    @staticmethod
    def get_properties(
            bond: Bond
    ) -> Dict[str, Union[bool, int, str]]:
        """
        Get the chemical compound bond properties.

        :parameter bond: The chemical compound bond `RDKit Bond` object.

        :returns: The chemical compound bond properties.
        """

        return {
            "index": bond.GetIndex(),
            "bond_type": str(bond.GetBondType()),
            "is_conjugated": bond.GetIsConjugated(),
            "is_in_ring": bond.IsInRing(),
            "is_aromatic": bond.GetIsAromatic(),
        }
