""" The ``ncsw_chemistry.utility.compound`` package ``sanitization`` module. """

from copy import deepcopy

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeMol


class CompoundSanitizationUtility:
    """ The chemical compound sanitization utility class. """

    @staticmethod
    def sanitize(
            compound_mol: Mol,
            deep_copy: bool = True
    ) -> Mol:
        """
        Sanitize a chemical compound.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit Mol` object should be
            constructed and modified.

        :returns: The sanitized chemical compound.
        """

        if deep_copy:
            compound_mol = deepcopy(
                x=compound_mol
            )

        SanitizeMol(compound_mol)

        return compound_mol
