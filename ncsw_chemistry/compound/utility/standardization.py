""" The ``ncsw_chemistry.compound.utility`` package ``standardization`` module. """

from functools import reduce
from typing import Collection, Dict, Optional

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeFlags, SanitizeMol


class CompoundStandardizationUtility:
    """ The chemical compound standardization utility class. """

    @staticmethod
    def get_compound_sanitization_operations() -> Dict[str, SanitizeFlags]:
        """
        Get the chemical compound sanitization operations.

        :returns: The chemical compound sanitization operations.
        """

        return {
            "adjust_hydrogens": SanitizeFlags.SANITIZE_ADJUSTHS,
            "all": SanitizeFlags.SANITIZE_ALL,
            "clean_up": SanitizeFlags.SANITIZE_CLEANUP,
            "clean_up_atropisomers": SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS,
            "clean_up_chirality": SanitizeFlags.SANITIZE_CLEANUPCHIRALITY,
            "clean_up_organometallics": SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS,
            "find_radicals": SanitizeFlags.SANITIZE_FINDRADICALS,
            "kekulize": SanitizeFlags.SANITIZE_KEKULIZE,
            "none": SanitizeFlags.SANITIZE_NONE,
            "properties": SanitizeFlags.SANITIZE_PROPERTIES,
            "set_aromaticity": SanitizeFlags.SANITIZE_SETAROMATICITY,
            "set_conjugation": SanitizeFlags.SANITIZE_SETCONJUGATION,
            "set_hybridization": SanitizeFlags.SANITIZE_SETHYBRIDIZATION,
            "symmetrize_rings": SanitizeFlags.SANITIZE_SYMMRINGS,
        }

    @staticmethod
    def sanitize_compound(
            compound_mol: Mol,
            compound_sanitization_operation_keys: Collection[str] = None,
            deep_copy: bool = True
    ) -> Optional[Mol]:
        """
        Sanitize a chemical compound.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter compound_sanitization_operation_keys: The keys of the chemical compound sanitization operations that
            should be performed. The value `None` indicates that all chemical compound sanitization operations should be
            performed.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound RDKit Mol object should be
            constructed and modified.

        :returns: The sanitized chemical compound.
        """

        if deep_copy:
            compound_mol = Mol(compound_mol)

        if compound_sanitization_operation_keys is None:
            SanitizeMol(
                mol=compound_mol
            )

        elif len(compound_sanitization_operation_keys) == 0:
            SanitizeMol(
                mol=compound_mol,
                sanitizeOps=SanitizeFlags.SANITIZE_NONE
            )

        else:
            compound_sanitization_operations = CompoundStandardizationUtility.get_compound_sanitization_operations()

            SanitizeMol(
                mol=compound_mol,
                sanitizeOps=reduce(
                    lambda compound_sanitization_operation_value_a, compound_sanitization_operation_value_b: (
                        compound_sanitization_operation_value_a | compound_sanitization_operation_value_b
                    ), [
                        compound_sanitization_operations[compound_sanitization_operation_key]
                        for compound_sanitization_operation_key in compound_sanitization_operation_keys
                        if compound_sanitization_operation_key in compound_sanitization_operations.keys()
                    ]
                )
            )

        return compound_mol
