""" The ``ncsw_chemistry.utility.compound`` package ``standardization`` module. """

from copy import deepcopy
from functools import reduce
from typing import Iterable

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeFlags, SanitizeMol


class CompoundStandardizationUtility:
    """ The chemical compound standardization utility class. """

    @staticmethod
    def sanitize(
            compound_mol: Mol,
            sanitization_operation_keys: Iterable[str] = None,
            deep_copy: bool = True
    ) -> Mol:
        """
        Sanitize a chemical compound.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter sanitization_operation_keys: The keys of the chemical compound sanitization operations that should be
            performed. The value `None` indicates that all chemical compound sanitization operations should be
            performed.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit Mol` object should be
            constructed and modified.

        :returns: The sanitized chemical compound.
        """

        if deep_copy:
            compound_mol = deepcopy(
                x=compound_mol
            )

        sanitization_operations = {
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

        sanitization_operation_keys = [
            sanitization_operation_key
            for sanitization_operation_key in sanitization_operation_keys
            if sanitization_operation_key in sanitization_operations.keys()
        ] if sanitization_operation_keys is not None else [
            "all",
        ]

        if len(sanitization_operation_keys) == 0 or "all" in sanitization_operation_keys:
            SanitizeMol(
                mol=compound_mol,
                sanitizeOps=SanitizeFlags.SANITIZE_ALL
            )

        elif "none" in sanitization_operation_keys:
            SanitizeMol(
                mol=compound_mol,
                sanitizeOps=SanitizeFlags.SANITIZE_NONE
            )

        else:
            SanitizeMol(
                mol=compound_mol,
                sanitizeOps=reduce(
                    function=(
                        lambda sanitization_operation_a, sanitization_operation_b:
                        sanitization_operation_a | sanitization_operation_b
                    ),
                    sequence=[
                        sanitization_operations[sanitization_operation_key]
                        for sanitization_operation_key in sanitization_operation_keys
                    ]
                )
            )

        return compound_mol
