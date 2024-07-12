""" The ``ncsw_chemistry.utility.reaction`` package ``standardization`` module. """

from copy import deepcopy
from functools import reduce
from typing import Iterable

from rdkit.Chem.rdChemReactions import ChemicalReaction, SanitizeFlags, SanitizeRxn


class ReactionStandardizationUtility:
    """ The chemical reaction standardization utility class. """

    @staticmethod
    def sanitize(
            reaction_rxn: ChemicalReaction,
            sanitization_operation_keys: Iterable[str] = None,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Sanitize a chemical reaction.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter sanitization_operation_keys: The keys of the chemical reaction sanitization operations that should be
            performed. The value `None` indicates that all chemical reaction sanitization operations should be
            performed.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit ChemicalReaction`
            object should be constructed and modified.

        :returns: The sanitized chemical reaction.
        """

        if deep_copy:
            reaction_rxn = deepcopy(
                x=reaction_rxn
            )

        sanitization_operations = {
            "adjust_reactants": SanitizeFlags.SANITIZE_ADJUST_REACTANTS,
            "all": SanitizeFlags.SANITIZE_ALL,
            "atom_map_numbers": SanitizeFlags.SANITIZE_ATOM_MAPS,
            "merge_hydrogens": SanitizeFlags.SANITIZE_MERGEHS,
            "none": SanitizeFlags.SANITIZE_NONE,
            "r_group_names": SanitizeFlags.SANITIZE_RGROUP_NAMES,
        }

        sanitization_operation_keys = [
            sanitization_operation_key
            for sanitization_operation_key in sanitization_operation_keys
            if sanitization_operation_key in sanitization_operations.keys()
        ] if sanitization_operation_keys is not None else [
            "all",
        ]

        if len(sanitization_operation_keys) == 0 or "all" in sanitization_operation_keys:
            SanitizeRxn(
                rxn=reaction_rxn,
                sanitizeOps=SanitizeFlags.SANITIZE_ALL
            )

        elif "none" in sanitization_operation_keys:
            SanitizeRxn(
                rxn=reaction_rxn,
                sanitizeOps=SanitizeFlags.SANITIZE_NONE
            )

        else:
            SanitizeRxn(
                rxn=reaction_rxn,
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

        return reaction_rxn
