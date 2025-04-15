""" The ``ncsw_chemistry.reaction.utility`` package ``standardization`` module. """

from functools import reduce
from typing import Collection, Dict

from rdkit.Chem.rdChemReactions import ChemicalReaction, SanitizeFlags, SanitizeRxn

from ncsw_chemistry.compound.utility.standardization import CompoundStandardizationUtility


class ReactionStandardizationUtility:
    """ The chemical reaction standardization utility class. """

    @staticmethod
    def get_reaction_sanitization_operations() -> Dict[str, SanitizeFlags]:
        """
        Get the chemical reaction sanitization operations.

        :returns: The chemical reaction sanitization operations.
        """

        return {
            "adjust_reactants": SanitizeFlags.SANITIZE_ADJUST_REACTANTS,
            "all": SanitizeFlags.SANITIZE_ALL,
            "atom_map_numbers": SanitizeFlags.SANITIZE_ATOM_MAPS,
            "merge_hydrogens": SanitizeFlags.SANITIZE_MERGEHS,
            "none": SanitizeFlags.SANITIZE_NONE,
            "r_group_names": SanitizeFlags.SANITIZE_RGROUP_NAMES,
        }

    @staticmethod
    def sanitize_reaction(
            reaction_rxn: ChemicalReaction,
            reaction_sanitization_operation_keys: Collection[str] = None,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Sanitize a chemical reaction.

        :parameter reaction_rxn: The RDKit ChemicalReaction object of the chemical reaction.
        :parameter reaction_sanitization_operation_keys: The keys of the chemical reaction sanitization operations that
            should be performed. The value `None` indicates that all chemical compound sanitization operations should be
            performed.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical reaction RDKit ChemicalReaction
            object should be constructed and modified.

        :returns: The sanitized chemical reaction.
        """

        if deep_copy:
            reaction_rxn = ChemicalReaction(reaction_rxn)

        if reaction_sanitization_operation_keys is None:
            SanitizeRxn(
                rxn=reaction_rxn
            )

        elif len(reaction_sanitization_operation_keys) == 0:
            SanitizeRxn(
                rxn=reaction_rxn,
                sanitizeOps=SanitizeFlags.SANITIZE_NONE
            )

        else:
            reaction_sanitization_operations = ReactionStandardizationUtility.get_reaction_sanitization_operations()

            SanitizeRxn(
                rxn=reaction_rxn,
                sanitizeOps=reduce(
                    lambda reaction_sanitization_operation_value_a, reaction_sanitization_operation_value_b: (
                        reaction_sanitization_operation_value_a | reaction_sanitization_operation_value_b
                    ), [
                        reaction_sanitization_operations[reaction_sanitization_operation_key]
                        for reaction_sanitization_operation_key in reaction_sanitization_operation_keys
                        if reaction_sanitization_operation_key in reaction_sanitization_operations.keys()
                    ]
                )
            )

        return reaction_rxn

    @staticmethod
    def sanitize_reaction_compounds(
            reaction_rxn: ChemicalReaction,
            reaction_compound_sanitization_operation_keys: Collection[str] = None,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Sanitize the compounds of a chemical reaction.

        :parameter reaction_rxn: The RDKit ChemicalReaction object of the chemical reaction.
        :parameter reaction_compound_sanitization_operation_keys: The keys of the chemical reaction compound
            sanitization operations that should be performed. The value `None` indicates that all chemical reaction
            compound sanitization operations should be performed.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical reaction RDKit ChemicalReaction
            object should be constructed and modified.

        :returns: The chemical reaction with sanitized compounds.
        """

        if deep_copy:
            reaction_rxn = ChemicalReaction(reaction_rxn)

        for reaction_compound_mols in (
            reaction_rxn.GetReactants(),
            reaction_rxn.GetAgents(),
            reaction_rxn.GetProducts(),
        ):
            for reaction_compound_mol in reaction_compound_mols:
                CompoundStandardizationUtility.sanitize_compound(
                    compound_mol=reaction_compound_mol,
                    compound_sanitization_operation_keys=reaction_compound_sanitization_operation_keys
                )

        return reaction_rxn
