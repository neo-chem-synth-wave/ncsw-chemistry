""" The ``ncsw_chemistry.utility.reaction`` package ``compound_standardization`` module. """

from copy import deepcopy
from typing import Iterable

from rdkit.Chem.rdChemReactions import ChemicalReaction

from ncsw_chemistry.utility.compound.standardization import CompoundStandardizationUtility


class ReactionCompoundStandardizationUtility:
    """ The chemical reaction compound standardization utility class. """

    @staticmethod
    def sanitize_compounds(
            reaction_rxn: ChemicalReaction,
            sanitization_operation_keys: Iterable[str] = None,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Sanitize the chemical reaction compounds.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter sanitization_operation_keys: The keys of the chemical reaction compound sanitization operations that
            should be performed. The value `None` indicates that all chemical reaction compound sanitization operations
            should be performed.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit ChemicalReaction`
            object should be constructed and modified.

        :returns: The sanitized chemical reaction.
        """

        if deep_copy:
            reaction_rxn = deepcopy(
                x=reaction_rxn
            )

        for reaction_compound_mols in (
            reaction_rxn.GetReactants(),
            reaction_rxn.GetAgents(),
            reaction_rxn.GetProducts(),
        ):
            for reaction_compound_mol in reaction_compound_mols:
                CompoundStandardizationUtility.sanitize(
                    compound_mol=reaction_compound_mol,
                    sanitization_operation_keys=sanitization_operation_keys,
                    deep_copy=False
                )

        return reaction_rxn
