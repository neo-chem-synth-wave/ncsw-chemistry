""" The ``ncsw_chemistry.utility.reaction`` package ``sanitization`` module. """

from copy import deepcopy

from rdkit.Chem.rdChemReactions import ChemicalReaction, SanitizeRxn


class ReactionSanitizationUtility:
    """ The chemical reaction sanitization utility class. """

    @staticmethod
    def sanitize_compounds(
            reaction_rxn: ChemicalReaction,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Sanitize a chemical reaction.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit ChemicalReaction`
            object should be constructed and modified.

        :returns: The sanitized chemical reaction.
        """

        if deep_copy:
            reaction_rxn = deepcopy(
                x=reaction_rxn
            )

        SanitizeRxn(reaction_rxn)

        return reaction_rxn
