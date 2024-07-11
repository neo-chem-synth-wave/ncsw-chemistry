""" The ``ncsw_chemistry.utility.reaction.compound`` package ``sanitization`` module. """

from copy import deepcopy

from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdmolops import SanitizeMol


class ReactionCompoundSanitizationUtility:
    """ The chemical reaction compound sanitization utility class. """

    @staticmethod
    def sanitize_compounds(
            reaction_rxn: ChemicalReaction,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Sanitize the chemical reaction compounds.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit ChemicalReaction`
            object should be constructed and modified.

        :returns: The chemical reaction with sanitized chemical compounds.
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
                SanitizeMol(reaction_compound_mol)

        return reaction_rxn
