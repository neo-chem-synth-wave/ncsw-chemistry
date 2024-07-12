""" The ``ncsw_chemistry.utility.reaction`` package ``compound_atom_map_number`` module. """

from copy import deepcopy

from rdkit.Chem.rdChemReactions import ChemicalReaction

from ncsw_chemistry.utility.compound.atom_map_number import CompoundAtomMapNumberUtility


class ReactionCompoundAtomMapNumberUtility:
    """ The chemical reaction compound atom map number utility class. """

    @staticmethod
    def remove_atom_map_numbers(
            reaction_rxn: ChemicalReaction,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Remove the chemical reaction compound atom map numbers.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit ChemicalReaction`
            object should be constructed and modified.

        :returns: The chemical reaction without the compound atom map numbers.
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
                CompoundAtomMapNumberUtility.remove_atom_map_numbers(
                    compound_mol=reaction_compound_mol,
                    deep_copy=False
                )

        return reaction_rxn
