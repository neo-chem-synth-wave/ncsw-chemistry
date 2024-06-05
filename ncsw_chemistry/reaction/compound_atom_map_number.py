""" The ``ncsw_chemistry.reaction`` package ``compound_atom_map_number`` module. """

from copy import deepcopy

from rdkit.Chem.rdChemReactions import ChemicalReaction

from ncsw_chemistry.compound.atom_map_number import CompoundAtomMapNumberUtilities


class ReactionCompoundAtomMapNumberUtilities:
    """ The chemical reaction compound atom map number utilities class. """

    @staticmethod
    def remove_compound_atom_map_numbers(
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

        for compound_mols in [
            reaction_rxn.GetReactants(),
            reaction_rxn.GetAgents(),
            reaction_rxn.GetProducts(),
        ]:
            for compound_mol in compound_mols:
                CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

        return reaction_rxn
