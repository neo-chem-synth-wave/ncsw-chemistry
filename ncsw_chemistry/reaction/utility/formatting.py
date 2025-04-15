""" The ``ncsw_chemistry.reaction.utility`` package ``formatting`` module. """

from typing import Optional

from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts, ReactionToSmarts, ReactionToSmiles

from ncsw_chemistry.reaction.utility.compound import ReactionCompoundUtility


class ReactionFormattingUtility:
    """ The chemical reaction formatting utility class. """

    @staticmethod
    def convert_reaction_rxn_to_smarts(
            reaction_rxn: ChemicalReaction,
            remove_reaction_compound_atom_map_numbers: bool = False
    ) -> Optional[str]:
        """
        Convert a chemical reaction RDKit ChemicalReaction object to a SMARTS string.

        :parameter reaction_rxn: The RDKit ChemicalReaction object of the chemical reaction.
        :parameter remove_reaction_compound_atom_map_numbers: The indicator of whether the chemical reaction compound
            atom map numbers should be removed.

        :returns: The SMARTS string of the chemical reaction.
        """

        if remove_reaction_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundUtility.remove_compound_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        return ReactionToSmarts(
            reaction=reaction_rxn
        )

    @staticmethod
    def convert_reaction_rxn_to_smiles(
            reaction_rxn: ChemicalReaction,
            remove_reaction_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical reaction RDKit ChemicalReaction object to a SMILES string.

        :parameter reaction_rxn: The RDKit ChemicalReaction object of the chemical reaction.
        :parameter remove_reaction_compound_atom_map_numbers: The indicator of whether the chemical reaction compound
            atom map numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionToSmiles` }.

        :returns: The SMILES string of the chemical reaction.
        """

        if remove_reaction_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundUtility.remove_compound_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        return ReactionToSmiles(
            reaction=reaction_rxn,
            **kwargs
        )

    @staticmethod
    def convert_reaction_smarts_to_rxn(
            reaction_smarts: str,
            remove_reaction_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction SMARTS string to a RDKit ChemicalReaction object.

        :parameter reaction_smarts: The SMARTS string of the chemical reaction.
        :parameter remove_reaction_compound_atom_map_numbers: The indicator of whether the chemical reaction compound
            atom map numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromSmarts` }.

        :returns: The RDKit ChemicalReaction object of the chemical reaction.
        """

        reaction_rxn = ReactionFromSmarts(
            SMARTS=reaction_smarts,
            **kwargs
        )

        if remove_reaction_compound_atom_map_numbers and reaction_rxn is not None:
            return ReactionCompoundUtility.remove_compound_atom_map_numbers(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        return reaction_rxn

    @staticmethod
    def convert_reaction_smiles_to_rxn(
            reaction_smiles: str,
            remove_reaction_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction SMILES string to a RDKit ChemicalReaction object.

        :parameter reaction_smiles: The SMILES string of the chemical reaction.
        :parameter remove_reaction_compound_atom_map_numbers: The indicator of whether the chemical reaction compound
            atom map numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromSmarts` }.

        :returns: The RDKit ChemicalReaction object of the chemical reaction.
        """

        return ReactionFormattingUtility.convert_reaction_smarts_to_rxn(
            reaction_smarts=reaction_smiles,
            remove_reaction_compound_atom_map_numbers=remove_reaction_compound_atom_map_numbers,
            **kwargs
        )
