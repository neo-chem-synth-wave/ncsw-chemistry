""" The ``ncsw_chemistry.utility.reaction`` package ``format_conversion`` module. """

from typing import Optional

from rdkit.Chem.rdChemReactions import (
    ChemicalReaction,
    ReactionFromRxnBlock,
    ReactionFromSmarts,
    ReactionToRxnBlock,
    ReactionToSmarts,
    ReactionToSmiles,
)

from ncsw_chemistry.utility.reaction.compound.atom_map_number import ReactionCompoundAtomMapNumberUtility


class ReactionFormatConversionUtility:
    """ The chemical reaction format conversion utility class. """

    @staticmethod
    def convert_rxn_block_to_rxn(
            reaction_rxn_block: str,
            remove_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction `MDL Rxn` block string to a `RDKit ChemicalReaction` object.

        :parameter reaction_rxn_block: The chemical reaction `MDL Rxn` block string.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromRxnBlock` }.

        :returns: The chemical reaction `RDKit ChemicalReaction` object.
        """

        reaction_rxn = ReactionFromRxnBlock(reaction_rxn_block, **kwargs)

        if remove_compound_atom_map_numbers and reaction_rxn is not None:
            return ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        return reaction_rxn

    @staticmethod
    def convert_rxn_to_rxn_block(
            reaction_rxn: ChemicalReaction,
            remove_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical reaction `RDKit ChemicalReaction` object to an `MDL Rxn` block string.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionToRxnBlock` }.

        :returns: The chemical reaction `MDL Rxn` block string.
        """

        if remove_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        return ReactionToRxnBlock(reaction_rxn, **kwargs)

    @staticmethod
    def convert_rxn_to_smarts(
            reaction_rxn: ChemicalReaction,
            remove_compound_atom_map_numbers: bool = False
    ) -> Optional[str]:
        """
        Convert a chemical reaction `RDKit ChemicalReaction` object to a `SMARTS` string.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.

        :returns: The chemical reaction `SMARTS` string.
        """

        if remove_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        return ReactionToSmarts(reaction_rxn)

    @staticmethod
    def convert_rxn_to_smiles(
            reaction_rxn: ChemicalReaction,
            remove_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical reaction `RDKit ChemicalReaction` object to a `SMILES` string.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionToSmiles` }.

        :returns: The chemical reaction `SMILES` string.
        """

        if remove_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        return ReactionToSmiles(reaction_rxn, **kwargs)

    @staticmethod
    def convert_smarts_to_rxn(
            reaction_smarts: str,
            remove_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction `SMARTS` string to a `RDKit ChemicalReaction` object.

        :parameter reaction_smarts: The chemical reaction `SMARTS` string.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromSmarts` }.

        :returns: The chemical reaction `RDKit ChemicalReaction` object.
        """

        reaction_rxn = ReactionFromSmarts(reaction_smarts, **kwargs)

        if remove_compound_atom_map_numbers and reaction_rxn is not None:
            return ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        return reaction_rxn

    @staticmethod
    def convert_smiles_to_rxn(
            reaction_smiles: str,
            remove_compound_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction `SMILES` string to a `RDKit ChemicalReaction` object.

        :parameter reaction_smiles: The chemical reaction `SMILES` string.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromSmarts` }.

        :returns: The chemical reaction `RDKit ChemicalReaction` object.
        """

        return ReactionFormatConversionUtility.convert_smarts_to_rxn(
            reaction_smarts=reaction_smiles,
            remove_compound_atom_map_numbers=remove_compound_atom_map_numbers,
            **kwargs
        )
