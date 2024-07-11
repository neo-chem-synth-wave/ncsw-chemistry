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
from ncsw_chemistry.utility.reaction.compound.sanitization import ReactionCompoundSanitizationUtility
from ncsw_chemistry.utility.reaction.sanitization import ReactionSanitizationUtility


class ReactionFormatConversionUtility:
    """ The chemical reaction format conversion utility class. """

    @staticmethod
    def convert_rxn_block_to_rxn(
            reaction_rxn_block: str,
            remove_compound_atom_map_numbers: bool = False,
            sanitize_compound_mols: bool = False,
            sanitize_reaction_rxn: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction `MDL Rxn` block string to a `RDKit ChemicalReaction` object.

        :parameter reaction_rxn_block: The chemical reaction `MDL Rxn` block string.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter sanitize_compound_mols: The indicator of whether the chemical reaction compound `RDKit Mol` objects
            should be sanitized.
        :parameter sanitize_reaction_rxn: The indicator of whether the chemical reaction `RDKit ChemicalReaction` object
            should be sanitized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromRxnBlock` }.

        :returns: The chemical reaction `RDKit ChemicalReaction` object.
        """

        reaction_rxn = ReactionFromRxnBlock(reaction_rxn_block, **kwargs)

        if remove_compound_atom_map_numbers and reaction_rxn is not None:
            ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        if sanitize_compound_mols and reaction_rxn is not None:
            ReactionCompoundSanitizationUtility.sanitize_compounds(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        if sanitize_reaction_rxn and reaction_rxn is not None:
            ReactionSanitizationUtility.sanitize(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        return reaction_rxn

    @staticmethod
    def convert_rxn_to_rxn_block(
            reaction_rxn: ChemicalReaction,
            remove_compound_atom_map_numbers: bool = False,
            sanitize_compound_mols: bool = False,
            sanitize_reaction_rxn: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical reaction `RDKit ChemicalReaction` object to an `MDL Rxn` block string.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter sanitize_compound_mols: The indicator of whether the chemical reaction compound `RDKit Mol` objects
            should be sanitized.
        :parameter sanitize_reaction_rxn: The indicator of whether the chemical reaction `RDKit ChemicalReaction` object
            should be sanitized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionToRxnBlock` }.

        :returns: The chemical reaction `MDL Rxn` block string.
        """

        if remove_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        if sanitize_compound_mols:
            reaction_rxn = ReactionCompoundSanitizationUtility.sanitize_compounds(
                reaction_rxn=reaction_rxn
            )

        if sanitize_reaction_rxn:
            reaction_rxn = ReactionSanitizationUtility.sanitize(
                reaction_rxn=reaction_rxn
            )

        return ReactionToRxnBlock(reaction_rxn, **kwargs)

    @staticmethod
    def convert_rxn_to_smarts(
            reaction_rxn: ChemicalReaction,
            remove_compound_atom_map_numbers: bool = False,
            sanitize_compound_mols: bool = False,
            sanitize_reaction_rxn: bool = False,
    ) -> Optional[str]:
        """
        Convert a chemical reaction `RDKit ChemicalReaction` object to a `SMARTS` string.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter sanitize_compound_mols: The indicator of whether the chemical reaction compound `RDKit Mol` objects
            should be sanitized.
        :parameter sanitize_reaction_rxn: The indicator of whether the chemical reaction `RDKit ChemicalReaction` object
            should be sanitized.

        :returns: The chemical reaction `SMARTS` string.
        """

        if remove_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        if sanitize_compound_mols:
            reaction_rxn = ReactionCompoundSanitizationUtility.sanitize_compounds(
                reaction_rxn=reaction_rxn
            )

        if sanitize_reaction_rxn:
            reaction_rxn = ReactionSanitizationUtility.sanitize(
                reaction_rxn=reaction_rxn
            )

        return ReactionToSmarts(reaction_rxn)

    @staticmethod
    def convert_rxn_to_smiles(
            reaction_rxn: ChemicalReaction,
            remove_compound_atom_map_numbers: bool = False,
            sanitize_compound_mols: bool = False,
            sanitize_reaction_rxn: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical reaction `RDKit ChemicalReaction` object to a `SMILES` string.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter sanitize_compound_mols: The indicator of whether the chemical reaction compound `RDKit Mol` objects
            should be sanitized.
        :parameter sanitize_reaction_rxn: The indicator of whether the chemical reaction `RDKit ChemicalReaction` object
            should be sanitized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionToSmiles` }.

        :returns: The chemical reaction `SMILES` string.
        """

        if remove_compound_atom_map_numbers:
            reaction_rxn = ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn
            )

        if sanitize_compound_mols:
            reaction_rxn = ReactionCompoundSanitizationUtility.sanitize_compounds(
                reaction_rxn=reaction_rxn
            )

        if sanitize_reaction_rxn:
            reaction_rxn = ReactionSanitizationUtility.sanitize(
                reaction_rxn=reaction_rxn
            )

        return ReactionToSmiles(reaction_rxn, **kwargs)

    @staticmethod
    def convert_smarts_to_rxn(
            reaction_smarts: str,
            remove_compound_atom_map_numbers: bool = False,
            sanitize_compound_mols: bool = False,
            sanitize_reaction_rxn: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction `SMARTS` string to a `RDKit ChemicalReaction` object.

        :parameter reaction_smarts: The chemical reaction `SMARTS` string.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter sanitize_compound_mols: The indicator of whether the chemical reaction compound `RDKit Mol` objects
            should be sanitized.
        :parameter sanitize_reaction_rxn: The indicator of whether the chemical reaction `RDKit ChemicalReaction` object
            should be sanitized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromSmarts` }.

        :returns: The chemical reaction `RDKit ChemicalReaction` object.
        """

        reaction_rxn = ReactionFromSmarts(reaction_smarts, **kwargs)

        if remove_compound_atom_map_numbers and reaction_rxn is not None:
            ReactionCompoundAtomMapNumberUtility.remove_atom_map_numbers(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        if sanitize_compound_mols and reaction_rxn is not None:
            ReactionCompoundSanitizationUtility.sanitize_compounds(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        if sanitize_reaction_rxn and reaction_rxn is not None:
            ReactionSanitizationUtility.sanitize(
                reaction_rxn=reaction_rxn,
                deep_copy=False
            )

        return reaction_rxn

    @staticmethod
    def convert_smiles_to_rxn(
            reaction_smiles: str,
            remove_compound_atom_map_numbers: bool = False,
            sanitize_compound_mols: bool = False,
            sanitize_reaction_rxn: bool = False,
            **kwargs
    ) -> Optional[ChemicalReaction]:
        """
        Convert a chemical reaction `SMILES` string to a `RDKit ChemicalReaction` object.

        :parameter reaction_smiles: The chemical reaction `SMILES` string.
        :parameter remove_compound_atom_map_numbers: The indicator of whether the chemical reaction compound atom map
            numbers should be removed.
        :parameter sanitize_compound_mols: The indicator of whether the chemical reaction compound `RDKit Mol` objects
            should be sanitized.
        :parameter sanitize_reaction_rxn: The indicator of whether the chemical reaction `RDKit ChemicalReaction` object
            should be sanitized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdChemReactions.ReactionFromSmarts` }.

        :returns: The chemical reaction `RDKit ChemicalReaction` object.
        """

        return ReactionFormatConversionUtility.convert_smarts_to_rxn(
            reaction_smarts=reaction_smiles,
            remove_compound_atom_map_numbers=remove_compound_atom_map_numbers,
            sanitize_compound_mols=sanitize_compound_mols,
            sanitize_reaction_rxn=sanitize_reaction_rxn,
            **kwargs
        )
