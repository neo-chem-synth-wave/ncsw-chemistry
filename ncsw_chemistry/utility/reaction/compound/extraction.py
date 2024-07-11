""" The ``ncsw_chemistry.utility.reaction.compound`` package ``extraction`` module. """

from typing import List, Tuple, Union

from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import Mol

from ncsw_chemistry.utility.compound.atom_map_number import CompoundAtomMapNumberUtility
from ncsw_chemistry.utility.compound.format_conversion import CompoundFormatConversionUtility


ReactionCompoundsTuple = Tuple[
    List[Tuple[str, Mol, str, Mol]],
    List[Tuple[str, Mol, str, Mol]],
    List[Tuple[str, Mol, str, Mol]],
]


class ReactionCompoundExtractionUtility:
    """ The chemical reaction compound extraction utility class. """

    @staticmethod
    def _extract_compounds_from_reaction_smiles(
            reaction_smiles: str
    ) -> ReactionCompoundsTuple:
        """
        Extract the chemical reaction compounds from a `SMILES` string.

        :parameter reaction_smiles: The chemical reaction `SMILES` string.

        :returns: The chemical reaction compounds.
        """

        reaction_compounds = (
            list(),
            list(),
            list(),
        )

        for reaction_compounds_index, reaction_compounds_smiles in enumerate(
            iterable=reaction_smiles.split(
                sep=">"
            )
        ):
            if reaction_compounds_smiles != "":
                for reaction_compound_smiles in reaction_compounds_smiles.split()[0].split(
                    sep="."
                ):
                    if reaction_compound_smiles != "":
                        reaction_compound_mol = CompoundFormatConversionUtility.convert_smiles_to_mol(
                            compound_smiles=reaction_compound_smiles
                        )

                        unmapped_reaction_compound_mol = CompoundAtomMapNumberUtility.remove_atom_map_numbers(
                            compound_mol=reaction_compound_mol
                        )

                        unmapped_reaction_compound_smiles = CompoundFormatConversionUtility.convert_mol_to_smiles(
                            compound_mol=unmapped_reaction_compound_mol
                        )

                        reaction_compounds[reaction_compounds_index].append((
                            reaction_compound_smiles,
                            reaction_compound_mol,
                            unmapped_reaction_compound_smiles,
                            unmapped_reaction_compound_mol,
                        ))

        return reaction_compounds

    @staticmethod
    def _extract_compounds_from_reaction_rxn(
            reaction_rxn: ChemicalReaction
    ) -> ReactionCompoundsTuple:
        """
        Extract the chemical reaction compounds from a `RDKit ChemicalReaction` object.

        :parameter reaction_rxn: The chemical reaction `RDKit ChemicalReaction` object.

        :returns: The chemical reaction compounds.
        """

        reaction_compounds = (
            list(),
            list(),
            list(),
        )

        for reaction_compounds_index, reaction_compound_mols in enumerate(
            iterable=(
                reaction_rxn.GetReactants(),
                reaction_rxn.GetAgents(),
                reaction_rxn.GetProducts(),
            )
        ):
            for reaction_compound_mol in reaction_compound_mols:
                reaction_compound_smiles = CompoundFormatConversionUtility.convert_mol_to_smiles(
                    compound_mol=reaction_compound_mol
                )

                unmapped_reaction_compound_mol = CompoundAtomMapNumberUtility.remove_atom_map_numbers(
                    compound_mol=reaction_compound_mol
                )

                unmapped_reaction_compound_smiles = CompoundFormatConversionUtility.convert_mol_to_smiles(
                    compound_mol=unmapped_reaction_compound_mol
                )

                reaction_compounds[reaction_compounds_index].append((
                    reaction_compound_smiles,
                    reaction_compound_mol,
                    unmapped_reaction_compound_smiles,
                    unmapped_reaction_compound_mol,
                ))

        return reaction_compounds

    @staticmethod
    def extract_compounds(
            reaction: Union[str, ChemicalReaction]
    ) -> ReactionCompoundsTuple:
        """
        Extract the chemical reaction compounds.

        :parameter reaction: The chemical reaction `SMILES` string or `RDKit ChemicalReaction` object.

        :returns: The chemical reaction compounds.
        """

        if isinstance(reaction, str):
            return ReactionCompoundExtractionUtility._extract_compounds_from_reaction_smiles(
                reaction_smiles=reaction
            )

        else:
            return ReactionCompoundExtractionUtility._extract_compounds_from_reaction_rxn(
                reaction_rxn=reaction
            )
