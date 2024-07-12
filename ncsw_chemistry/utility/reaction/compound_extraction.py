""" The ``ncsw_chemistry.utility.reaction`` package ``compound_extraction`` module. """

from rdkit.Chem.rdChemReactions import ChemicalReaction

from ncsw_chemistry.utility.compound.atom_map_number import CompoundAtomMapNumberUtility
from ncsw_chemistry.utility.compound.format_conversion import CompoundFormatConversionUtility
from ncsw_chemistry.utility.reaction.typing import ReactionCompoundsTuple


class ReactionCompoundExtractionUtility:
    """ The chemical reaction compound extraction utility class. """

    @staticmethod
    def extract_compounds_from_reaction_smiles(
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
    def extract_compounds_from_reaction_rxn(
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
