""" The ``ncsw_chemistry.reaction`` package ``compound_extraction`` module. """

from logging import Logger
from typing import Optional

from ncsw_chemistry.compound.atom_map_number import CompoundAtomMapNumberUtilities
from ncsw_chemistry.compound.format_conversion import CompoundFormatConversionUtilities
from ncsw_chemistry.reaction.typing import ReactionCompoundsTuple


class ReactionCompoundExtractionUtilities:
    """ The chemical reaction compound extraction utilities class. """

    @staticmethod
    def extract_compounds(
            reaction_smiles: str,
            logger: Optional[Logger] = None
    ) -> Optional[ReactionCompoundsTuple]:
        """
        Extract the chemical reaction compounds.

        :parameter reaction_smiles: The chemical reaction `SMILES` string.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.

        :returns: The chemical reaction compounds.
        """

        try:
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
                            reaction_compound_mol = \
                                CompoundFormatConversionUtilities.convert_smiles_to_mol(
                                    compound_smiles=reaction_compound_smiles,
                                    logger=logger
                                )

                            reaction_compound_canonical_smiles = \
                                CompoundFormatConversionUtilities.convert_mol_to_smiles(
                                    compound_mol=reaction_compound_mol,
                                    logger=logger
                                )

                            unmapped_reaction_compound_mol = \
                                CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                                    compound_mol=reaction_compound_mol
                                )

                            unmapped_reaction_compound_canonical_smiles = \
                                CompoundFormatConversionUtilities.convert_mol_to_smiles(
                                    compound_mol=unmapped_reaction_compound_mol,
                                    logger=logger
                                )

                            reaction_compounds[reaction_compounds_index].append((
                                reaction_compound_smiles,
                                reaction_compound_canonical_smiles,
                                reaction_compound_mol,
                                unmapped_reaction_compound_canonical_smiles,
                                unmapped_reaction_compound_mol,
                            ))

            return reaction_compounds

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None
