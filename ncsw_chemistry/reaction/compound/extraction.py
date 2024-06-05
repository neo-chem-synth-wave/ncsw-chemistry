""" The ``ncsw_chemistry.reaction.compound`` package ``extraction`` module. """

from typing import List, Optional, Tuple

from rdkit.Chem.rdchem import Mol

from ncsw_chemistry.compound.format_conversion import CompoundFormatConversionUtilities
from ncsw_chemistry.compound.atom.map_number import CompoundAtomMapNumberUtilities


ReactionCompoundsTuple = Tuple[
    List[Tuple[str, Optional[str], Optional[Mol], Optional[str], Optional[Mol]]],
    List[Tuple[str, Optional[str], Optional[Mol], Optional[str], Optional[Mol]]],
    List[Tuple[str, Optional[str], Optional[Mol], Optional[str], Optional[Mol]]],
]


class ReactionCompoundExtractionUtilities:
    """ The chemical reaction compound extraction utilities class. """

    @staticmethod
    def extract_compounds(
            reaction_smiles: str
    ) -> Optional[ReactionCompoundsTuple]:
        """
        Extract the chemical reaction compounds.

        :parameter reaction_smiles: The chemical reaction `SMILES` string.

        :returns: The chemical reaction compounds.
        """

        reaction_compounds = (
            list(),
            list(),
            list(),
        )

        for compounds_index, compounds_smiles in enumerate(
            iterable=reaction_smiles.split(
                sep=">"
            )
        ):
            if compounds_smiles != "":
                for compound_smiles in compounds_smiles.split()[0].split(
                    sep="."
                ):
                    if compound_smiles != "":
                        compound_mol = CompoundFormatConversionUtilities.convert_smiles_to_mol(
                            compound_smiles=compound_smiles
                        )

                        compound_canonical_smiles = CompoundFormatConversionUtilities.convert_mol_to_smiles(
                            compound_mol=compound_mol
                        )

                        unmapped_compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                            compound_mol=compound_mol
                        )

                        unmapped_compound_canonical_smiles = CompoundFormatConversionUtilities.convert_mol_to_smiles(
                            compound_mol=unmapped_compound_mol,
                        )

                        reaction_compounds[compounds_index].append((
                            compound_smiles,
                            compound_canonical_smiles,
                            compound_mol,
                            unmapped_compound_canonical_smiles,
                            unmapped_compound_mol,
                        ))

        return reaction_compounds
