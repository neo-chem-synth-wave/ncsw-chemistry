""" The ``ncsw_chemistry.utility.reaction.compound`` package ``extraction`` module. """

from typing import List, Tuple

from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import Mol

from ncsw_chemistry.utility.compound.atom_map_number import CompoundAtomMapNumberUtility
from ncsw_chemistry.utility.compound.format_conversion import CompoundFormatConversionUtility


ReactionCompoundsTuple = Tuple[
    List[Tuple[Mol, str, Mol, str]],
    List[Tuple[Mol, str, Mol, str]],
    List[Tuple[Mol, str, Mol, str]],
]


class ReactionCompoundExtractionUtility:
    """ The chemical reaction compound extraction utility class. """

    @staticmethod
    def extract_compounds(
            reaction_rxn: ChemicalReaction
    ) -> ReactionCompoundsTuple:
        """
        Extract the chemical reaction compounds.

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
                    reaction_compound_mol,
                    reaction_compound_smiles,
                    unmapped_reaction_compound_mol,
                    unmapped_reaction_compound_smiles,
                ))

        return reaction_compounds
