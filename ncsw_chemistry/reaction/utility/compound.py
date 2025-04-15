""" The ``ncsw_chemistry.reaction.utility`` package ``compound`` module. """

from typing import List, Tuple

from rdkit.Chem.rdChemReactions import ChemicalReaction

from ncsw_chemistry.compound.utility.atom import CompoundAtomUtility


class ReactionCompoundUtility:
    """ The chemical reaction compound utility class. """

    @staticmethod
    def remove_compound_atom_map_numbers(
            reaction_rxn: ChemicalReaction,
            deep_copy: bool = True
    ) -> ChemicalReaction:
        """
        Remove the compound atom map numbers from a chemical reaction.

        :parameter reaction_rxn: The RDKit ChemicalReaction object of the chemical reaction.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical reaction RDKit ChemicalReaction
            object should be constructed and modified.

        :returns: The chemical reaction without the compound atom map numbers.
        """

        if deep_copy:
            reaction_rxn = ChemicalReaction(reaction_rxn)

        for compound_mols in (
            reaction_rxn.GetReactants(),
            reaction_rxn.GetAgents(),
            reaction_rxn.GetProducts(),
        ):
            for compound_mol in compound_mols:
                CompoundAtomUtility.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

        return reaction_rxn

    @staticmethod
    def extract_compound_smiles_or_smarts(
            reaction_smiles_or_smarts: str
    ) -> Tuple[List[str], List[str], List[str]]:
        """
        Extract the compound SMILES or SMARTS strings from a chemical reaction.

        :parameter reaction_smiles_or_smarts: The SMILES or SMARTS string of the chemical reaction.

        :returns: The compound SMILES or SMARTS strings from the chemical reaction.
        """

        compound_smiles_or_smarts_strings = (list(), list(), list(), )

        for compounds_index, compounds_smiles_or_smarts in enumerate(reaction_smiles_or_smarts.split(
            sep=">"
        )):
            if compounds_smiles_or_smarts != "":
                for compound_smiles_or_smarts in compounds_smiles_or_smarts.split()[0].split(
                    sep="."
                ):
                    if compound_smiles_or_smarts != "":
                        compound_smiles_or_smarts_strings[compounds_index].append(
                            compound_smiles_or_smarts
                        )

        return compound_smiles_or_smarts_strings
