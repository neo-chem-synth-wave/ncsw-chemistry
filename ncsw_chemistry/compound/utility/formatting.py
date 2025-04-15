""" The ``ncsw_chemistry.compound.utility`` package ``formatting`` module. """

from typing import Optional

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles, MolToSmarts, MolToSmiles

from ncsw_chemistry.compound.utility.atom import CompoundAtomUtility


class CompoundFormattingUtility:
    """ The chemical compound formatting utility class. """

    @staticmethod
    def convert_compound_mol_to_smarts(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound RDKit Mol object to a SMARTS string.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter remove_atom_map_numbers: The indicator of whether the atom map numbers should be removed from the
            chemical compound.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToSmarts` }.

        :returns: The SMARTS string of the chemical compound.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomUtility.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToSmarts(
            mol=compound_mol,
            **kwargs
        )

    @staticmethod
    def convert_compound_mol_to_smiles(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound RDKit Mol object to a SMILES string.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter remove_atom_map_numbers: The indicator of whether the atom map numbers should be removed from the
            chemical compound.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToSmiles` }.

        :returns: The SMILES string of the chemical compound.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomUtility.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToSmiles(
            mol=compound_mol,
            **kwargs
        )

    @staticmethod
    def convert_compound_smarts_to_mol(
            compound_smarts: str,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound SMARTS string to a RDKit Mol object.

        :parameter compound_smarts: The SMARTS string of the chemical compound.
        :parameter remove_atom_map_numbers: The indicator of whether the atom map numbers should be removed from the
            chemical compound.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromSmarts` }.

        :returns: The RDKit Mol object of the chemical compound.
        """

        compound_mol = MolFromSmarts(
            SMARTS=compound_smarts,
            **kwargs
        )

        if remove_atom_map_numbers and compound_mol is not None:
            return CompoundAtomUtility.remove_atom_map_numbers(
                compound_mol=compound_mol,
                deep_copy=False
            )

        return compound_mol

    @staticmethod
    def convert_compound_smiles_to_mol(
            compound_smiles: str,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound SMILES string to a RDKit Mol object.

        :parameter compound_smiles: The SMILES string of the chemical compound.
        :parameter remove_atom_map_numbers: The indicator of whether the atom map numbers should be removed from the
            chemical compound.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromSmiles` }.

        :returns: The RDKit Mol object of the chemical compound.
        """

        compound_mol = MolFromSmiles(
            SMILES=compound_smiles,
            **kwargs
        )

        if remove_atom_map_numbers and compound_mol is not None:
            return CompoundAtomUtility.remove_atom_map_numbers(
                compound_mol=compound_mol,
                deep_copy=False
            )

        return compound_mol
