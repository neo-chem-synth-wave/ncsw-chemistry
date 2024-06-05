""" The ``ncsw_chemistry.compound`` package ``format_conversion`` module. """

from typing import Optional

from rdkit.Chem.inchi import MolFromInchi, MolToInchi, MolToInchiKey
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import (
    MolFromMolBlock,
    MolFromSmarts,
    MolFromSmiles,
    MolToMolBlock,
    MolToSmarts,
    MolToSmiles,
)

from ncsw_chemistry.compound.atom.map_number import CompoundAtomMapNumberUtilities


class CompoundFormatConversionUtilities:
    """ The chemical compound format conversion utilities class. """

    @staticmethod
    def convert_inchi_to_mol(
            compound_inchi: str,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `InChI` string to a `RDKit Mol` object.

        :parameter compound_inchi: The chemical compound `InChI` string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.inchi.MolFromInchi` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        compound_mol = MolFromInchi(compound_inchi, **kwargs)

        if remove_atom_map_numbers and compound_mol is not None:
            return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol,
                deep_copy=False
            )

        return compound_mol

    @staticmethod
    def convert_mol_block_to_mol(
            compound_mol_block: str,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `MDL Mol` block string to a `RDKit Mol` object.

        :parameter compound_mol_block: The chemical compound `MDL Mol` block string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromMolBlock` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        compound_mol = MolFromMolBlock(compound_mol_block, **kwargs)

        if remove_atom_map_numbers and compound_mol is not None:
            return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol,
                deep_copy=False
            )

        return compound_mol

    @staticmethod
    def convert_mol_to_inchi(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `InChI` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.inchi.MolToInchi` }.

        :returns: The chemical compound `InChI` string.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToInchi(compound_mol, **kwargs)

    @staticmethod
    def convert_mol_to_inchi_key(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `InChI Key` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.inchi.MolToInchiKey` }.

        :returns: The chemical compound `InChI Key` string.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToInchiKey(compound_mol, **kwargs)

    @staticmethod
    def convert_mol_to_mol_block(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `MDL Mol` block string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToMolBlock` }.

        :returns: The chemical compound `MDL Mol` block string.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToMolBlock(compound_mol, **kwargs)

    @staticmethod
    def convert_mol_to_smarts(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `SMARTS` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToSmarts` }.

        :returns: The chemical compound `SMARTS` string.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToSmarts(compound_mol, **kwargs)

    @staticmethod
    def convert_mol_to_smiles(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `SMILES` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToSmiles` }.

        :returns: The chemical compound `SMILES` string.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        return MolToSmiles(compound_mol, **kwargs)

    @staticmethod
    def convert_smarts_to_mol(
            compound_smarts: str,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `SMARTS` string to a `RDKit Mol` object.

        :parameter compound_smarts: The chemical compound `SMARTS` string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromSmarts` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        compound_mol = MolFromSmarts(compound_smarts, **kwargs)

        if remove_atom_map_numbers and compound_mol is not None:
            return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol,
                deep_copy=False
            )

        return compound_mol

    @staticmethod
    def convert_smiles_to_mol(
            compound_smiles: str,
            remove_atom_map_numbers: bool = False,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `SMILES` string to a `RDKit Mol` object.

        :parameter compound_smiles: The chemical compound `SMILES` string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromSmiles` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        compound_mol = MolFromSmiles(compound_smiles, **kwargs)

        if remove_atom_map_numbers and compound_mol is not None:
            return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                compound_mol=compound_mol,
                deep_copy=False
            )

        return compound_mol
