""" The ``ncsw_chemistry.compound`` package ``format_conversion`` module. """

from logging import Logger
from typing import Optional

from rdkit.Chem.inchi import MolFromInchi, MolToInchi, MolToInchiKey
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import (
    MolFromMol2Block,
    MolFromMolBlock,
    MolFromMrvBlock,
    MolFromSmarts,
    MolFromSmiles,
    MolToCXSmarts,
    MolToCXSmiles,
    MolToMolBlock,
    MolToMrvBlock,
    MolToSmarts,
    MolToSmiles,
    MolToV3KMolBlock,
)

from ncsw_chemistry.compound.atom_map_number import CompoundAtomMapNumberUtilities


class CompoundFormatConversionUtilities:
    """ The chemical compound format conversion utilities class. """

    @staticmethod
    def convert_inchi_to_mol(
            compound_inchi: str,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `InChI` string to a `RDKit Mol` object.

        :parameter compound_inchi: The chemical compound `InChI` string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.inchi.MolFromInchi` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        try:
            compound_mol = MolFromInchi(compound_inchi, **kwargs)

            if remove_atom_map_numbers and compound_mol is not None:
                return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

            return compound_mol

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_block_to_mol(
            compound_mol_block: str,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `MDL Mol` block string to a `RDKit Mol` object.

        :parameter compound_mol_block: The chemical compound `MDL Mol` block string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromMolBlock` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        try:
            compound_mol = MolFromMolBlock(compound_mol_block, **kwargs)

            if remove_atom_map_numbers and compound_mol is not None:
                return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

            return compound_mol

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol2_block_to_mol(
            compound_mol2_block: str,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `Tripos Mol2` block string to a `RDKit Mol` object.

        :parameter compound_mol2_block: The chemical compound `Tripos Mol2` block string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromMol2Block` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        try:
            compound_mol = MolFromMol2Block(compound_mol2_block, **kwargs)

            if remove_atom_map_numbers and compound_mol is not None:
                return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

            return compound_mol

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mrv_block_to_mol(
            compound_mrv_block: str,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `Marvin` block string to a `RDKit Mol` object.

        :parameter compound_mrv_block: The chemical compound `Marvin` block string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromMrvBlock` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        try:
            compound_mol = MolFromMrvBlock(compound_mrv_block, **kwargs)

            if remove_atom_map_numbers and compound_mol is not None:
                return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

            return compound_mol

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_cxsmarts(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `CXSMARTS` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToCXSmarts` }.

        :returns: The chemical compound `CXSMARTS` string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToCXSmarts(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_cxsmiles(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `CXSMILES` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToCXSmiles` }.

        :returns: The chemical compound `CXSMILES` string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToCXSmiles(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_inchi(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `InChI` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.inchi.MolToInchi` }.

        :returns: The chemical compound `InChI` string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToInchi(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_inchi_key(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `InChI Key` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.inchi.MolToInchiKey` }.

        :returns: The chemical compound `InChI Key` string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToInchiKey(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_mol_block(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `MDL Mol` block string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToMolBlock` }.

        :returns: The chemical compound `MDL Mol` block string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToMolBlock(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_mrv_block(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `Marvin` block string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToMrvBlock` }.

        :returns: The chemical compound `Marvin` block string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToMrvBlock(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_smarts(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `SMARTS` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToSmarts` }.

        :returns: The chemical compound `SMARTS` string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToSmarts(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_smiles(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to a `SMILES` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToSmiles` }.

        :returns: The chemical compound `SMILES` string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToSmiles(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_smarts_to_mol(
            compound_smarts: str,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `SMARTS` string to a `RDKit Mol` object.

        :parameter compound_smarts: The chemical compound `SMARTS` string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromSmarts` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        try:
            compound_mol = MolFromSmarts(compound_smarts, **kwargs)

            if remove_atom_map_numbers and compound_mol is not None:
                return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

            return compound_mol

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_smiles_to_mol(
            compound_smiles: str,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[Mol]:
        """
        Convert a chemical compound `SMILES` string to a `RDKit Mol` object.

        :parameter compound_smiles: The chemical compound `SMILES` string.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFromSmiles` }.

        :returns: The chemical compound `RDKit Mol` object.
        """

        try:
            compound_mol = MolFromSmiles(compound_smiles, **kwargs)

            if remove_atom_map_numbers and compound_mol is not None:
                return CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol,
                    deep_copy=False
                )

            return compound_mol

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None

    @staticmethod
    def convert_mol_to_v3000_mol_block(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = False,
            logger: Optional[Logger] = None,
            **kwargs
    ) -> Optional[str]:
        """
        Convert a chemical compound `RDKit Mol` object to an `MDL V3000 Mol` block string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter logger: The logger. The value `None` indicates that the logger should not be utilized.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolToV3KMolBlock` }.

        :returns: The chemical compound `MDL V3000 Mol` block string.
        """

        try:
            if remove_atom_map_numbers:
                compound_mol = CompoundAtomMapNumberUtilities.remove_atom_map_numbers(
                    compound_mol=compound_mol
                )

            return MolToV3KMolBlock(compound_mol, **kwargs)

        except Exception as exception_handle:
            if logger is not None:
                logger.debug(
                    msg=exception_handle,
                    exc_info=True
                )

            return None
