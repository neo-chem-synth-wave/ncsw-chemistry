""" The ``ncsw_chemistry.compound.utility`` package ``substructure`` module. """

from collections import Counter
from typing import Container, Optional, Sequence

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFragmentToSmarts, MolFragmentToSmiles

from ncsw_chemistry.compound.utility.atom import CompoundAtomUtility
from ncsw_chemistry.compound.utility.bond import CompoundBondUtility
from ncsw_chemistry.compound.utility.typing_ import CompoundSubstructurePropertyIDTuple


class CompoundSubstructureUtility:
    """ The chemical compound substructure utility class. """

    @staticmethod
    def get_substructure_property_id(
            compound_mol: Mol,
            substructure_atom_indices: Container[int],
            substructure_atom_property_keys: Optional[Sequence[str]] = None,
            substructure_bond_property_keys: Optional[Sequence[str]] = None
    ) -> CompoundSubstructurePropertyIDTuple:
        """
        Get the property ID of a chemical compound substructure.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter substructure_atom_indices: The indices of the chemical compound substructure atoms.
        :parameter substructure_atom_property_keys: The keys of the chemical compound substructure atom properties that
            should be utilized in the property ID. The value `None` indicates that all chemical compound substructure
            atom properties should be utilized in the property ID.
        :parameter substructure_bond_property_keys: The keys of the chemical compound substructure bond properties that
            should be utilized in the property ID. The value `None` indicates that all chemical compound substructure
            bond properties should be utilized in the property ID.

        :returns: The property ID of the chemical compound substructure.
        """

        substructure_atom_index_to_property_id, substructure_bond_atom_indices_to_property_id = dict(), dict()

        for bond in compound_mol.GetBonds():
            if (
                bond.GetBeginAtomIdx() in substructure_atom_indices and
                bond.GetEndAtomIdx() in substructure_atom_indices
            ):
                if bond.GetBeginAtomIdx() not in substructure_atom_index_to_property_id.keys():
                    substructure_atom_index_to_property_id[
                        bond.GetBeginAtomIdx()
                    ] = CompoundAtomUtility.get_atom_property_id(
                        atom=bond.GetBeginAtom(),
                        atom_property_keys=substructure_atom_property_keys
                    )

                if bond.GetEndAtomIdx() not in substructure_atom_index_to_property_id.keys():
                    substructure_atom_index_to_property_id[
                        bond.GetEndAtomIdx()
                    ] = CompoundAtomUtility.get_atom_property_id(
                        atom=bond.GetEndAtom(),
                        atom_property_keys=substructure_atom_property_keys
                    )

                substructure_bond_atom_indices_to_property_id[
                    frozenset({bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), })
                ] = CompoundBondUtility.get_bond_property_id(
                    bond=bond,
                    bond_atom_property_keys=substructure_atom_property_keys,
                    bond_property_keys=substructure_bond_property_keys
                )

        return (
            frozenset(Counter(substructure_atom_index_to_property_id.values()).items()),
            frozenset(Counter(substructure_bond_atom_indices_to_property_id.values()).items()),
        )

    @staticmethod
    def get_substructure_smarts(
            compound_mol: Mol,
            substructure_atom_indices: Container[int],
            **kwargs
    ) -> Optional[str]:
        """
        Get the SMARTS string of a chemical compound substructure.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter substructure_atom_indices: The indices of the chemical compound substructure atoms.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { ``rdkit.Chem.rdmolfiles.MolFragmentToSmarts`` }.

        :returns: The SMARTS string of the chemical compound substructure.
        """

        return MolFragmentToSmarts(
            mol=compound_mol,
            atomsToUse=substructure_atom_indices,
            **kwargs
        )

    @staticmethod
    def get_substructure_smiles(
            compound_mol: Mol,
            substructure_atom_indices: Container[int],
            **kwargs
    ) -> Optional[str]:
        """
        Get the SMILES string of a chemical compound substructure.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter substructure_atom_indices: The indices of the chemical compound substructure atoms.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { ``rdkit.Chem.rdmolfiles.MolFragmentToSmiles`` }.

        :returns: The SMILES string of the chemical compound substructure.
        """

        return MolFragmentToSmiles(
            mol=compound_mol,
            atomsToUse=substructure_atom_indices,
            **kwargs
        )
