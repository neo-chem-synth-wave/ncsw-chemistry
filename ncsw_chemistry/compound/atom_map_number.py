""" The ``ncsw_chemistry.compound`` package ``atom_map_number`` module. """

from copy import deepcopy
from typing import Container, Optional

from rdkit.Chem.rdchem import Mol


class CompoundAtomMapNumberUtilities:
    """ The chemical compound atom map number utilities class. """

    @staticmethod
    def remove_atom_map_numbers(
            compound_mol: Mol,
            atom_indices: Optional[Container[int]] = None,
            deep_copy: bool = True
    ) -> Mol:
        """
        Remove the chemical compound atom map numbers.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter atom_indices: The indices of the chemical compound atoms from which the map numbers should be
            removed. The value `None` indicates that the map numbers should be removed from all chemical compound
            atoms.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound `RDKit Mol` object should
            be constructed and modified.

        :returns: The chemical compound without the atom map numbers.
        """

        if deep_copy:
            compound_mol = deepcopy(
                x=compound_mol
            )

        for atom in compound_mol.GetAtoms():
            if atom_indices is None or atom.GetIdx() in atom_indices:
                atom.ClearProp(
                    key="molAtomMapNumber"
                )

        return compound_mol
