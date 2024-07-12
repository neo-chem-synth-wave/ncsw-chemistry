""" The ``ncsw_chemistry.utility.compound`` package ``bond_property`` module. """

from typing import Container, Dict, FrozenSet, Optional, Tuple, Union

from rdkit.Chem.rdchem import Bond

from ncsw_chemistry.utility.compound.atom_property import CompoundAtomPropertyUtility


class CompoundBondPropertyUtility:
    """ The chemical compound bond property utility class. """

    @staticmethod
    def get_bond_properties(
            bond: Bond
    ) -> Dict[str, Union[bool, str]]:
        """
        Get the chemical compound bond properties.

        :parameter bond: The chemical compound bond `RDKit Bond` object.

        :returns: The chemical compound bond properties.
        """

        return {
            "direction": str(bond.GetBondDir()),
            "is_aromatic": bond.GetIsAromatic(),
            "is_conjugated": bond.GetIsConjugated(),
            "is_in_ring": bond.IsInRing(),
            "stereo_configuration": str(bond.GetStereo()),
            "type": str(bond.GetBondType()),
        }

    @staticmethod
    def construct_bond_property_identification_tag(
            bond: Bond,
            atom_property_keys: Optional[Container[str]] = None,
            bond_property_keys: Optional[Container[str]] = None
    ) -> Tuple[FrozenSet[Tuple[Union[bool, int, float, str], ...]], Tuple[Union[bool, str], ...]]:
        """
        Construct the chemical compound bond property identification tag.

        :parameter bond: The chemical compound bond `RDKit Bond` object.
        :parameter atom_property_keys: The keys of the chemical compound atom properties that should be utilized to
            construct the identification tag. The value `None` indicates that all chemical compound atom properties
            should be utilized to construct the identification tag.
        :parameter bond_property_keys: The keys of the chemical compound bond properties that should be utilized to
            construct the identification tag. The value `None` indicates that all chemical compound bond properties
            should be utilized to construct the identification tag.

        :returns: The chemical compound bond property identification tag.
        """

        return (
            frozenset({
                CompoundAtomPropertyUtility.construct_atom_property_identification_tag(
                    atom=bond.GetBeginAtom(),
                    atom_property_keys=atom_property_keys
                ),
                CompoundAtomPropertyUtility.construct_atom_property_identification_tag(
                    atom=bond.GetEndAtom(),
                    atom_property_keys=atom_property_keys
                ),
            }),
            tuple(
                bond_property_value
                for bond_property_key, bond_property_value in CompoundBondPropertyUtility.get_bond_properties(
                    bond=bond
                ).items() if bond_property_keys is None or bond_property_key in bond_property_keys
            ),
        )
