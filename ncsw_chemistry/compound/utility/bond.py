""" The ``ncsw_chemistry.compound.utility`` package ``bond`` module. """

from typing import Dict, Optional, Sequence, Union

from rdkit.Chem.rdchem import Bond

from ncsw_chemistry.compound.utility.atom import CompoundAtomUtility
from ncsw_chemistry.compound.utility.typing_ import CompoundBondPropertyIDTuple


class CompoundBondUtility:
    """ The chemical compound bond utility class. """

    @staticmethod
    def get_bond_properties(
            bond: Bond,
            bond_property_keys: Optional[Sequence[str]] = None
    ) -> Dict[str, Union[bool, str]]:
        """
        Get the properties of a chemical compound bond.

        :parameter bond: The RDKit Bond object of the chemical compound bond.
        :parameter bond_property_keys: The keys of the chemical compound bond properties that should be retrieved. The
            value `None` indicates that all chemical compound bond properties should be retrieved.

        :returns: The properties of the chemical compound bond.
        """

        bond_property_getters = {
            "direction": lambda: str(bond.GetBondDir()),
            "is_aromatic": bond.GetIsAromatic,
            "is_conjugated": bond.GetIsConjugated,
            "is_in_ring": bond.IsInRing,
            "stereo_configuration": lambda: str(bond.GetStereo()),
            "type": lambda: str(bond.GetBondType()),
        }

        return {
            bond_property_key: bond_property_getters[bond_property_key]()
            for bond_property_key in (
                bond_property_getters.keys() if bond_property_keys is None else bond_property_keys
            )
        }

    @staticmethod
    def get_bond_property_id(
            bond: Bond,
            bond_atom_property_keys: Optional[Sequence[str]] = None,
            bond_property_keys: Optional[Sequence[str]] = None
    ) -> CompoundBondPropertyIDTuple:
        """
        Get the property ID of a chemical compound bond.

        :parameter bond: The RDKit Bond object of the chemical compound bond.
        :parameter bond_atom_property_keys: The keys of the chemical compound bond atom properties that should be
            utilized in the property ID. The value `None` indicates that all chemical compound bond atom properties
            should be utilized in the property ID.
        :parameter bond_property_keys: The keys of the chemical compound bond properties that should be utilized in the
            property ID. The value `None` indicates that all chemical compound bond properties should be utilized in the
            property ID.

        :returns: The property ID of the chemical compound bond.
        """

        return (
            frozenset({
                CompoundAtomUtility.get_atom_property_id(
                    atom=bond.GetBeginAtom(),
                    atom_property_keys=bond_atom_property_keys
                ),
                CompoundAtomUtility.get_atom_property_id(
                    atom=bond.GetEndAtom(),
                    atom_property_keys=bond_atom_property_keys
                ),
            }),
            tuple(
                CompoundBondUtility.get_bond_properties(
                    bond=bond,
                    bond_property_keys=bond_property_keys
                ).values()
            ),
        )
