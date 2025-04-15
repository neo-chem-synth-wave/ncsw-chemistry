""" The ``ncsw_chemistry.compound.utility`` package ``atom`` module. """

from typing import Container, Dict, Optional, Sequence, Union

from rdkit.Chem.rdchem import Atom, Mol

from ncsw_chemistry.compound.utility.typing_ import CompoundAtomPropertyIDTuple


class CompoundAtomUtility:
    """ The chemical compound atom utility class. """

    @staticmethod
    def remove_atom_map_numbers(
            compound_mol: Mol,
            atom_indices: Optional[Container[int]] = None,
            deep_copy: bool = True
    ) -> Mol:
        """
        Remove the atom map numbers from a chemical compound.

        :parameter compound_mol: The RDKit Mol object of the chemical compound.
        :parameter atom_indices: The indices of the chemical compound atoms from which the map numbers should be
            removed. The value `None` indicates that the map numbers should be removed from all chemical compound atoms.
        :parameter deep_copy: The indicator of whether a deep copy of the chemical compound RDKit Mol object should be
            constructed and modified.

        :returns: The chemical compound without the atom map numbers.
        """

        if deep_copy:
            compound_mol = Mol(compound_mol)

        for atom in compound_mol.GetAtoms():
            if atom_indices is None or atom.GetIdx() in atom_indices:
                atom.ClearProp(
                    key="molAtomMapNumber"
                )

        return compound_mol

    @staticmethod
    def get_atom_properties(
            atom: Atom,
            atom_property_keys: Optional[Sequence[str]] = None
    ) -> Dict[str, Union[bool, int, float, str]]:
        """
        Get the properties of a chemical compound atom.

        :parameter atom: The RDKit Atom object of the chemical compound atom.
        :parameter atom_property_keys: The keys of the chemical compound atom properties that should be retrieved. The
            value `None` indicates that all chemical compound atom properties should be retrieved.

        :returns: The properties of the chemical compound atom.
        """

        atom_property_getters = {
            "atomic_number": atom.GetAtomicNum,
            "chiral_tag": lambda: str(atom.GetChiralTag()),
            "degree": atom.GetDegree,
            "explicit_valence": atom.GetExplicitValence,
            "formal_charge": atom.GetFormalCharge,
            "hybridization": lambda: str(atom.GetHybridization()),
            "implicit_valence": atom.GetImplicitValence,
            "is_aromatic": atom.GetIsAromatic,
            "is_in_ring": atom.IsInRing,
            "isotope": atom.GetIsotope,
            "mass": atom.GetMass,
            "number_of_explicit_hydrogen_atoms": atom.GetNumExplicitHs,
            "number_of_implicit_hydrogen_atoms": atom.GetNumImplicitHs,
            "number_of_radical_electrons": atom.GetNumRadicalElectrons,
            "symbol": atom.GetSymbol,
            "total_degree": atom.GetTotalDegree,
            "total_number_of_hydrogen_atoms": atom.GetTotalNumHs,
            "total_valence": atom.GetTotalValence,
        }

        return {
            atom_property_key: atom_property_getters[atom_property_key]()
            for atom_property_key in (
                atom_property_getters.keys() if atom_property_keys is None else atom_property_keys
            )
        }

    @staticmethod
    def get_atom_property_id(
            atom: Atom,
            atom_property_keys: Optional[Sequence[str]] = None
    ) -> CompoundAtomPropertyIDTuple:
        """
        Get the property ID of a chemical compound atom.

        :parameter atom: The RDKit Atom object of the chemical compound atom.
        :parameter atom_property_keys: The keys of the chemical compound atom properties that should be utilized in the
            property ID. The value `None` indicates that all chemical compound atom properties should be utilized in the
            property ID.

        :returns: The property ID of the chemical compound atom.
        """

        return tuple(
            CompoundAtomUtility.get_atom_properties(
                atom=atom,
                atom_property_keys=atom_property_keys
            ).values()
        )
