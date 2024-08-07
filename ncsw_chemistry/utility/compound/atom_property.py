""" The ``ncsw_chemistry.utility.compound`` package ``atom_property`` module. """

from typing import Container, Dict, Optional, Tuple, Union

from rdkit.Chem.rdchem import Atom


class CompoundAtomPropertyUtility:
    """ The chemical compound atom property utility class. """

    @staticmethod
    def get_atom_properties(
            atom: Atom
    ) -> Dict[str, Union[bool, int, float, str]]:
        """
        Get the chemical compound atom properties.

        :parameter atom: The chemical compound atom `RDKit Atom` object.

        :returns: The chemical compound atom properties.
        """

        return {
            "atomic_number": atom.GetAtomicNum(),
            "chiral_tag": str(atom.GetChiralTag()),
            "degree": atom.GetDegree(),
            "explicit_valence": atom.GetExplicitValence(),
            "formal_charge": atom.GetFormalCharge(),
            "hybridization": str(atom.GetHybridization()),
            "implicit_valence": atom.GetImplicitValence(),
            "is_aromatic": atom.GetIsAromatic(),
            "is_in_ring": atom.IsInRing(),
            "isotope": atom.GetIsotope(),
            "mass": atom.GetMass(),
            "number_of_explicit_hydrogens": atom.GetNumExplicitHs(),
            "number_of_implicit_hydrogens": atom.GetNumImplicitHs(),
            "number_of_radical_electrons": atom.GetNumRadicalElectrons(),
            "symbol": atom.GetSymbol(),
            "total_degree": atom.GetTotalDegree(),
            "total_number_of_hydrogens": atom.GetTotalNumHs(),
            "total_valence": atom.GetTotalValence(),
        }

    @staticmethod
    def construct_atom_property_identification_tag(
            atom: Atom,
            atom_property_keys: Optional[Container[str]] = None
    ) -> Tuple[Union[bool, int, float, str], ...]:
        """
        Construct the chemical compound atom property identification tag.

        :parameter atom: The chemical compound atom `RDKit Atom` object.
        :parameter atom_property_keys: The keys of the chemical compound atom properties that should be utilized to
            construct the identification tag. The value `None` indicates that all chemical compound atom properties
            should be utilized to construct the identification tag.

        :returns: The chemical compound atom property identification tag.
        """

        return tuple(
            atom_property_value
            for atom_property_key, atom_property_value in CompoundAtomPropertyUtility.get_atom_properties(
                atom=atom
            ).items() if atom_property_keys is None or atom_property_key in atom_property_keys
        )
