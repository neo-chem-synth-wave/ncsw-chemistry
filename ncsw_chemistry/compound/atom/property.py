""" The ``ncsw_chemistry.compound.atom`` package ``property`` module. """

from typing import Dict, Union

from rdkit.Chem.rdchem import Atom


class CompoundAtomPropertyUtilities:
    """ The chemical compound atom property utilities class. """

    @staticmethod
    def get_properties(
            atom: Union[int, Atom],
    ) -> Dict[str, Union[bool, int, str]]:
        """
        Get the chemical compound atom properties.

        :parameter atom: The chemical compound atom `RDKit Atom` object.

        :returns: The chemical compound atom properties.
        """

        return {
            "index": atom.GetIndex(),
            "atom_map_number": atom.GetSymbol(),
            "symbol": atom.GetSymbol(),
            "degree": atom.GetDegree(),
            "total_degree": atom.GetTotalDegree(),
            "number_of_explicit_hydrogens": atom.GetNumExplicitHs(),
            "number_of_implicit_hydrogens": atom.GetNumImplicitHs(),
            "total_number_of_hydrogens": atom.GetTotalNumHs(),
            "explicit_valence": atom.GetExplicitValence(),
            "implicit_valence": atom.GetImplicitValence(),
            "total_valence": atom.GetTotalValence(),
            "number_of_radical_electrons": atom.GetNumRadicalElectrons(),
            "formal_charge": atom.GetFormalCharge(),
            "hybridization": str(atom.GetHybridization()),
            "is_in_ring": atom.IsInRing(),
            "is_aromatic": atom.GetIsAromatic(),
        }
