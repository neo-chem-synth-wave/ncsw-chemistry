""" The ``ncsw_chemistry.utility.compound`` package ``substructure`` module. """

from typing import Dict, Iterable, List, Union

from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
from rdkit.Chem.rdFMCS import FindMCS
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFragmentToSmarts, MolFragmentToSmiles

from ncsw_chemistry.utility.compound.atom_map_number import CompoundAtomMapNumberUtility


class CompoundSubstructureUtility:
    """ The chemical compound substructure utility class. """

    @staticmethod
    def get_substructure_smarts(
            compound_mol: Mol,
            atom_indices: List[int],
            **kwargs
    ) -> str:
        """
        Get the chemical compound substructure `SMARTS` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter atom_indices: The indices of the chemical compound substructure atoms.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFragmentToSmarts` }.

        :returns: The chemical compound substructure `SMARTS` string.
        """

        return MolFragmentToSmarts(
            mol=compound_mol,
            atomsToUse=atom_indices,
            **kwargs
        )

    @staticmethod
    def get_substructure_smiles(
            compound_mol: Mol,
            atom_indices: List[int],
            **kwargs
    ) -> str:
        """
        Get the chemical compound substructure `SMILES` string.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter atom_indices: The indices of the chemical compound substructure atoms.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdmolfiles.MolFragmentToSmiles` }.

        :returns: The chemical compound substructure `SMILES` string.
        """

        return MolFragmentToSmiles(
            mol=compound_mol,
            atomsToUse=atom_indices,
            **kwargs
        )

    @staticmethod
    def get_maximum_common_substructure(
            compound_mols: Iterable[Mol],
            **kwargs
    ) -> Dict[str, Union[bool, int, str, Mol, Dict]]:
        """
        Get the maximum common substructure of the chemical compounds.

        :parameter compound_mols: The chemical compound `RDKit Mol` objects.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdkit.Chem.rdFMCS.FindMCS` }.

        :returns: The maximum common substructure of the chemical compounds.
        """

        maximum_common_substructure = FindMCS(
            mols=compound_mols,
            **kwargs
        )

        return {
            "smarts": maximum_common_substructure.smartsString,
            "mol": maximum_common_substructure.queryMol,
            "number_of_atoms": maximum_common_substructure.numAtoms,
            "number_of_bonds": maximum_common_substructure.numBonds,
            "canceled": maximum_common_substructure.canceled,
            "degenerate_smarts_and_mols": maximum_common_substructure.degenerateSmartsQueryMolDict,
        }

    @staticmethod
    def get_bemis_murcko_scaffold(
            compound_mol: Mol,
            remove_atom_map_numbers: bool = True,
            generic_bemis_murcko_scaffold: bool = True
    ) -> Mol:
        """
        Get the Bemis-Murcko scaffold of a chemical compound.

        :parameter compound_mol: The chemical compound `RDKit Mol` object.
        :parameter remove_atom_map_numbers: The indicator of whether the chemical compound atom map numbers should be
            removed.
        :parameter generic_bemis_murcko_scaffold: The indicator of whether the Bemis-Murcko scaffold should be generic.

        :returns: The Bemis-Murcko scaffold of the chemical compound.
        """

        if remove_atom_map_numbers:
            compound_mol = CompoundAtomMapNumberUtility.remove_atom_map_numbers(
                compound_mol=compound_mol
            )

        bemis_murcko_scaffold_mol = GetScaffoldForMol(
            mol=compound_mol
        )

        if generic_bemis_murcko_scaffold:
            return MakeScaffoldGeneric(
                mol=bemis_murcko_scaffold_mol
            )

        return bemis_murcko_scaffold_mol
