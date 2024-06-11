""" The ``ncsw_chemistry.utility.reaction`` package ``reactivity`` module. """

from itertools import chain
from typing import Container, Dict, Optional, Sequence, Set, Tuple

from rdkit.Chem.rdchem import Mol

from ncsw_chemistry.utility.compound.atom.property import CompoundAtomPropertyUtility
from ncsw_chemistry.utility.compound.bond.property import CompoundBondPropertyUtility


class ReactionReactivityUtility:
    """ The chemical reaction reactivity utility class. """

    @staticmethod
    def get_synthon_atom_map_numbers(
            mapped_reactant_compound_mol: Mol,
            mapped_product_compound_mol: Mol,
            atom_property_keys: Optional[Container[str]] = None,
            bond_property_keys: Optional[Container[str]] = None
    ) -> Set[int]:
        """
        Get the chemical reaction reactant and product compound synthon atom map numbers.

        :parameter mapped_reactant_compound_mol: The mapped chemical reaction reactant compound `RDKit Mol` object.
        :parameter mapped_product_compound_mol: The mapped chemical reaction product compound `RDKit Mol` object.
        :parameter atom_property_keys: The keys of the chemical compound atom properties that should be utilized to
            construct the identification tags. The value `None` indicates that all chemical compound atom properties
            should be utilized to construct the identification tags.
        :parameter bond_property_keys: The keys of the chemical compound bond properties that should be utilized to
            construct the identification tags. The value `None` indicates that all chemical compound bond properties
            should be utilized to construct the identification tags.

        :returns: The chemical reaction reactant and product compound synthon atom map numbers.
        """

        reactant_bond_atom_map_numbers = set()

        for bond in mapped_reactant_compound_mol.GetBonds():
            if bond.GetBeginAtom().HasProp(
                key="molAtomMapNumber"
            ) and bond.GetEndAtom().HasProp(
                key="molAtomMapNumber"
            ):
                reactant_bond_atom_map_numbers.add(
                    frozenset({
                        bond.GetBeginAtom().GetAtomMapNum(),
                        bond.GetEndAtom().GetAtomMapNum(),
                    })
                )

        product_bond_atom_map_numbers = set()

        for bond in mapped_product_compound_mol.GetBonds():
            if bond.GetBeginAtom().HasProp(
                key="molAtomMapNumber"
            ) and bond.GetEndAtom().HasProp(
                key="molAtomMapNumber"
            ):
                product_bond_atom_map_numbers.add(
                    frozenset({
                        bond.GetBeginAtom().GetAtomMapNum(),
                        bond.GetEndAtom().GetAtomMapNum(),
                    })
                )

        non_synthon_bond_atom_map_numbers = set(
            chain.from_iterable(
                reactant_bond_atom_map_numbers.symmetric_difference(
                    product_bond_atom_map_numbers
                )
            )
        )

        synthon_bond_atom_map_numbers = set()

        for reactant_bond in mapped_reactant_compound_mol.GetBonds():
            if reactant_bond.GetBeginAtom().HasProp(
                key="molAtomMapNumber"
            ) and reactant_bond.GetEndAtom().HasProp(
                key="molAtomMapNumber"
            ):
                for product_bond in mapped_product_compound_mol.GetBonds():
                    if product_bond.GetBeginAtom().HasProp(
                        key="molAtomMapNumber"
                    ) and product_bond.GetEndAtom().HasProp(
                        key="molAtomMapNumber"
                    ):
                        if {
                            reactant_bond.GetBeginAtom().GetAtomMapNum(),
                            reactant_bond.GetEndAtom().GetAtomMapNum(),
                        } == {
                            product_bond.GetBeginAtom().GetAtomMapNum(),
                            product_bond.GetEndAtom().GetAtomMapNum(),
                        }:
                            if CompoundBondPropertyUtility.construct_property_identification_tag(
                                bond=reactant_bond,
                                atom_property_keys=atom_property_keys,
                                bond_property_keys=bond_property_keys
                            ) == CompoundBondPropertyUtility.construct_property_identification_tag(
                                bond=product_bond,
                                atom_property_keys=atom_property_keys,
                                bond_property_keys=bond_property_keys
                            ):
                                synthon_bond_atom_map_numbers.add(
                                    reactant_bond.GetBeginAtom().GetAtomMapNum()
                                    # product_bond.GetBeginAtom().GetAtomMapNum() is equivalent.
                                )

                                synthon_bond_atom_map_numbers.add(
                                    reactant_bond.GetEndAtom().GetAtomMapNum()
                                    # product_bond.GetEndAtom().GetAtomMapNum() is equivalent.
                                )

                            else:
                                non_synthon_bond_atom_map_numbers.add(
                                    reactant_bond.GetBeginAtom().GetAtomMapNum()
                                    # product_bond.GetBeginAtom().GetAtomMapNum() is equivalent.
                                )

                                non_synthon_bond_atom_map_numbers.add(
                                    reactant_bond.GetEndAtom().GetAtomMapNum()
                                    # product_bond.GetEndAtom().GetAtomMapNum() is equivalent.
                                )

        synthon_atom_map_numbers, non_synthon_atom_map_numbers = set(), set()

        for reactant_atom in mapped_reactant_compound_mol.GetAtoms():
            if reactant_atom.HasProp(
                key="molAtomMapNumber"
            ) and reactant_atom.GetAtomMapNum() not in synthon_bond_atom_map_numbers and \
                    reactant_atom.GetAtomMapNum() not in non_synthon_bond_atom_map_numbers:
                for product_atom in mapped_product_compound_mol.GetAtoms():
                    if product_atom.HasProp(
                        key="molAtomMapNumber"
                    ) and product_atom.GetAtomMapNum() not in synthon_bond_atom_map_numbers and \
                            product_atom.GetAtomMapNum() not in non_synthon_bond_atom_map_numbers:
                        if reactant_atom.GetAtomMapNum() == product_atom.GetAtomMapNum():
                            if CompoundAtomPropertyUtility.construct_property_identification_tag(
                                atom=reactant_atom,
                                atom_property_keys=atom_property_keys
                            ) == CompoundAtomPropertyUtility.construct_property_identification_tag(
                                atom=product_atom,
                                atom_property_keys=atom_property_keys
                            ):
                                synthon_atom_map_numbers.add(
                                    reactant_atom.GetAtomMapNum()
                                    # product_atom.GetAtomMapNum() is equivalent.
                                )

                            else:
                                non_synthon_atom_map_numbers.add(
                                    reactant_atom.GetAtomMapNum()
                                    # product_atom.GetAtomMapNum() is equivalent.
                                )

        synthon_atom_map_numbers.update(
            synthon_bond_atom_map_numbers
        )

        non_synthon_atom_map_numbers.update(
            non_synthon_bond_atom_map_numbers
        )

        synthon_atom_map_numbers.difference_update(
            non_synthon_bond_atom_map_numbers
        )

        return synthon_atom_map_numbers

    @staticmethod
    def extract_reactive_sites_and_synthons_using_atom_map_numbers(
            mapped_reactant_compound_mols: Sequence[Mol],
            mapped_product_compound_mols: Sequence[Mol],
            atom_property_keys: Optional[Container[str]] = None,
            bond_property_keys: Optional[Container[str]] = None
    ) -> Dict[int, Tuple[Dict[int, Tuple[Set[int], Dict[int, int]]], Set[int]]]:
        """
        Extract the chemical reaction reactant and product compound reactive sites and synthons using atom map numbers.

        :parameter mapped_reactant_compound_mols: The mapped chemical reaction reactant compound `RDKit Mol` objects.
        :parameter mapped_product_compound_mols: The mapped chemical reaction product compound `RDKit Mol` objects.
        :parameter atom_property_keys: The keys of the chemical compound atom properties that should be utilized to
            construct the identification tags. The value `None` indicates that all chemical compound atom properties
            should be utilized to construct the identification tags.
        :parameter bond_property_keys: The keys of the chemical compound bond properties that should be utilized to
            construct the identification tags. The value `None` indicates that all chemical compound bond properties
            should be utilized to construct the identification tags.

        :returns: The chemical reaction reactant and product compound reactive sites and synthons.
        """

        reactant_atom_map_number_to_index_dictionaries, product_atom_map_number_to_index_dictionaries = dict(), dict()

        product_reactive_sites_and_synthons = dict()

        for product_index, mapped_product_compound_mol in enumerate(
            iterable=mapped_product_compound_mols
        ):
            if product_index not in product_atom_map_number_to_index_dictionaries.keys():
                product_atom_map_number_to_index_dictionaries[product_index] = {
                    product_atom.GetAtomMapNum(): product_atom.GetIdx()
                    for product_atom in mapped_product_compound_mol.GetAtoms()
                }

            reactant_reactive_sites_and_synthons = dict()

            for reactant_index, mapped_reactant_compound_mol in enumerate(
                iterable=mapped_reactant_compound_mols
            ):
                if reactant_index not in reactant_atom_map_number_to_index_dictionaries.keys():
                    reactant_atom_map_number_to_index_dictionaries[reactant_index] = {
                        reactant_atom.GetAtomMapNum(): reactant_atom.GetIdx()
                        for reactant_atom in mapped_reactant_compound_mol.GetAtoms()
                    }

                synthon_atom_map_numbers = ReactionReactivityUtility.get_synthon_atom_map_numbers(
                    mapped_reactant_compound_mol=mapped_reactant_compound_mol,
                    mapped_product_compound_mol=mapped_product_compound_mol,
                    atom_property_keys=atom_property_keys,
                    bond_property_keys=bond_property_keys
                )

                reactant_reactive_site_atom_indices, reactant_synthon_atom_indices = set(), dict()

                for reactant_atom in mapped_reactant_compound_mol.GetAtoms():
                    if reactant_atom.GetAtomMapNum() not in synthon_atom_map_numbers:
                        reactant_reactive_site_atom_indices.add(
                            reactant_atom.GetIdx()
                        )

                for synthon_atom_map_number in synthon_atom_map_numbers:
                    if synthon_atom_map_number in reactant_atom_map_number_to_index_dictionaries[reactant_index].keys():
                        reactant_synthon_atom_indices[
                            reactant_atom_map_number_to_index_dictionaries[reactant_index][synthon_atom_map_number]
                        ] = product_atom_map_number_to_index_dictionaries[reactant_index][synthon_atom_map_number]

                reactant_reactive_sites_and_synthons[reactant_index] = (
                    reactant_reactive_site_atom_indices,
                    reactant_synthon_atom_indices,
                )

            product_reactive_site_atom_indices = set()

            for product_atom in mapped_product_compound_mol.GetAtoms():
                if not any(
                    product_atom.GetIdx() in reactant_synthon_atom_indices.values()
                    for _, reactant_synthon_atom_indices in reactant_reactive_sites_and_synthons.values()
                ):
                    product_reactive_site_atom_indices.add(
                        product_atom.GetIdx()
                    )

            product_reactive_sites_and_synthons[product_index] = (
                reactant_reactive_sites_and_synthons,
                product_reactive_site_atom_indices,
            )

        return product_reactive_sites_and_synthons
