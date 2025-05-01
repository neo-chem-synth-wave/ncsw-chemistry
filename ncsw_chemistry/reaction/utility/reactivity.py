""" The ``ncsw_chemistry.reaction.utility`` package ``reactivity`` module. """

from itertools import chain
from typing import Dict, List, Optional, Sequence, Set, Tuple

from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun
from rdchiral.template_extractor import extract_from_reaction

from rdkit.Chem.rdchem import Mol

from ncsw_chemistry.compound.utility.atom import CompoundAtomUtility
from ncsw_chemistry.compound.utility.bond import CompoundBondUtility


class ReactionReactivityUtility:
    """ The chemical reaction reactivity utility class. """

    @staticmethod
    def extract_retro_template_using_rdchiral(
            mapped_reactant_compound_smiles_strings: Sequence[str],
            mapped_product_compound_smiles: str
    ) -> Optional[str]:
        """
        Extract the retro template from a chemical reaction using the RDChiral library.

        :parameter mapped_reactant_compound_smiles_strings: The SMILES strings of the mapped chemical reaction reactant
            compounds.
        :parameter mapped_product_compound_smiles: The SMILES string of the mapped chemical reaction product compound.

        :returns: The chemical reaction retro template.
        """

        return extract_from_reaction({
            "_id": None,
            "reactants": ".".join(mapped_reactant_compound_smiles_strings),
            "products": mapped_product_compound_smiles,
        }).get("reaction_smarts", None)

    @staticmethod
    def apply_retro_template_using_rdchiral(
            retro_template_smarts: str,
            compound_smiles: str,
            **kwargs
    ) -> Optional[List[str]]:
        """
        Apply a chemical reaction retro template on a chemical compound using the RDChiral library.

        :parameter retro_template_smarts: The chemical reaction retro template SMARTS string.
        :parameter compound_smiles: The SMILES string of the chemical compound.
        :parameter kwargs: The keyword arguments for the adjustment of the following underlying functions:
            { `rdchiral.main.rdchiralRun` }.

        :returns: The outcomes of the application of the chemical reaction retro template on the chemical compound.
        """

        return rdchiralRun(
            rxn=rdchiralReaction(
                reaction_smarts=retro_template_smarts
            ),
            reactants=rdchiralReactants(
                reactant_smiles=compound_smiles
            ),
            **kwargs
        )

    @staticmethod
    def get_synthon_atom_map_numbers(
            mapped_reactant_compound_mol: Mol,
            mapped_product_compound_mol: Mol,
            atom_property_keys: Optional[Sequence[str]] = None,
            bond_atom_property_keys: Optional[Sequence[str]] = None,
            bond_property_keys: Optional[Sequence[str]] = None
    ) -> Set[int]:
        """
        Get the synthon atom map numbers of the mapped chemical reaction reactant and product compounds.

        :parameter mapped_reactant_compound_mol: The RDKit Mol object of the mapped chemical reaction reactant compound.
        :parameter mapped_product_compound_mol: The RDKit Mol object of mapped chemical reaction product compound.
        :parameter atom_property_keys: The keys of the chemical reaction compound atom properties that should be
            utilized in the property ID. The value `None` indicates that all chemical reaction compound atom properties
            should be utilized in the property ID.
        :parameter bond_atom_property_keys: The keys of the chemical reaction compound bond atom properties that should
            be utilized in the property ID. The value `None` indicates that all chemical reaction compound bond atom
            properties should be utilized in the property ID.
        :parameter bond_property_keys: The keys of the chemical reaction compound bond properties that should be
            utilized in the property ID. The value `None` indicates that all chemical reaction compound bond properties
            should be utilized in the property ID.

        :returns: The synthon atom map numbers of the mapped chemical reaction reactant and product compounds.
        """

        reactant_compound_bond_atom_map_numbers = set()

        for reactant_compound_bond in mapped_reactant_compound_mol.GetBonds():
            if reactant_compound_bond.GetBeginAtom().HasProp(
                key="molAtomMapNumber"
            ) and reactant_compound_bond.GetEndAtom().HasProp(
                key="molAtomMapNumber"
            ):
                reactant_compound_bond_atom_map_numbers.add(
                    frozenset({
                        reactant_compound_bond.GetBeginAtom().GetAtomMapNum(),
                        reactant_compound_bond.GetEndAtom().GetAtomMapNum(),
                    })
                )

        product_compound_bond_atom_map_numbers = set()

        for product_compound_bond in mapped_product_compound_mol.GetBonds():
            if product_compound_bond.GetBeginAtom().HasProp(
                key="molAtomMapNumber"
            ) and product_compound_bond.GetEndAtom().HasProp(
                key="molAtomMapNumber"
            ):
                product_compound_bond_atom_map_numbers.add(
                    frozenset({
                        product_compound_bond.GetBeginAtom().GetAtomMapNum(),
                        product_compound_bond.GetEndAtom().GetAtomMapNum(),
                    })
                )

        non_synthon_bond_atom_map_numbers = set(
            chain.from_iterable(
                reactant_compound_bond_atom_map_numbers.symmetric_difference(
                    product_compound_bond_atom_map_numbers
                )
            )
        )

        synthon_bond_atom_map_numbers = set()

        for reactant_compound_bond in mapped_reactant_compound_mol.GetBonds():
            if reactant_compound_bond.GetBeginAtom().HasProp(
                key="molAtomMapNumber"
            ) and reactant_compound_bond.GetEndAtom().HasProp(
                key="molAtomMapNumber"
            ):
                for product_compound_bond in mapped_product_compound_mol.GetBonds():
                    if product_compound_bond.GetBeginAtom().HasProp(
                        key="molAtomMapNumber"
                    ) and product_compound_bond.GetEndAtom().HasProp(
                        key="molAtomMapNumber"
                    ):
                        if {
                            reactant_compound_bond.GetBeginAtom().GetAtomMapNum(),
                            reactant_compound_bond.GetEndAtom().GetAtomMapNum(),
                        } == {
                            product_compound_bond.GetBeginAtom().GetAtomMapNum(),
                            product_compound_bond.GetEndAtom().GetAtomMapNum(),
                        }:
                            if CompoundBondUtility.get_bond_property_id(
                                bond=reactant_compound_bond,
                                bond_atom_property_keys=bond_atom_property_keys,
                                bond_property_keys=bond_property_keys
                            ) == CompoundBondUtility.get_bond_property_id(
                                bond=product_compound_bond,
                                bond_atom_property_keys=bond_atom_property_keys,
                                bond_property_keys=bond_property_keys
                            ):
                                synthon_bond_atom_map_numbers.add(
                                    reactant_compound_bond.GetBeginAtom().GetAtomMapNum()
                                    # product_compound_bond.GetBeginAtom().GetAtomMapNum()
                                )

                                synthon_bond_atom_map_numbers.add(
                                    reactant_compound_bond.GetEndAtom().GetAtomMapNum()
                                    # product_compound_bond.GetEndAtom().GetAtomMapNum()
                                )

                            else:
                                non_synthon_bond_atom_map_numbers.add(
                                    reactant_compound_bond.GetBeginAtom().GetAtomMapNum()
                                    # product_compound_bond.GetBeginAtom().GetAtomMapNum()
                                )

                                non_synthon_bond_atom_map_numbers.add(
                                    reactant_compound_bond.GetEndAtom().GetAtomMapNum()
                                    # product_compound_bond.GetEndAtom().GetAtomMapNum()
                                )

        synthon_atom_map_numbers, non_synthon_atom_map_numbers = set(), set()

        for reactant_compound_atom in mapped_reactant_compound_mol.GetAtoms():
            if (
                reactant_compound_atom.HasProp(
                    key="molAtomMapNumber"
                ) and reactant_compound_atom.GetAtomMapNum() not in synthon_bond_atom_map_numbers and
                reactant_compound_atom.GetAtomMapNum() not in non_synthon_bond_atom_map_numbers
            ):
                for product_compound_atom in mapped_product_compound_mol.GetAtoms():
                    if (
                        product_compound_atom.HasProp(
                            key="molAtomMapNumber"
                        ) and product_compound_atom.GetAtomMapNum() not in synthon_bond_atom_map_numbers and
                        product_compound_atom.GetAtomMapNum() not in non_synthon_bond_atom_map_numbers
                    ):
                        if reactant_compound_atom.GetAtomMapNum() == product_compound_atom.GetAtomMapNum():
                            if CompoundAtomUtility.get_atom_property_id(
                                atom=reactant_compound_atom,
                                atom_property_keys=atom_property_keys
                            ) == CompoundAtomUtility.get_atom_property_id(
                                atom=product_compound_atom,
                                atom_property_keys=atom_property_keys
                            ):
                                synthon_atom_map_numbers.add(
                                    reactant_compound_atom.GetAtomMapNum()
                                    # product_compound_atom.GetAtomMapNum()
                                )

                            else:
                                non_synthon_atom_map_numbers.add(
                                    reactant_compound_atom.GetAtomMapNum()
                                    # product_compound_atom.GetAtomMapNum()
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
    def extract_reactive_sites_and_synthons(
            mapped_reactant_compound_mols: Sequence[Mol],
            mapped_product_compound_mols: Sequence[Mol],
            atom_property_keys: Optional[Sequence[str]] = None,
            bond_atom_property_keys: Optional[Sequence[str]] = None,
            bond_property_keys: Optional[Sequence[str]] = None
    ) -> Dict[int, Tuple[Dict[int, Tuple[Set[int], Dict[int, int]]], Set[int]]]:
        """
        Extract the reactive sites and synthons of the chemical reaction reactant and product compounds.

        :parameter mapped_reactant_compound_mols: The RDKit Mol objects of the mapped chemical reaction reactant
            compounds.
        :parameter mapped_product_compound_mols: The RDKit Mol objects of the mapped chemical reaction product
            compounds.
        :parameter atom_property_keys: The keys of the chemical reaction compound atom properties that should be
            utilized in the property ID. The value `None` indicates that all chemical reaction compound atom properties
            should be utilized in the property ID.
        :parameter bond_atom_property_keys: The keys of the chemical reaction compound bond atom properties that should
            be utilized in the property ID. The value `None` indicates that all chemical reaction compound bond atom
            properties should be utilized in the property ID.
        :parameter bond_property_keys: The keys of the chemical reaction compound bond properties that should be
            utilized in the property ID. The value `None` indicates that all chemical reaction compound bond properties
            should be utilized in the property ID.

        :returns: The reactive sites and synthons of the chemical reaction reactant and product compounds.
        """

        product_compound_atom_map_number_to_index_dictionaries = dict()
        reactant_compound_atom_map_number_to_index_dictionaries = dict()

        product_compound_reactive_sites_and_synthons = dict()

        for product_compound_index, product_compound_mol in enumerate(mapped_product_compound_mols):
            if product_compound_index not in product_compound_atom_map_number_to_index_dictionaries.keys():
                product_compound_atom_map_number_to_index_dictionaries[product_compound_index] = {
                    product_compound_atom.GetAtomMapNum(): product_compound_atom.GetIdx()
                    for product_compound_atom in product_compound_mol.GetAtoms()
                }

            reactant_compound_reactive_sites_and_synthons = dict()

            for reactant_compound_index, reactant_compound_mol in enumerate(mapped_reactant_compound_mols):
                if reactant_compound_index not in reactant_compound_atom_map_number_to_index_dictionaries.keys():
                    reactant_compound_atom_map_number_to_index_dictionaries[reactant_compound_index] = {
                        reactant_compound_atom.GetAtomMapNum(): reactant_compound_atom.GetIdx()
                        for reactant_compound_atom in reactant_compound_mol.GetAtoms()
                    }

                synthon_atom_map_numbers = ReactionReactivityUtility.get_synthon_atom_map_numbers(
                    mapped_reactant_compound_mol=reactant_compound_mol,
                    mapped_product_compound_mol=product_compound_mol,
                    atom_property_keys=atom_property_keys,
                    bond_atom_property_keys=bond_atom_property_keys,
                    bond_property_keys=bond_property_keys
                )

                reactant_compound_reactive_site_atom_indices, reactant_compound_synthon_atom_indices = set(), dict()

                for reactant_compound_atom in reactant_compound_mol.GetAtoms():
                    if reactant_compound_atom.GetAtomMapNum() not in synthon_atom_map_numbers:
                        reactant_compound_reactive_site_atom_indices.add(
                            reactant_compound_atom.GetIdx()
                        )

                for synthon_atom_map_number in synthon_atom_map_numbers:
                    if synthon_atom_map_number in reactant_compound_atom_map_number_to_index_dictionaries[
                        reactant_compound_index
                    ].keys():
                        reactant_compound_synthon_atom_indices[
                            reactant_compound_atom_map_number_to_index_dictionaries[
                                reactant_compound_index
                            ][synthon_atom_map_number]
                        ] = product_compound_atom_map_number_to_index_dictionaries[
                            product_compound_index
                        ][synthon_atom_map_number]

                reactant_compound_reactive_sites_and_synthons[reactant_compound_index] = (
                    reactant_compound_reactive_site_atom_indices,
                    reactant_compound_synthon_atom_indices,
                )

            product_compound_reactive_site_atom_indices = set()

            for product_compound_atom in product_compound_mol.GetAtoms():
                if not any(
                    product_compound_atom.GetIdx() in reactant_compound_synthon_atom_indices.values()
                    for _, reactant_compound_synthon_atom_indices in
                    reactant_compound_reactive_sites_and_synthons.values()
                ):
                    product_compound_reactive_site_atom_indices.add(
                        product_compound_atom.GetIdx()
                    )

            product_compound_reactive_sites_and_synthons[product_compound_index] = (
                reactant_compound_reactive_sites_and_synthons,
                product_compound_reactive_site_atom_indices,
            )

        return product_compound_reactive_sites_and_synthons
