""" The ``ncsw_chemistry.utility.reaction`` package ``template`` module. """

from typing import Any, Dict, List, Optional, Sequence, Tuple

from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun
from rdchiral.template_extractor import extract_from_reaction

from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import Mol


class ReactionTemplateUtility:
    """ The chemical reaction template utility class. """

    @staticmethod
    def extract_retro_template_using_rdchiral(
            mapped_reactant_compound_smiles_strings: Sequence[str],
            mapped_product_compound_smiles: str
    ) -> Optional[Dict[str, Any]]:
        """
        Extract the chemical reaction retro template using the `RDChiral` library.

        :parameter mapped_reactant_compound_smiles_strings: The mapped chemical reaction reactant compound `SMILES`
            strings.
        :parameter mapped_product_compound_smiles: The mapped chemical reaction product compound SMILES string.

        :returns: The chemical reaction retro template.
        """

        return extract_from_reaction({
            "_id": None,
            "reactants": ".".join(mapped_reactant_compound_smiles_strings),
            "products": mapped_product_compound_smiles,
        })

    @staticmethod
    def apply_retro_template_using_rdchiral(
            retro_template_smarts: str,
            compound_smiles: str,
            **kwargs
    ) -> Optional[List[str]]:
        """
        Apply a chemical reaction retro template on a chemical compound using the `RDChiral` library.

        :parameter retro_template_smarts: The chemical reaction retro template `SMARTS` string.
        :parameter compound_smiles: The chemical compound `SMILES` string.
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
    def apply_template_using_rdkit(
            reaction_template_rxn: ChemicalReaction,
            compound_mols: Tuple[Mol],
            maximum_number_of_outcomes: int = 1000
    ) -> Tuple[Tuple[Mol, ...]]:
        """
        Apply a chemical reaction template on the chemical compounds using the `RDKit` library.

        :parameter reaction_template_rxn: The chemical reaction template `ChemicalReaction` object.
        :parameter compound_mols: The chemical compound `RDKit Mol` objects.
        :parameter maximum_number_of_outcomes: The maximum number of outcomes.

        :returns: The outcomes of the application of the chemical reaction template on the chemical compounds.
        """

        if not reaction_template_rxn.IsInitialized():
            reaction_template_rxn.Initialize()

        return reaction_template_rxn.RunReactants(
            reactants=compound_mols,
            maxProducts=maximum_number_of_outcomes
        )
