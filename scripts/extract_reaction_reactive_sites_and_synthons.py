""" The ``scripts`` directory ``extract_reaction_reactive_sites_and_synthons`` script. """

from argparse import ArgumentParser, Namespace
from logging import Formatter, Logger, StreamHandler, getLogger

from ncsw_chemistry.utility.reaction import ReactionFormatConversionUtility, ReactionReactivityUtility
from ncsw_chemistry.utility.reaction.compound import ReactionCompoundExtractionUtility, ReactionCompoundSanitizationUtility


def get_script_arguments() -> Namespace:
    """
    Get the script arguments.

    :returns: The script arguments.
    """

    argument_parser = ArgumentParser()

    argument_parser.add_argument(
        "-mrs",
        "--mapped_reaction_smiles",
        type=str,
        required=False,  # True
        help="The mapped chemical reaction SMILES string."
    )

    return argument_parser.parse_args()


def get_script_logger() -> Logger:
    """
    Get the script logger.

    :returns: The script logger.
    """

    logger = getLogger(
        name=__name__
    )

    logger.setLevel(
        level="INFO"
    )

    formatter = Formatter(
        fmt="[{asctime:s}] {levelname:s}: \"{message:s}\"",
        style="{"
    )

    stream_handler = StreamHandler()

    stream_handler.setLevel(
        level="INFO"
    )

    stream_handler.setFormatter(
        fmt=formatter
    )

    logger.addHandler(
        hdlr=stream_handler
    )

    return logger


if __name__ == "__main__":
    script_logger = get_script_logger()

    try:
        script_arguments = get_script_arguments()

        script_arguments.mapped_reaction_smiles = "Cl[C:1]([CH3:2])=[O:3].[NH:4]=[S:5](=[O:6])([CH2:7][CH2:8][OH:9])[c:10]1[cH:11][cH:12][c:13]([C:14](=[O:15])[NH:16][c:17]2[cH:18][cH:19][c:20]([Cl:21])[cH:22][c:23]2[C:24](=[O:25])[NH:26][c:27]2[cH:28][cH:29][c:30]([Cl:31])[cH:32][n:33]2)[cH:34][cH:35]1>>[C:1]([CH3:2])(=[O:3])[N:4]=[S:5](=[O:6])([CH2:7][CH2:8][OH:9])[c:10]1[cH:11][cH:12][c:13]([C:14](=[O:15])[NH:16][c:17]2[cH:18][cH:19][c:20]([Cl:21])[cH:22][c:23]2[C:24](=[O:25])[NH:26][c:27]2[cH:28][cH:29][c:30]([Cl:31])[cH:32][n:33]2)[cH:34][cH:35]1"

        reaction_rxn = ReactionFormatConversionUtility.convert_smiles_to_rxn(
            reaction_smiles=script_arguments.mapped_reaction_smiles
        )

        ReactionCompoundSanitizationUtility.sanitize_compounds(
            reaction_rxn=reaction_rxn,
            deep_copy=False
        )

        reaction_compounds = ReactionCompoundExtractionUtility.extract_compounds(
            reaction_rxn=reaction_rxn
        )

        print([
                reactant_compound[0]
                for reactant_compound in reaction_compounds[0]
            ])

        reaction_reactive_sites_and_synthons = ReactionReactivityUtility.extract_reactive_sites_and_synthons_using_atom_map_numbers(
            mapped_reactant_compound_mols=[
                reactant_compound[0]
                for reactant_compound in reaction_compounds[0]
            ],
            mapped_product_compound_mols=[
                product_compound[0]
                for product_compound in reaction_compounds[2]
            ],
            atom_property_keys=[
                "atomic_number",
                "degree",
                "formal_charge",
                "hybridization",
                "implicit_valence",
                "is_aromatic",
                "is_in_ring",
                "number_of_explicit_hydrogens",
                "number_of_implicit_hydrogens",
                "number_of_radical_electrons",
                "symbol",
                "total_degree",
                "total_number_of_hydrogens",
                "total_valence",
            ]
        )

        print()
        print(reaction_rxn)
        print()
        print(reaction_compounds)
        print()
        print(reaction_reactive_sites_and_synthons)

    except Exception as exception_handle:
        script_logger.error(
            msg=exception_handle
        )

        raise
