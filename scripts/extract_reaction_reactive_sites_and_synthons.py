""" The ``scripts`` directory ``extract_reaction_reactive_sites_and_synthons`` script. """

from argparse import ArgumentParser, Namespace
from logging import Formatter, Logger, StreamHandler, getLogger

from ncsw_chemistry.utility.reaction import ReactionCompoundExtractionUtility, ReactionReactivityUtility


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
        required=False,
        help="The mapped chemical reaction SMILES string."
    )

    return argument_parser.parse_args()


def get_script_logger() -> Logger:
    """
    Get the script logger.

    :returns: The script logger.
    """

    logger = getLogger(
        name="extract_reaction_reactive_sites_and_synthons"
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

        reaction_compounds = ReactionCompoundExtractionUtility.extract_compounds_from_reaction_smiles(
            reaction_smiles=script_arguments.mapped_reaction_smiles
        )

        print({
            "mapped_reaction_smiles": script_arguments.mapped_reaction_smiles,
        })

        print({
            "reaction_reactive_sites_and_synthons": (
                ReactionReactivityUtility.extract_reactive_sites_and_synthons_using_atom_map_numbers(
                    mapped_reactant_compound_mols=[
                        reactant_compound[1]
                        for reactant_compound in reaction_compounds[0]
                    ],
                    mapped_product_compound_mols=[
                        product_compound[1]
                        for product_compound in reaction_compounds[2]
                    ]
                )
            )
        })

    except Exception as exception_handle:
        script_logger.error(
            msg=exception_handle
        )

        raise
