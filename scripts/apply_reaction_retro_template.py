""" The ``scripts`` directory ``apply_reaction_retro_template`` script. """

from argparse import ArgumentParser, Namespace
from logging import Formatter, Logger, StreamHandler, getLogger

from ncsw_chemistry.utility.reaction import ReactionTemplateUtility


def get_script_arguments() -> Namespace:
    """
    Get the script arguments.

    :returns: The script arguments.
    """

    argument_parser = ArgumentParser()

    argument_parser.add_argument(
        "-rrts",
        "--reaction_retro_template_smarts",
        type=str,
        required=True,
        help="The chemical reaction retro template SMARTS string."
    )

    argument_parser.add_argument(
        "-cs",
        "--compound_smiles",
        type=str,
        required=True,
        help="The chemical compound SMILES string."
    )

    return argument_parser.parse_args()


def get_script_logger() -> Logger:
    """
    Get the script logger.

    :returns: The script logger.
    """

    logger = getLogger(
        name="apply_reaction_retro_template"
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

        print({
            "reaction_retro_template_smarts": script_arguments.reaction_retro_template_smarts,
            "compound_smiles": script_arguments.compound_smiles,
        })

        print({
            "precursor_compounds_smiles_strings": ReactionTemplateUtility.apply_retro_template_using_rdchiral(
                retro_template_smarts=script_arguments.reaction_retro_template_smarts,
                compound_smiles=script_arguments.compound_smiles
            )
        })

    except Exception as exception_handle:
        script_logger.error(
            msg=exception_handle
        )

        raise
