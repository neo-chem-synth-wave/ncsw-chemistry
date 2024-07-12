""" The ``scripts`` directory ``extract_reaction_retro_template`` script. """

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
        "-mrcs",
        "--mapped_reactant_compounds_smiles",
        type=str,
        required=True,
        help="The mapped chemical reaction reactant compounds SMILES string."
    )

    argument_parser.add_argument(
        "-mpcs",
        "--mapped_product_compound_smiles",
        type=str,
        required=True,
        help="The mapped chemical reaction product compound SMILES string."
    )

    return argument_parser.parse_args()


def get_script_logger() -> Logger:
    """
    Get the script logger.

    :returns: The script logger.
    """

    logger = getLogger(
        name="extract_reaction_retro_template"
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
            "mapped_reaction_smiles": ">>".join([
                script_arguments.mapped_reactant_compounds_smiles,
                script_arguments.mapped_product_compound_smiles,
            ]),
            "reaction_retro_template_smarts": ReactionTemplateUtility.extract_retro_template_using_rdchiral(
                mapped_reactant_compounds_smiles=script_arguments.mapped_reactant_compounds_smiles,
                mapped_product_compound_smiles=script_arguments.mapped_product_compound_smiles
            )
        })

    except Exception as exception_handle:
        script_logger.error(
            msg=exception_handle
        )

        raise
