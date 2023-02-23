import logging

import click
from rxn.chemutils.tokenization import detokenize_smiles
from rxn.utilities.files import (
    dump_list_to_file,
    is_path_creatable,
    iterate_lines_from_file,
)
from rxn.utilities.logging import setup_console_logger

from rxn_standardization.utils import canonicalize

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.command()
@click.option(
    "--input_file",
    "-i",
    type=str,
    required=True,
    help="Path to the input file to detokenize, with one SMILES string per row.",
)
@click.option(
    "--output_file",
    "-o",
    type=str,
    required=True,
    help="Path to the output file.",
)
@click.option(
    "--canonicalize_output/--no_canonicalize",
    "-c/-n_c",
    default=False,
    help="Whether to canonicalize the SMILES.",
)
def main(
    input_file: str,
    output_file: str,
    canonicalize_output: bool,
):
    "Detokenize SMILES."
    setup_console_logger()

    if not is_path_creatable(f"{output_file}"):
        raise ValueError(f'Permissions insufficient to create file "{output_file}".')

    # Read txt
    tokenized_smiles = iterate_lines_from_file(input_file)

    # Detokenize SMILES
    detokenized_smiles = (detokenize_smiles(smi) for smi in tokenized_smiles)

    # Prepend token
    if canonicalize_output:
        logger.info("Canonicalizing SMILES...")
        detokenized_smiles = (canonicalize(smi) for smi in detokenized_smiles)

    dump_list_to_file(detokenized_smiles, output_file)


if __name__ == "__main__":
    main()
