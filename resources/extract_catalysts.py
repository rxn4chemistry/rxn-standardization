import json
import logging
from typing import List

import click
import pandas as pd
from rxn.utilities.logging import setup_console_logger

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.command()
@click.option(
    "--input_file",
    "-i",
    type=str,
    required=True,
    help="Path to the JSON data file with standardized compounds.",
)
@click.option(
    "--output_file",
    "-o",
    type=str,
    required=True,
    help="Path to output CSV file containing 2 columns: src, tgt.",
)
def main(
    input_file: str,
    output_file: str,
) -> None:
    """
    Extract SMILES strings from JSON file to CSV file with two columns (src, tgt).
    """
    setup_console_logger()

    src_list: List[str] = []
    tgt_list: List[str] = []

    logger.info(f'Reading file "{input_file}"')
    with open(input_file, "r") as f:
        data = json.load(f)

    # using get() to avoid KeyError if key is not in dictionary
    for compound_dict in data:
        if compound_dict.get("decision") == "accept":
            src_smiles = compound_dict["original_smiles"]
            tgt_smiles = compound_dict["updated_smiles"]
            if tgt_smiles is None:
                tgt_smiles = src_smiles
            src_list.append(src_smiles)
            tgt_list.append(tgt_smiles)

    df = pd.DataFrame({"src": src_list, "tgt": tgt_list})

    df.to_csv(output_file, index=False, header=["src", "tgt"])


if __name__ == "__main__":
    main()
