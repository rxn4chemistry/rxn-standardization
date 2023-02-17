import logging
from typing import Optional

import os
import click
import pandas as pd
from rxn.utilities.files import is_path_creatable
from rxn.utilities.logging import setup_console_logger
from rxn.chemutils.tokenization import tokenize_smiles

from rxn_standardization.some_module import process_input

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

@click.command()
@click.option(
    "--input_csv",
    "-i",
    type=str,
    requrired=True,
    help="Path to the input file with src, tgt columns."
)
@click.option(
    "--save_dir",
    "-s",
    type=str,
    requrired=True,
    help="Path to the output directory."
)
@click.option(
    "--prepend_token",
    "-p",
    type=str,
    default=None,
    help="Whether to prepend token to source SMILES and what token that is. Example format: [PUBCHEM]."
)
@click.option(
    "--train_frac",
    "-t",
    type=float,
    default=0.9,
    help="Fraction of dataset to use for training."
)

def main(
    input_csv: str,
    save_dir: str,
    prepend_token: Optional[str],
    train_frac: float,
):
    "Tokenize SMILES, split dataset and generate source and target files."
    setup_console_logger()

    if not is_path_creatable(f"{save_dir}/src-train.txt"):
        raise ValueError(f'Permissions insufficient to create file in "{save_dir}".')

    # Read csv
    df: pd.DataFrame = pd.read_csv(input_csv)
    df.columns = ["src", "tgt"]

    # Tokenize SMILES
    df["src"] = [process_input(tokenize_smiles(smi)) for smi in df["src"].values]
    df["tgt"] = [process_input(tokenize_smiles(smi)) for smi in df["tgt"].values]

    # Prepend token
    if prepend_token is not None:
        df["src"] = [f"{prepend_token} {smi}" for smi in df["src"].values]

    # Split data into train/test/validation sets
    train = df.sample(frac=train_frac, random_state=42)
    test_valid = df.drop(train.index)
    test = test_valid.sample(frac=0.6, random_state=42)
    valid = test.drop(test.index)

    # Save files
    if not save_dir.exists():
        os.mkdir(save_dir)
    train["src"].to_csv(f"{save_dir}/src-train.txt", header=False, index=False)
    train["tgt"].to_csv(f"{save_dir}/tgt-train.txt", header=False, index=False)
    test["src"].to_csv(f"{save_dir}/src-test.txt", header=False, index=False)
    test["tgt"].to_csv(f"{save_dir}/tgt-test.txt", header=False, index=False)
    valid["src"].to_csv(f"{save_dir}/src-valid.txt", header=False, index=False)
    valid["tgt"].to_csv(f"{save_dir}/tgt-valid.txt", header=False, index=False)

if __name__ == "__main__":
    main()
