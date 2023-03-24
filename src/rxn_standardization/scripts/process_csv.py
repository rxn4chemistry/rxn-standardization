import logging
import os
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from rxn.chemutils.tokenization import TokenizationError, tokenize_smiles
from rxn.utilities.files import is_path_creatable
from rxn.utilities.logging import setup_console_logger

from rxn_standardization.utils import process_input

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def smiles_to_tokens(smiles: str) -> Optional[str]:
    """
    Tokenizes SMILES and raises a warning in case of TokenizationError.
    """
    try:
        return process_input(tokenize_smiles(smiles))
    except TokenizationError as e:
        logger.warning(
            f"Error during tokenizing {smiles}: {e.title}, {e.detail}. Skipping this entry."
        )
        return None
    except TypeError:
        logger.warning(f"Error during tokenizing {smiles}. Skipping this entry.")
        return None


@click.command(context_settings={"show_default": True})
@click.option(
    "--input_csv",
    "-i",
    type=str,
    required=True,
    help="Path to the input file with non-standardized and standardized SMILES.",
)
@click.option(
    "--src_col",
    "-sc",
    type=str,
    default="src",
    help="Name of column holding source SMILES strings.",
)
@click.option(
    "--tgt_col",
    "-tc",
    type=str,
    default="tgt",
    help="Name of column holding target SMILES strings.",
)
@click.option(
    "--save_dir",
    "-s",
    type=str,
    required=True,
    help="Path to the output directory.",
)
@click.option(
    "--prepend_token",
    "-p",
    type=str,
    default=None,
    help="Whether to prepend token to source SMILES and what token that is. Example format: [PUBCHEM].",
)
@click.option(
    "--train_frac",
    "-t",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.9,
    help="Fraction of dataset to use for training.",
)
def main(
    input_csv: str,
    src_col: str,
    tgt_col: str,
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

    # Tokenize SMILES
    df[src_col] = [smiles_to_tokens(smi) for smi in df[src_col].values]
    df[tgt_col] = [smiles_to_tokens(smi) for smi in df[tgt_col].values]
    df.dropna(inplace=True)  # Drop the rows where smiles_to_tokens returned None

    # Prepend token
    if prepend_token is not None:
        df[src_col] = [f"{prepend_token} {smi}" for smi in df[src_col].values]

    # Split data into train/test/validation sets
    train = df.sample(frac=train_frac, random_state=42)
    test_valid = df.drop(train.index)
    test = test_valid.sample(frac=0.6, random_state=42)
    valid = test_valid.drop(test.index)

    # Save files
    if not Path(save_dir).exists():
        os.mkdir(save_dir)
    train[src_col].to_csv(f"{save_dir}/src-train.txt", header=False, index=False)
    train[tgt_col].to_csv(f"{save_dir}/tgt-train.txt", header=False, index=False)
    test[src_col].to_csv(f"{save_dir}/src-test.txt", header=False, index=False)
    test[tgt_col].to_csv(f"{save_dir}/tgt-test.txt", header=False, index=False)
    valid[src_col].to_csv(f"{save_dir}/src-valid.txt", header=False, index=False)
    valid[tgt_col].to_csv(f"{save_dir}/tgt-valid.txt", header=False, index=False)


if __name__ == "__main__":
    main()
