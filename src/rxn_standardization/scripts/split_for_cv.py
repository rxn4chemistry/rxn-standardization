import logging
import os
import random
from pathlib import Path
from typing import Optional

import click
import pandas as pd
from rxn.chemutils.tokenization import TokenizationError, tokenize_smiles
from rxn.utilities.files import dump_list_to_file, load_list_from_file
from rxn.utilities.logging import setup_console_logger
from sklearn.model_selection import KFold

from rxn_standardization.utils import augment, process_input

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


@click.command(context_settings={"show_default": True})
@click.option(
    "--input_csv",
    "-i",
    type=str,
    required=True,
    help="Path to the input file with non-standardized and standardized SMILES.",
)
@click.option(
    "--save_dir",
    "-s",
    type=str,
    required=True,
    help="Path to the directory where the data will be saved.",
)
@click.option(
    "--prepend_token",
    "-p",
    type=str,
    default=None,
    help="Whether to prepend token to source SMILES and what token that is. Example format: [PUBCHEM].",
)
@click.option(
    "--augmentation",
    "-a",
    type=bool,
    default=True,
    help="Whether to augment SMILES (for training set).",
)
@click.option(
    "--test_size",
    "-t",
    type=int,
    required=True,
    help="Size of held-out test set.",
)

def main(
    input_csv: str,
    save_dir: str,
    test_size: int,
    prepend_token: Optional[str],
    augmentation: bool,
):
    setup_console_logger()

    all_smiles = load_list_from_file(input_csv)
    random.shuffle(all_smiles)
    test_set = all_smiles[:test_size]
    src_test = [smiles_to_tokens(s.split(",")[0]) for s in test_set]
    tgt_test = [smiles_to_tokens(s.split(",")[1]) for s in test_set]

    train_valid_sets = all_smiles[test_size:]

    kf = KFold(n_splits=5)
    for i, (train_index, test_index) in enumerate(kf.split(train_valid_sets)):
        train_set = [train_valid_sets[k] for k in train_index]
        valid_set = [train_valid_sets[k] for k in test_index]

        # Separate into src, tgt and tokenize
        src_train = [smiles_to_tokens(s.split(",")[0]) for s in train_set]
        tgt_train = [smiles_to_tokens(s.split(",")[1]) for s in train_set]

        # Augment SMILES
        if augmentation:
            src_train = augment(src_train, detokenize=True)
            tgt_train = [val for val in tgt_train for _ in (0, 1, 2)]
            # Shuffle
            df = pd.DataFrame({"src": src_train, "tgt": tgt_train})
            df = df.sample(frac=1)
            src_train = df["src"].to_list()
            tgt_train = df["tgt"].to_list()

        src_valid = [smiles_to_tokens(s.split(",")[0]) for s in valid_set]
        tgt_valid = [smiles_to_tokens(s.split(",")[1]) for s in valid_set]

        # Prepend token
        if prepend_token is not None:
            src_train = [f"{prepend_token} {smi}" for smi in src_train]
            src_valid = [f"{prepend_token} {smi}" for smi in src_valid]
            src_test = [f"{prepend_token} {smi}" for smi in src_test]

        # Save files
        if not Path(f"{save_dir}-{i}").exists():
            os.mkdir(f"{save_dir}-{i}")
        dump_list_to_file(src_train, f"{save_dir}-{i}/src-train.txt")
        dump_list_to_file(tgt_train, f"{save_dir}-{i}/tgt-train.txt")
        dump_list_to_file(src_valid, f"{save_dir}-{i}/src-valid.txt")
        dump_list_to_file(tgt_valid, f"{save_dir}-{i}/tgt-valid.txt")
        dump_list_to_file(src_test, f"{save_dir}-{i}/src-test.txt")
        dump_list_to_file(tgt_test, f"{save_dir}-{i}/tgt-test.txt")


if __name__ == "__main__":
    main()