import logging
from typing import Optional

import click
from rxn.chemutils.tokenization import detokenize_smiles
from rxn.metrics.metrics import top_n_accuracy
from rxn.utilities.containers import chunker
from rxn.utilities.files import load_list_from_file
from rxn.utilities.logging import setup_console_logger

from rxn_standardization.utils import (
    canonicalize,
    get_sequence_multiplier,
    remove_stereochemistry,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.command()
@click.option(
    "--pred_file",
    "-p",
    type=str,
    required=True,
    help="File with predictions, one SMILES string per line.",
)
@click.option(
    "--tgt_file",
    "-t",
    type=str,
    required=True,
    help="File with ground truth outputs, one SMILES string per line.",
)
@click.option(
    "--src_file",
    "-s",
    type=str,
    required=False,
    help="File with source SMILES strings, only provide if you want to calculate accuracy for compounds that get modified only.",
)
@click.option(
    "--remove_stereo/--keep_stereo",
    "-r/-k",
    default=False,
    help="Whether to remove stereochemistry before comparing predictions to targets.",
)
@click.option(
    "--modified_score/--overall_score",
    "-m/-o",
    default=False,
    help="Whether to calculate accuracy for compounds that get modified during standardization only.",
)
@click.option(
    "--canonicalize_pred/--no_canonicalize",
    "-c/-n_c",
    default=True,
    help="Whether to canonicalize predictions and targets before comparing.",
)
def main(
    pred_file: str,
    tgt_file: str,
    src_file: Optional[str],
    remove_stereo: bool,
    modified_score: bool,
    canonicalize_pred: bool,
):
    setup_console_logger()
    predictions = load_list_from_file(pred_file)
    targets = load_list_from_file(tgt_file)

    if modified_score:
        if src_file is None:
            raise ValueError(
                "The source file must be provided to calculate the modified score"
            )
        source = load_list_from_file(src_file)
        multiplier = get_sequence_multiplier(
            targets, predictions
        )  # In case of top-n accuracy with n>1:
        prediction_chuncks = chunker(predictions, chunk_size=multiplier)
        # In case of top-n accuracy with n>1, ensure that all n predictions are selected:
        predictions = [
            p
            for p_chunk, t, s in zip(prediction_chuncks, targets, source)
            if t != s
            for p in p_chunk
        ]
        targets = [t for t, s in zip(targets, source) if t != s]

    if remove_stereo:
        logger.info("Removing stereochemistry...")
        predictions = [remove_stereochemistry(smi) for smi in predictions]
        targets = [remove_stereochemistry(smi) for smi in targets]

    if canonicalize_pred:
        logger.info("Canonicalizing SMILES...")
        predictions = [canonicalize(detokenize_smiles(smi)) for smi in predictions]
        targets = [canonicalize(detokenize_smiles(smi)) for smi in targets]

    print(top_n_accuracy(targets, predictions))


if __name__ == "__main__":
    main()
