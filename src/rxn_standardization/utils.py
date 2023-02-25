import logging
import re
from typing import Sequence, TypeVar

from rdkit.Chem import RemoveStereochemistry
from rxn.chemutils.conversion import canonicalize_smiles, mol_to_smiles, smiles_to_mol
from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.tokenization import detokenize_smiles
from rxn.utilities.misc import get_multiplier
from rxn.utilities.regex import capturing, optional

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

charge_block = optional("[+-][0-9]*", capture_group=True)
element = "[A-Z][a-z]?"
full_regex = r"\[" + capturing(element) + charge_block + r"\]"

T = TypeVar("T")


def process_token(token: str) -> str:
    """
    Adds a space before and after elements captured by the regex.
    """
    match = re.match(full_regex, token)
    if match is None:
        return token
    captured_new_tokens = [t for t in match.groups() if t is not None]
    return " ".join(["[", *captured_new_tokens, "]"])


def process_input(tokenized_smiles: str) -> str:
    """
    Takes as input tokenized SMILES and returns tokenized SMILES with added spaces in between atoms and their charges,
    as well as before and after  (to emphasise the element charge).
    """
    tokens = tokenized_smiles.split(" ")
    return " ".join(process_token(token) for token in tokens)


def remove_stereochemistry(
    smi: str,
) -> str:
    """
    Remove stereochemistry from SMILES.
    """
    try:
        mol = smiles_to_mol(detokenize_smiles(smi), sanitize=True)
    except InvalidSmiles as e:
        logger.warning(
            f'Invalid SMILES "{e.smiles}"; cannot remove stereochemistry and leaving as is.'
        )
        return smi
    RemoveStereochemistry(mol)
    return mol_to_smiles(mol)


def canonicalize(
    smi: str,
) -> str:
    """Canonicalize SMILES and raise warning if SMILES is invalid."""
    try:
        can_smi = canonicalize_smiles(smi)
    except InvalidSmiles as e:
        logger.warning(
            f'Invalid SMILES "{e.smiles}"; cannot canonicalize and leaving as is.'
        )
        return smi
    return can_smi


def get_sequence_multiplier(ground_truth: Sequence[T], predictions: Sequence[T]) -> int:
    """
    Get the multiplier for the number of predictions by ground truth sample.
    Raises:
        ValueError: if the lists have inadequate sizes (possibly forwarded
            from get_multiplier).
    """
    n_gt = len(ground_truth)
    n_pred = len(predictions)

    return get_multiplier(n_gt, n_pred)
