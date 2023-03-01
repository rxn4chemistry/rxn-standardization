import logging
from typing import Any, Optional, Tuple

import click
import pandas as pd
from rdkit import Chem
from rxn.utilities.logging import setup_console_logger

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def remove_atom_mapping(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def first_tautomer_is_src(
    log_K: Any, tautomer_1_percent: str, preferred: str
) -> Optional[bool]:
    """Determine, from values given in Tautobase, which tautomer is the source
    and which is the target.

    Returns:
        True if the first tautomer should be the source, False if the second
        tautomer should be the source, None if this cannot be determined.
    """
    if log_K != 0:  # if log_K is given
        try:
            return float(log_K) < 0
        except ValueError:
            try:
                # if log_K is of format ">4"
                return float(log_K[1:]) < 0
            except ValueError:
                # if log_K is of format "<<0"
                return float(log_K[2:]) < 0
    elif tautomer_1_percent != 0:  # if percentage of tautomer 1 is given
        try:
            return float(tautomer_1_percent) >= 50
        except ValueError:
            return float(tautomer_1_percent[1:]) >= 50
    elif preferred != 0:  # if the preferred tautomer is specified
        if preferred == "1":
            return True
        elif preferred == "2":
            return False
        elif preferred == "Both":
            return None

    return None


def process_tautobase_entry(
    smirks: str, log_K: Any, tautomer_1_percent: str, preferred: str
) -> Optional[Tuple[str, str]]:
    """
    From a tautobase entry, get the source and target tautomers.

    Returns:
        Tuple: source tautomer, target tautomer. None if this cannot be determined.
    """
    tautomer_1, tautomer_2 = smirks.split(">>")

    try:
        tautomer_1 = remove_atom_mapping(tautomer_1)
        tautomer_2 = remove_atom_mapping(tautomer_2)
    except AttributeError:
        logger.warning(
            f"Error during conversion of RDKit Mol; cannot remove atom mapping from {smirks}. Skipping this entry."
        )
        return None

    first_is_src = first_tautomer_is_src(
        log_K=log_K, tautomer_1_percent=tautomer_1_percent, preferred=preferred
    )
    if first_is_src is None:
        return None
    if first_is_src:
        return tautomer_1, tautomer_2
    else:
        return tautomer_2, tautomer_1


@click.command()
@click.option(
    "--input_file",
    "-i",
    type=str,
    help="Path to TXT file containing tautomers SMIRKS and information on their ratios in solution.",
)
@click.option(
    "--output_file",
    "-o",
    type=str,
    help="Path to output CSV file containing 2 columns: src, tgt.",
)
@click.option(
    "--seed",
    type=int,
    default=42,
    help="Random seed for shuffling the resulting file",
)
def main(
    input_file: str,
    output_file: str,
    seed: float,
) -> None:
    """
    Extract SMILES strings from TXT file containing tautomers as SMIRKS and various columns with info on their ratios in solution.
    Not all columns are populated for each compound and major tautomer is either determined by log_K, percentage of tautomer 1,
    or denoted preferred tautomer. A "solvent" column specifies the solvent for which the major tautomer was determined.
    """
    setup_console_logger()

    df = pd.read_csv(
        input_file,
        delimiter="\t",
        header=0,
        usecols=[0, 2, 3, 4, 5],
        names=["smirks", "log_K", "percent_tautomer_1", "preferred", "solvent"],
    )
    df.fillna(0, inplace=True)  # replace NaN values with None
    print(df)
    src_list = []
    tgt_list = []

    for _, row in df.iterrows():
        smirks, log_K, tautomer_1_percent, preferred, solvent = (
            row["smirks"],
            row["log_K"],
            row["percent_tautomer_1"],
            row["preferred"],
            row["solvent"],
        )

        if solvent != "Water":
            continue

        processed = process_tautobase_entry(
            smirks=smirks,
            log_K=log_K,
            tautomer_1_percent=tautomer_1_percent,
            preferred=preferred,
        )
        if processed is None:
            continue

        src_list.append(processed[0])
        tgt_list.append(processed[1])

    # Create and save dataframe with src, tgt values
    new_df = pd.DataFrame({"src": src_list, "tgt": tgt_list})
    new_df.drop_duplicates(inplace=True)
    new_df = new_df.sample(frac=1, random_state=seed)
    new_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
