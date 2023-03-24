import logging

import click
import pandas as pd
from rxn.chemutils.conversion import mol_to_smiles
from rxn.utilities.logging import setup_console_logger
from tqdm import tqdm

from rxn_standardization.utils import remove_stereochemistry

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.command()
@click.option(
    "--sid_map_file",
    "-sm",
    type=str,
    help="Path to SID-Map ASCII file. First column is SID, fourth column is CID.",
)
@click.option(
    "--cid_smiles_file",
    "-c",
    type=str,
    help="Path to CID-SMILES ASCII file. First column is CID, second column is SMILES (standardized).",
)
@click.option(
    "--sid_smiles_file",
    "-ss",
    type=str,
    help="Path to SID-SMILES ASCII file. First column is SID, second column is SMILES (non-standardized).",
)
@click.option(
    "--output_file",
    "-o",
    type=str,
    help="Path to output CSV file containing 2 columns: src, tgt.",
)
def main(
    sid_map_file: str,
    cid_smiles_file: str,
    sid_smiles_file: str,
    output_file: str,
):
    """
    Extract src, tgt SMILES from PubChem files. 3 relevant ASCII files are downloaded from https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/ (src)
    and https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/ (tgt):
    SID-Map: maps substance ids to compound ids
    CID-SMILES: maps compound ids to their SMILES
    SID-SMILES: maps substance ids to their SMILES
    """
    setup_console_logger()

    sid_map = pd.concat(
        [
            chunk
            for chunk in tqdm(
                pd.read_csv(
                    sid_map_file,
                    sep="\t",
                    usecols=[0, 3],
                    header=None,
                    names=["sid", "cid"],
                    engine="python",
                    dtype={"sid": int, "cid": "Int32"},
                    chunksize=1000,
                ),
                desc="Loading SID-MAP",
            )
        ]
    )
    sid_map.dropna(inplace=True)
    cid_smiles = pd.concat(
        [
            chunk
            for chunk in tqdm(
                pd.read_csv(
                    cid_smiles_file,
                    sep="\t",
                    header=None,
                    names=["cid", "smiles"],
                    dtype={"cid": int, "smiles": str},
                    chunksize=1000,
                ),
                desc="Loading CID-SMILES",
            )
        ]
    )
    sid_smiles = pd.concat(
        [
            chunk
            for chunk in tqdm(
                pd.read_csv(
                    sid_smiles_file,
                    header=None,
                    names=["sid", "smiles"],
                    dtype={"sid": int, "smiles": str},
                    chunksize=1000,
                ),
                desc="Loading SID-SMILES",
            )
        ]
    )

    sid_map_dict = dict(zip(sid_map.sid, sid_map.cid))
    cid_smiles_dict = dict(zip(cid_smiles.cid, cid_smiles.smiles))
    sid_smiles_dict = dict(zip(sid_smiles.sid, sid_smiles.smiles))

    substance_compound_dict = {
        sid_smiles_dict[key]: cid_smiles_dict[sid_map_dict[key]]
        for key in sid_smiles_dict.keys()
        if key in sid_map_dict.keys()
    }
    substance_compound_df = pd.DataFrame(
        {"src": substance_compound_dict.keys(), "tgt": substance_compound_dict.values()}
    )
    substance_compound_df.dropna(inplace=True)
    
    # Remove stereochemistry
    logger.info("Removing stereochemistry...")
    substance_compound_df["src"] = substance_compound_df["src"].apply(remove_stereochemistry)
    substance_compound_df["tgt"] = substance_compound_df["tgt"].apply(remove_stereochemistry)

    substance_compound_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
