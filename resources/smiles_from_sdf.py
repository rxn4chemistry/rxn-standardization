import logging

import click
import pandas as pd
from rdkit.Chem import PandasTools
from rxn.utilities.logging import setup_console_logger

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.command()
@click.option(
    "--input_file",
    type=str,
    required=True,
    help="Path to input SDF file.",
)
@click.option(
    "--output_file",
    type=str,
    required=True,
    help="Path to output CSV file.",
)
def main(
    input_file: str,
    output_file: str,
) -> None:
    """
    Retrieve substance IDs and SMILES from PubChem Substance SDF files.
    """
    setup_console_logger()

    df = PandasTools.LoadSDF(
        input_file, molColName=None, smilesName="SMILES", embedProps=True
    )

    new_df = pd.DataFrame({"SID": df.PUBCHEM_SUBSTANCE_ID, "smiles": df.SMILES})
    new_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
