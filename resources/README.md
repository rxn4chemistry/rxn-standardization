# Resources

This folder contains links to the data used in the manuscript, as well as scripts to process the files and extract source and target SMILES strings. Ultimately, these scripts create CSV files with two columns: non-standardized (`src`) and standardized (`tgt`) SMILES. These are then used by the `rxn-std-process-csv` script.

## Catalysts

The catalyst dataset can be found in the `resources/data` folder as JSON files. To process the file, run:
```bash
python extract_catalysts.py --input_file <file_path> --output_file <file_path>
```

## Tautomers

The Tautobase dataset was downloaded from [here](https://acs.figshare.com/articles/dataset/_i_Tautobase_i_An_Open_Tautomer_Database/11768304). To process the file, run:
```bash
python extract_tautomers_water.py --input_file <file_path> --output_file <file_path>
```

## PubChem

Data was downloaded from the [Substance](https://ftp.ncbi.nlm.nih.gov/pubchem/Substance/) (for `src` SMILES) and [Compound](https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/) (for `tgt` SMILES) PubChem repositories. SID-Map ASCII files map substance IDs to compound ids, CID-SMILES maps compound IDs to their SMILES, Substance SDF files contain information about substances, including their IDs and SMILES.

First, the Substance SDF files must be processed to extract substance IDs and SMILES in a CSV file:
```bash
python smiles_from_sdf.py --input_file <file_path> --output_file <file_path>
```

To extract 'src' and 'tgt' SMILES strings in a CSV file, run:
```bash
python extract_pubchem.py --sid_map_file <file_path> --cid_smiles_file <file_path> --sid_smiles_file <file_path> --output_file <file_path>
```
## ChEMBL

We generated target SMILES strings for the ChEMBL protocol using the [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline/tree/87afedd453e388cb4759ed03259169fa6c324415).

# Tanimoto split

To make train/test splits based on Tanimoto indices, we used the pipeline developed by [Kovacs et al](https://github.com/davkovacs/MTExplainer/tree/master/data/tanimoto_splits).
