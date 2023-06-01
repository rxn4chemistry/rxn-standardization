# Chemical compounds standardization with language models

This repository contains the code for ["Standardizing chemical compounds with language models"](https://chemrxiv.org/engage/chemrxiv/article-details/6409e14fcc600523a3eb545a).

# Installation guide

Create a dedicated Conda environment for the package:
```bash
conda create -n rxn-standardization python=3.7
conda activate rxn-standardization
```

Install the package and its dependencies via the following command:
```bash
pip install -e .
```

# Training the transformer model for standardization

This section explains how to preprocess input data, and train, test and evaluate the translation model for standardization.

## Preprocessing

Note: More information on how to acquire the datasets used in the manuscript and extract the relevant information into CSV files can be found [here](./resources/README.md).

For simplicity, set the following environment variable:
```bash
export DATA_DIR="$(pwd)/data"
```
`DATA_DIR` can be changed to any other location containing the data to train and test on.

First, split your dataset and generate source and target files by running:
```bash
rxn-std-process-csv --input_csv <input_file_path> --save_dir $DATA_DIR
```
To see details on the required input file format and all options available when running `process-csv`, run:
```bash
rxn-std-process-csv --help
```

`DATA_DIR` will then contain the following files, with *tokenized* SMILES:
```bash
src-test.txt    src-train.txt   src-valid.txt   tgt-test.txt    tgt-train.txt   tgt-valid.txt
```

To perform multiple dataset splits for cross-validation, run:
```bash
rxn-std-split-for-cv --input_csv <input_file_path> --save_dir $DATA_DIR
```
`DATA_DIR` will then contain 5 src/tgt files with different splits, with *tokenized* SMILES. To see all options available when performing the splits (augmentation, prepending tokens, specifying test size), run:
```bash
rxn-std-split-for-cv --help
```

## Training

Convert the data to the format required by OpenNMT:
```bash
onmt_preprocess \
  -train_src $DATA_DIR/src-train.txt -train_tgt $DATA_DIR/tgt-train.txt \
  -valid_src $DATA_DIR/src-valid.txt -valid_tgt $DATA_DIR/tgt-valid.txt \
  -save_data $DATA_DIR/preprocessed -src_seq_length 500 -tgt_seq_length 500 \
  -src_vocab_size 500 -tgt_vocab_size 500 -share_vocab
```

To then train the transformer model with OpenNMT: 
```bash
onmt_train \
  -data $DATA_DIR/preprocessed  -save_model  $DATA_DIR/models/model  \
  -seed 42 -save_checkpoint_steps 10000 -keep_checkpoint -1 \
  -train_steps 200000 -param_init 0  -param_init_glorot -max_generator_batches 32 \
  -batch_size 4096 -batch_type tokens -normalization tokens -max_grad_norm 0  -accum_count 4 \
  -optim adam -adam_beta1 0.9 -adam_beta2 0.998 -decay_method noam -warmup_steps 8000  \
  -learning_rate 2 -label_smoothing 0.0 -report_every 1000  -valid_batch_size 8 \
  -layers 4 -rnn_size 256 -word_vec_size 256 -encoder_type transformer -decoder_type transformer \
  -dropout 0.1 -position_encoding -share_embeddings \
  -global_attention general -global_attention_function softmax -self_attn_type scaled-dot \
  -heads 8 -transformer_ff 2048
```
Training the model can take up to a one day in a GPU-enabled environment for ~200k data points. To enable GPU usage, run the same command with `-gpu_ranks 0`

## Finetuning

For finetuning, first preprocess the data using the `onmt_preprocess` command above. Assuming that the preprocessed data for finetuning has been saved with the prefix `$DATA_DIR/preprocessed_finetuning`, first extend the pretrained model's vocabulary using:

```bash
rxn-extend-model-with-vocab --model_path $DATA_DIR/models_finetuning/model_step_120000.pt --vocab_path $DATA_DIR/preprocessed_finetuning/preprocessed.vocab.pt --output_path $DATA_DIR/models_finetuning/model_step_120000_extended_vocab.pt
```

Then, use the same training command with slightly different parameters:
```bash
onmt_train \
  -data $DATA_DIR/preprocessed_finetuning  \
  -train_from $DATA_DIR/models_finetuning/model_step_120000_extended_vocab.pt \
  -save_model  $DATA_DIR/models/model  \
  -seed 42 -save_checkpoint_steps 5000 -keep_checkpoint -1 \
  -train_steps 50000 -param_init 0  -param_init_glorot -max_generator_batches 32 \
  -batch_size 4096 -batch_type tokens -normalization tokens -max_grad_norm 0  -accum_count 4 \
  -optim adam -adam_beta1 0.9 -adam_beta2 0.998 -decay_method noam -warmup_steps 8000  \
  -learning_rate 2 -label_smoothing 0.0 -report_every 200  -valid_batch_size 8 \
  -layers 4 -rnn_size 256 -word_vec_size 256 -encoder_type transformer -decoder_type transformer \
  -dropout 0.1 -position_encoding -share_embeddings \
  -global_attention general -global_attention_function softmax -self_attn_type scaled-dot \
  -heads 8 -transformer_ff 2048
```

## Standardize SMILES with the transformer model

SMILES strings can be translated to their standardized format with the following:
```bash
# Update the path to the OpenNMT model as required
export MODEL="$DATA_DIR/models/model_step_120000.pt"
# Specify whether to report top-1/top-2/etc. accuracy
export TOPN=1

onmt_translate -model $MODEL -src $DATA_DIR/src-test.txt -tgt $DATA_DIR/tgt-test.txt -output $DATA_DIR/pred.txt -log_probs -n_best $TOPN -beam_size 10 -max_length 300 -batch_size 10

# To detokenize the predictions file:
rxn-std-process-output --input_file $DATA_DIR/pred.txt --output_file $DATA_DIR/pred_detok.txt --canonicalize_output
```

## Evaluation

To print the metrics on the predictions, the following command can be used:
```bash
rxn-std-score-predictions --pred_file $DATA_DIR/pred.txt --tgt_file $DATA_DIR/tgt-test.txt 
```

To see all options available when obtaining the metrics, run:
```bash
rxn-std-score-predictions --help
```
