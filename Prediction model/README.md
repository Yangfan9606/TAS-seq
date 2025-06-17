# TAS-seq Prediction Model

## Overview
This repository contains tools and a pre-trained model for sequence classification in TAS-seq analysis. The core model was trained on randomly sampled biological sequences.

## Key Files
- **`run.sh`**  
  Records the execution pipeline and commands used for running the training/prediction workflow.
  
- **`ma_classifier_model.pth`**  
  Pre-trained PyTorch model for sequence classification. Trained on **200,000 randomly selected sequences** from HEK293T NES data.

## Additional Scripts
| File | Status | Description |
|------|--------|-------------|
| `GetSeqFromFa.pl` | _Ready to use_ | FASTA sequence extractor |
| `Negative_A_sites_from_HEK293T_NES` | _Ready to use_ | Negative modificated sites |
| `Positive_A_sites_from_HEK293T_NES` | _Ready to use_ | Positive modificated sites |
| `Sequence_classifier.py` | _Ready to use_ | Model training script |
| `Sequence_predictor.py` | _Ready to use_ | Prediction script |
| `fna.NCBI.Chr_Rename.pl` | _Ready to use_ | Chromosome name normalizer |
| `yangfan.pm` | _Ready to use_ | Perl utility modules |

---
*Developed by Yangfan • Model trained on HEK293T NES dataset*
