# TAS-seq: TadA-8e Deaminase-Assisted Subcellular Structural Adenosine Sequencing

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.crmeth.2025.101087-blue)](https://doi.org/10.1016/j.crmeth.2025.101087)

Repository for the TAS-seq computational pipeline and analysis scripts from the publication:

> **TadA-8e deaminase-assisted subcellular structural adenosine sequencing**  
> *DOI: 10.1016/j.crmeth.2025.101087*

---

## Repository Structure
├── Perl_scripts/ # Perl analysis scripts
├── Prediction_model/ # Machine learning prediction models
│ └── README.md # Model documentation
├── R_scripts/ # R analysis and visualization scripts
├── run.sh # Main pipeline execution script
└── README.md # This overview document


## Requirements
- Unix/Linux environment
- Perl 5.10+
- R 4.0+
- Python 3.10+ (for prediction models)

## Quick Start
Execute the main pipeline:
```bash
chmod +x run.sh  # Set execute permission
./run.sh
