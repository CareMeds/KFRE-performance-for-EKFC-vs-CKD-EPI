# KFRE Performance with EKFC vs CKD-EPI eGFR Equations

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Authors

This work was co-first authored by ([**Antoine Créon**](https://github.com/Antoine-Creon)) and **Malou Magnani**, supervised by [**Edouard Fu**](https://edouard-fu.github.io/) (Leiden University Medical Center) in collaboration with the [Cardio-Renal Epidemiology group](https://ki.se/en/research/research-areas-centres-and-networks/research-groups/cardio-renal-epidemiology-juan-jesus-carreros-research-group) of **Juan Jesus Carrero** at Karolinska Institutet.

## Overview

This repository contains the analysis code for a study evaluating the performance and clinical utility of the Kidney Failure Risk Equation (KFRE) when using European Kidney Function Consortium (EKFC) equations compared to the CKD-EPI 2009-2012 equations for estimating glomerular filtration rate (eGFR).

## Repository Structure

```
├── Code/
│   ├── 1. Cohort derivation.R       # Define study cohort and inclusion criteria
│   ├── 2. Covariate derivation.R    # Derive baseline covariates
│   ├── 3. Outcome derivation.R      # Define kidney failure and death outcomes
│   ├── 4. Predictions.R             # Calculate KFRE predictions using different eGFR equations
│   ├── 5. Validation.R              # Model validation (discrimination, calibration, Brier scores)
│   ├── 6. Figures and Tables.R      # Generate manuscript figures and tables
│   ├── Functions eGFR equations.R   # eGFR calculation functions (CKD-EPI, EKFC)
│   ├── Functions for analyses.R     # Helper functions for statistical analyses
└── README.md
```

## Requirements

### R Version
- R ≥ 4.0

### Required Packages

```r
# Data manipulation
install.packages(c("dplyr", "tidyr", "purrr", "stringr", "forcats"))

# Survival analysis
install.packages(c("survival", "riskRegression", "timeROC", "cmprsk"))

# Utilities
install.packages(c("here", "readxl", "scales"))
```

## Usage

Scripts are numbered to indicate execution order:

1. Run scripts `0.` through `3.` to prepare the analysis dataset
2. Run script `4.` to generate KFRE predictions
3. Run script `5.` to perform model validation

**Note:** Raw data from the SCREAM cohort are not publicly available due to privacy regulations. The code is provided for transparency purposes.

## eGFR Equations

This study compares the following eGFR equations:

| Equation | Biomarker(s) |
|----------|--------------|
| CKD-EPI 2009 | Creatinine |
| CKD-EPI 2012 | Cystatin C |
| CKD-EPI 2012 | Creatinine + Cystatin C |
| EKFC | Creatinine |
| EKFC | Cystatin C |
| EKFC | Creatinine + Cystatin C |

## Data Availability

This study uses data from the Stockholm CREAtinine Measurement (SCREAM) cohort. Due to Swedish and European privacy regulations, individual-level data cannot be shared publicly. Researchers interested in accessing SCREAM data should contact the study principal investigators.

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt this material, provided you give appropriate credit to the original authors.

## Citation

If you use this code, please cite:

> XXX (ADD CITATION ONCE PUBLISHED)

## Contact

For questions regarding the code or analysis, please open an issue in this repository.

## Acknowledgments

This study was conducted using data from the Stockholm CREAtinine Measurement (SCREAM) project.
