# Similar Predictive Performance and Clinical Utility of the Kidney Failure Risk Equation Using EKFC Or CKD-EPI Estimated Glomerular Filtration Rate

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Authors

This work was co-first authored by [**Antoine Créon**](https://github.com/Antoine-Creon) and **Malou Magnani**, and supervised by [**Edouard Fu**](https://edouard-fu.github.io/) (Leiden University Medical Center) in collaboration with the [Cardio-Renal Epidemiology group](https://ki.se/en/research/research-areas-centres-and-networks/research-groups/cardio-renal-epidemiology-juan-jesus-carreros-research-group) of **Juan Jesus Carrero** at Karolinska Institutet.

## Overview

This repository contains the analysis code for the study *"Similar Predictive Performance and Clinical Utility of the Kidney Failure Risk Equation Using EKFC Or CKD-EPI Estimated Glomerular Filtration Rate"*, published in *[JOURNAL NAME]* in *[PUBLICATION DATE]*.

*[LINK TO THE PUBLISHED MANUSCRIPT]*

## Repository Structure

```
├── Code/
│   ├── 1. Cohort derivation.R       # Define study cohort and inclusion criteria
│   ├── 2. Covariate derivation.R    # Derive baseline covariates
│   ├── 3. Outcome derivation.R      # Define kidney failure and death outcomes
│   ├── 4. Predictions.R             # Calculate KFRE predictions using different eGFR equations
│   ├── 5. Validation.R              # Model validation (discrimination, calibration, Brier scores)
│   ├── 6.1 Sensitivity analysis - restricted to two eGFRs.R              # Sensitivity analysis
│   ├── 6.2 Sensitivity analysis - subgroups.R              # Sensitivity analysis
│   ├── Functions eGFR equations.R   # eGFR calculation functions (CKD-EPI, EKFC)
│   ├── Functions for analyses.R     # Helper functions for statistical analyses
└── README.md
```

## Requirements

### R Version
- R 4.5.1

### Required Packages

| Package | Version | Purpose |
|---------|---------|---------|
| here | 1.0.2 | Project-relative paths |
| dplyr | 1.1.4 | Data manipulation |
| tidyr | 1.3.1 | Data reshaping |
| purrr | 1.1.0 | Functional programming |
| stringr | 1.5.2 | String manipulation |
| forcats | 1.0.1 | Factor handling |
| survival | 3.8.3 | Survival analysis |
| cmprsk | 2.2.12 | Competing risks analysis |
| riskRegression | 2025.9.17 | Risk prediction models |
| timeROC | 0.4 | Time-dependent ROC curves |
| dcurves | 0.5.1 | Decision curve analysis |
| scales | 1.4.0 | Axis scaling |


## Data Availability

This study uses data from the Stockholm CREAtinine Measurement (SCREAM) cohort. Due to Swedish and European privacy regulations, individual-level data cannot be shared publicly. The code is provided for transparency purposes. Researchers interested in accessing SCREAM data should contact the Steering Committee of the SCREAM project (juan.jesus.carrero@ki.se).


## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt this material, provided you give appropriate credit to the original authors.

## Citation

If you use this code, please cite:

> XXX (ADD CITATION ONCE PUBLISHED)
