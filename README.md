# Similar Predictive Performance and Clinical Utility of the Kidney Failure Risk Equation Using EKFC Or CKD-EPI Estimated Glomerular Filtration Rate

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

This repository contains the analysis code for the study: [Créon A*, Magnani M*, Maas CHM, Van Diepen M, Caldinelli A, Russel WA, Dekker FW, Carrero JJ, Fu EL. Similar predictive performance and clinical utility of the kidney failure risk equation using EKFC or CKD-EPI estimated glomerular filtration rate. _Clinical Kidney Journal_. Published online June 8, 2026:sfag187. doi:10.1093/ckj/sfag187](https://doi.org/10.1093/ckj/sfag187)

The work was also presented as an oral communication at the [**European Renal Association 2026 Congress**](https://www.era-online.org/events/glasgow-2026/scientific-programme/), where it received the *Best Abstract Presented by Young Authors* award in the Chronic Kidney Disease category.

## Authors

This study was co-first authored by [**Antoine Créon**](https://github.com/Antoine-Creon) and **Malou Magnani**, and supervised by [**Edouard Fu**](https://edouard-fu.github.io/) (Leiden University Medical Center) in collaboration with the [Cardio-Renal Epidemiology group](https://ki.se/en/research/research-areas-centres-and-networks/research-groups/cardio-renal-epidemiology-juan-jesus-carreros-research-group) of **Juan Jesus Carrero** at Karolinska Institutet.

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

The forest plots were generated using an adapted version of the [forester](https://github.com/rdboyes/forester) package.

## Data Availability

This study uses data from the Stockholm CREAtinine Measurement (SCREAM) cohort. Due to Swedish and European privacy regulations, individual-level data cannot be shared publicly. The code is provided for transparency purposes. Researchers interested in accessing SCREAM data should contact the Steering Committee of the SCREAM project (juan.jesus.carrero@ki.se).


## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt this material, provided you give appropriate credit to the original authors.