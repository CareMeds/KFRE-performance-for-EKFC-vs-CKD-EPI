################################################################################
# FILENAME: 4. Predictions.R
# PROJECT: Performance and Clinical Utility of the KFRE Using EKFC Versus CKD-EPI eGFR
# PURPOSE: Predict KF probabilities for each individual
# AUTHOR: Malou Magnani
# CREATED: 2025-02
# UPDATED:
# - 2025-10-22 (Antoine Creon)
# - 2025-11-05 (Antoine Creon)

# R VERSION: 4.5
################################################################################

# remove history
rm(list = ls(all.names = TRUE))

# set seed for reproducibility
set.seed(27)

# load libraries
library(cmprsk) # competing risk
library(dplyr) # data manipulation
library(here)

# load functions
source(here("Code", "Functions for analyses.R"))

# load data sets
load(here("Data", "clean", "cohort_outcomes.RData"))

# # List eGFRs
egfrs <- c(
  "ckd_epi_2009_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2012_cr_cys",
  "ekfc_cr",
  "ekfc_cys",
  "ekfc_cr_cys"
)

################################################################################
# Cumulative incidence rates ###################################################
################################################################################
cin_2y <- cmprsk::cuminc(
  ftime = cohort$time_to_event_2y,
  fstatus = cohort$outcome_2y,
  rho = 0,
  cencode = 0
)
cin_5y <- cmprsk::cuminc(
  ftime = cohort$time_to_event_5y,
  fstatus = cohort$outcome_5y,
  rho = 0,
  cencode = 0
)


################################################################################
# Prognostic index according to different equations ############################
################################################################################

cohort <- cohort |>
  mutate(
    male = ifelse(female == 0, 1, 0), # create male variable
    alb = ifelse(alb == 0, 1, alb)
  ) |>
  mutate(across(
    all_of(egfrs),
    ~ PI_KFRE(age, male, .x, alb),
    .names = "PI_{.col}"
  ))


################################################################################
# Predicted 2-year KFRE risk according to different equations ##################
################################################################################
# Name columns
PIs <- glue::glue("PI_{egfrs}")

cohort <- cohort |>
  mutate(across(
    all_of(PIs),
    ~ KFRE_risk_2y(.x),
    .names = "risk_2y_{egfrs}"
  ))


################################################################################
# Predicted 5-year KFRE risk according to different equations ##################
################################################################################

cohort <- cohort |>
  mutate(across(
    all_of(PIs),
    ~ KFRE_risk_5y(.x),
    .names = "risk_5y_{egfrs}"
  ))

# save cohort with predicted probabilities
save(
  cohort,
  cin_2y,
  cin_5y,
  file = here("Data", "clean", "cohort_predictions.RData")
)
