################################################################################
# FILENAME: 6.1 Sensitivity analysis - restricted to two eGFRs.R
# PROJECT: Performance and Clinical Utility of the KFRE Using EKFC Versus CKD-EPI eGFR
# PURPOSE: Sensitivity analysis restricted to two eGFR values <60 ml/min/1.73m2
# AUTHOR: Antoine Creon
# CREATED: 2026-03-10

# R VERSION: 4.5
################################################################################

# remove history
rm(list = ls(all.names = TRUE))

# set seed for reproducibility
set.seed(27)

# load libraries
library(dplyr)
library(here)
library(purrr)
library(survival) # time-to-event analyses
library(stringr)
library(timeROC)
library(dcurves)


# load data
load(here("Data", "raw", "albuminuria_lab_test_clean.Rda"))
load(here("Data", "raw", "death.Rda"))
load(here("Data", "raw", "demographics.Rda"))
load(here("Data", "raw", "lab_values.Rda"))
load(here("Data", "raw", "rrt.Rda"))

# load functions
source(here("Code", "Functions eGFR equations.R"))
source(here("Code", "Functions for analyses.R"))

# Load study population
load(file = here("Data", "clean", "cohort_outcomes.RData"))

################################################################################
# Filter patients with a history of eGFR<60 ml/min/1.73m2 ######################
################################################################################

## collect outpatient creatinine -----------------------------------------------

prev_creatinine <- lab_values |>
  dplyr::filter(test == "crea") |>
  dplyr::mutate(datum = as.Date(datum, format = "%Y-%m-%d")) |>
  dplyr::rename(date_prev_creat = datum) |>
  dplyr::rename(prev_creat = result) |>
  dplyr::filter(ip == 0) |>
  dplyr::select(lopnr, date_prev_creat, prev_creat) |>
  dplyr::arrange(lopnr, date_prev_creat)

## Add eGFR measurements before the qualifying eGFR to the dataset -------------

# Add creatinine, age and sex, and compute previous eGFR using CKD-EPI 2009
creat_history <- cohort |>
  left_join(
    prev_creatinine,
    join_by(lopnr == lopnr, date_creat > date_prev_creat)
  ) |>
  left_join(demo) |>
  ungroup() |>
  mutate(
    dob = as.Date(dob, format = "%Y-%m-%d"),
    age_prev_creat = as.numeric(difftime(
      date_prev_creat,
      dob,
      units = "days"
    )) /
      365.25,
    prev_eGFR = ckd_epi_2009_cr(prev_creat, age_prev_creat, female)
  ) |>
  select(
    lopnr,
    date_prev_creat,
    date_creat,
    index_dt,
    prev_eGFR,
    ckd_epi_2009_cr
  )

# Identify all the patients with at least one eGFR < 60 ml/min/1.73m2
# at least 3 months before the qualifying eGFR
confirmed_CKD_ids <- creat_history |>
  filter(
    as.numeric(date_creat - date_prev_creat, units = "secs") >
      lubridate::dmonths(3)
  ) |>
  group_by(lopnr) |>
  mutate(egfr_below_60 = prev_eGFR < 60) |>
  filter(egfr_below_60 == TRUE) |>
  pull(lopnr) |>
  unique()

# Filter the cohort to include only patients with confirmed CKD
confirmed_CKD <- cohort |>
  filter(lopnr %in% confirmed_CKD_ids)

# Save
save(confirmed_CKD, file = here("Data", "clean", "confirmed_CKD.RData"))


################################################################################
# Predictions ##################################################################
################################################################################

# # List eGFRs
egfrs <- c(
  "ckd_epi_2009_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2012_cr_cys",
  "ekfc_cr",
  "ekfc_cys",
  "ekfc_cr_cys"
)

## Prognostic index according to different equations ---------------------------

confirmed_CKD <- confirmed_CKD |>
  mutate(
    male = ifelse(female == 0, 1, 0), # create male variable
    alb = ifelse(alb == 0, 1, alb)
  ) |>
  mutate(across(
    all_of(egfrs),
    ~ PI_KFRE(age, male, .x, alb),
    .names = "PI_{.col}"
  ))


## Predicted 2-year KFRE risk according to different equations -----------------

# Name columns
PIs <- glue::glue("PI_{egfrs}")

confirmed_CKD <- confirmed_CKD |>
  mutate(across(
    all_of(PIs),
    ~ KFRE_risk_2y(.x),
    .names = "risk_2y_{egfrs}"
  ))

## Predicted 5-year KFRE risk according to different equations -----------------

confirmed_CKD <- confirmed_CKD |>
  mutate(across(
    all_of(PIs),
    ~ KFRE_risk_5y(.x),
    .names = "risk_5y_{egfrs}"
  ))

# save confirmed_CKD with predicted probabilities
save(
  confirmed_CKD,
  file = here("Data", "clean", "confirmed_CKD_predictions.RData")
)

################################################################################
# Validation ###################################################################
################################################################################

## Prepapre analyses -----------------------------------------------------------

equations <- c(
  "ckd_epi_2009_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2012_cr_cys",
  "ekfc_cr",
  "ekfc_cys",
  "ekfc_cr_cys"
)

equation_names <- c(
  "CKD-EPI[cr]~2009",
  "CKD-EPI[cys]~2012",
  "CKD-EPI[cr-cys]~2012",
  "EKFC[cr]",
  "EKFC[cys]",
  "EKFC[cr-cys]"
)

horizons <- c(2, 5)

# Make a table of results
results <- tidyr::crossing(equations, horizons)


## Area under the curve --------------------------------------------------------

AUCs <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(data = confirmed_CKD, horizon = .x, equation = .y)
) |>
  bind_rows()

AUCs <- cbind(results, AUCs)


# Save the results tables (long computation time)
save(AUCs, file = here("Data", "analyses", "confirmed_CKD_AUCs.rda"))

## Calibration plots -----------------------------------------------------------

calibration_plot_tbl <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(data = confirmed_CKD, equation = .x, horizon = .y)
) |>
  list_rbind()

calibration_plot_tbl <- calibration_plot_tbl %>%
  mutate(
    color_group = paste0(model, "_", horizon),
    plot_group = case_when(
      model %in%
        c("ckd_epi_2009_cr", "ekfc_cr", "revised_lm_cr") ~ "Creatinine",
      model %in%
        c("ckd_epi_2012_cys", "ekfc_cys", "revised_lm_cys") ~ "Cystatin C",
      model %in%
        c(
          "ckd_epi_2012_cr_cys",
          "ekfc_cr_cys",
          "revised_lm_cr_cys"
        ) ~ "Creatinine and Cystatin C",
    )
  )

save(
  calibration_plot_tbl,
  file = here("Data", "analyses", "confirmed_CKD_calibration_plots_tbl.rda")
)

## Decision curve analysis -----------------------------------------------------

# Format dataset appropriately
confirmed_CKD <-
  confirmed_CKD %>%
  mutate(
    outcome_2y_dca = factor(
      outcome_2y,
      levels = 0:2,
      labels = c("censor", "KFRT", "Death")
    ),
    outcome_5y_dca = factor(
      outcome_5y,
      levels = 0:2,
      labels = c("censor", "KFRT", "Death")
    )
  )

dca_combinations <- tidyr::crossing(
  horizon = c(2, 5),
  biomarkers = c("Creatinine", "Cystatin C", "Creatinine and Cystatin C")
) |>
  mutate(
    equation1 = case_when(
      biomarkers == "Creatinine" ~ "ckd_epi_2009_cr",
      biomarkers == "Cystatin C" ~ "ckd_epi_2012_cys",
      biomarkers == "Creatinine and Cystatin C" ~ "ckd_epi_2012_cr_cys"
    ),
    equation2 = case_when(
      biomarkers == "Creatinine" ~ "ekfc_cr",
      biomarkers == "Cystatin C" ~ "ekfc_cys",
      biomarkers == "Creatinine and Cystatin C" ~ "ekfc_cr_cys"
    )
  )

# 2-y KFRE in individuals with reduced eGFR
DCAs_2y_tbl <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = confirmed_CKD |>
        filter(ckd_epi_2009_cr >= 10, ckd_epi_2009_cr < 30),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.45
    )
  ) |>
  bind_rows()

# Format results
DCAs_2y_tbl <- DCAs_2y_tbl |>
  mutate(
    biomarker = case_when(
      str_detect(equations, "cr_cys") ~ "bold(eGFR[cr~'-'~cys])",
      str_detect(equations, "cr") ~ "bold(eGFR[cr])",
      str_detect(equations, "cys") ~ "bold(eGFR[cys])"
    ),
    biomarker = factor(
      biomarker,
      levels = c("bold(eGFR[cr])", "bold(eGFR[cys])", "bold(eGFR[cr~'-'~cys])")
    ),
    updated_label = case_when(
      str_detect(label, "ekfc") ~ "Use EKFC",
      str_detect(label, "ckd_epi") ~ "Use CKD-EPI",
      str_detect(label, "Treat All") ~ "Refer All",
      str_detect(label, "Treat None") ~ "Refer None",
      TRUE ~ label,
    )
  )

# 5-y KFRE in individuals with higher eGFR
DCAs_5y_tbl <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = confirmed_CKD |>
        filter(ckd_epi_2009_cr >= 30, ckd_epi_2009_cr < 60),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.2
    )
  ) |>
  bind_rows()

# Format results
DCAs_5y_tbl <- DCAs_5y_tbl |>
  mutate(
    biomarker = case_when(
      str_detect(equations, "cr_cys") ~ "bold(eGFR[cr~'-'~cys])",
      str_detect(equations, "cr") ~ "bold(eGFR[cr])",
      str_detect(equations, "cys") ~ "bold(eGFR[cys])"
    ),
    biomarker = factor(
      biomarker,
      levels = c("bold(eGFR[cr])", "bold(eGFR[cys])", "bold(eGFR[cr~'-'~cys])")
    ),
    updated_label = case_when(
      str_detect(label, "ekfc") ~ "Use EKFC",
      str_detect(label, "ckd_epi") ~ "Use CKD-EPI",
      str_detect(label, "Treat All") ~ "Refer All",
      str_detect(label, "Treat None") ~ "Refer None",
      TRUE ~ label,
    )
  )

save(
  DCAs_2y_tbl,
  DCAs_5y_tbl,
  file = here("Data", "analyses", "confirmed_CKD_DCA.rda")
)
