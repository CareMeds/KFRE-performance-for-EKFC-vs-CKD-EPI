################################################################################
# FILENAME: 5. Validation.R
# PROJECT: Performance and Clinical Utility of the KFRE Using EKFC Versus CKD-EPI eGFR
# PURPOSE: Validate KF predictions
# AUTHOR: Malou Magnani
# CREATED: 2025-02
# UPDATED:
# - 2025-10-22 (Antoine Creon)

# R VERSION: 4.5
################################################################################
# remove history
rm(list = ls(all.names = TRUE))

# load libraries
library(survival) # time-to-event analyses
library(patchwork) # combine plots
library(here) # combine plots
library(purrr) # combine plots
library(dplyr) # combine plots
library(geepack)
library(boot)
library(parallel)
library(stringr)
library(timeROC)
library(dcurves)

# set seed for reproducibility
B <- 500
seed <- 27
set.seed(seed)

# load data sets
load(here("Data", "clean", "cohort_predictions.RData"))

# load functions
source(here("Code", "Functions for analyses.R"))

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

################################################################################
# For each equation, calculate performance measures ############################
################################################################################

## Area under the curve --------------------------------------------------------

AUCs <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(data = cohort, horizon = .x, equation = .y)
) |>
  bind_rows()

# Add them to the results table
results <- cbind(results, AUCs)

# Save the results tables (long computation time)
save(results, file = here("Data", "analyses", "res_table.rda"))

## Calibration intercept and slope ---------------------------------------------
# load(file = here("Data", "analyses", "res_table.rda"))

# Compute calibration intercept, slope and their CI for each combination of equation and horizon
calibration <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ cal_int_slope(data = cohort, equation = .x, horizon = .y)
) |>
  list_rbind()

# Add calibration to the results table
results <- cbind(results, calibration)

## O/E ratio -------------------------------------------------------------------
# load(file = here("Data", "analyses", "res_table.rda"))

# Put CIF in function, change output to mimick the two above

OE <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ oe_ratio(equation = .x, horizon = .y, bootstrap = TRUE, B = B, seed = seed)
) |>
  list_rbind()

results <- cbind(results, OE)

## Brier scores ----------------------------------------------------------------

Brier_scores_tbl <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ brier_scores(data = cohort, horizon = .x, equation = .y)
) |>
  bind_rows()

results <- cbind(results, Brier_scores_tbl)

# Save the results tables
save(results, file = here("Data", "analyses", "res_table.rda"))

################################################################################
# Compute comparison metrics between equations #################################
################################################################################

delta_brier_tbl <- tidyr::crossing(
  biomarker = c("eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]"),
  horizon = c(2, 5)
) |>
  mutate(
    biomarker = factor(
      biomarker,
      levels = c("eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]")
    )
  ) |>
  arrange(biomarker, horizon) |>
  mutate(
    equation_ref = case_when(
      biomarker == "eGFR[cr]" ~ "ckd_epi_2009_cr",
      biomarker == "eGFR[cys]" ~ "ckd_epi_2012_cys",
      biomarker == "eGFR[cr-cys]" ~ "ckd_epi_2012_cr_cys"
    ),
    equation_test = case_when(
      biomarker == "eGFR[cr]" ~ "ekfc_cr",
      biomarker == "eGFR[cys]" ~ "ekfc_cys",
      biomarker == "eGFR[cr-cys]" ~ "ekfc_cr_cys"
    )
  )

delta_scaled_brier_res <- delta_brier_tbl %>%
  pmap(
    ~ delta_scaled_brier_boot(
      data = cohort,
      equation_ref = ..3,
      equation_test = ..4,
      horizon = ..2,
      B = 500,
      seed = 19920903
    )
  ) |>
  bind_rows()

delta_scaled_brier_res <- cbind(delta_brier_tbl, delta_scaled_brier_res)

save(
  delta_scaled_brier_res,
  file = here("Data", "analyses", "delta_scaled_brier_res.rda")
)


################################################################################
# Combined calibration plot for all models for both horizons ###################
################################################################################

calibration_plot_tbl <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(data = cohort, equation = .x, horizon = .y)
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
  file = here("Data", "analyses", "calibration_plots_tbl.rda")
)
# load(file = here("Data", "analyses", "calibration_plots_tbl.rda"))

################################################################################
# Decision curve analysis ######################################################
################################################################################

# Format dataset appropriately
cohort <-
  cohort %>%
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
      data = cohort |> filter(ckd_epi_2009_cr >= 10, ckd_epi_2009_cr < 30),
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
      data = cohort |> filter(ckd_epi_2009_cr >= 30, ckd_epi_2009_cr < 60),
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
  file = here("Data", "analyses", "decision_curve_analysis_main.rda")
)

# load(here::here("Data", "analyses", "decision_curve_analysis_main.rda"))

################################################################################
# Risk distribution plots ######################################################
################################################################################

equation_order_2y <- paste0("risk_2y_", equations)

## 2 years ---------------------------------------------------------------------

# long format
cohort_risk_2y_long <- cohort |>
  mutate(
    KFRT = if_else(
      outcome_2y == 1,
      "KFRT",
      "No KFRT"
    ),
    KFRT = factor(
      KFRT,
      levels = c("No KFRT", "KFRT")
    )
  ) |>
  dplyr::select(
    lopnr,
    KFRT,
    all_of(equation_order_2y)
  ) |>
  tidyr::pivot_longer(
    cols = starts_with("risk_"),
    names_to = "Equation",
    values_to = "Risk"
  ) |>
  dplyr::mutate(
    Equation = factor(
      Equation,
      levels = equation_order_2y,
      labels = equation_names
    )
  )

save(
  cohort_risk_2y_long,
  equation_order_2y,
  file = here("Data", "analyses", "risk_2y_long.rda")
)

## 5 years ---------------------------------------------------------------------

equation_order_5y <- paste0("risk_5y_", equations)

# long format
cohort_risk_5y_long <- cohort |>
  mutate(
    KFRT = if_else(
      outcome_5y == 1,
      "KFRT",
      "No KFRT"
    ),
    KFRT = factor(
      KFRT,
      levels = c("No KFRT", "KFRT")
    )
  ) |>
  dplyr::select(
    lopnr,
    KFRT,
    all_of(equation_order_5y)
  ) |>
  tidyr::pivot_longer(
    cols = starts_with("risk_"),
    names_to = "Equation",
    values_to = "Risk"
  ) |>
  dplyr::mutate(
    Equation = factor(
      Equation,
      levels = equation_order_5y,
      labels = equation_names
    )
  )

save(
  cohort_risk_5y_long,
  equation_order_5y,
  file = here("Data", "analyses", "risk_5y_long.rda")
)


################################################################################
# Risk differences #############################################################
################################################################################

## 2 years ---------------------------------------------------------------------

cohort_risk_diff_2y_long <- cohort_risk_2y_long |>
  filter(!str_detect(Equation, "R-LM")) %>%
  mutate(
    method = case_when(
      str_detect(Equation, "CKD-EPI") ~ "CKD-EPI",
      str_detect(Equation, "EKFC") ~ "EKFC",
      TRUE ~ NA_character_
    ),
    biomarker = case_when(
      str_detect(Equation, "cr-cys\\]") ~ "eGFR[cr-cys]",
      str_detect(Equation, "cr\\]") ~ "eGFR[cr]",
      str_detect(Equation, "cys\\]") ~ "eGFR[cys]",
      TRUE ~ NA_character_
    ),
    biomarker = factor(
      biomarker,
      levels = c("eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]")
    )
  ) %>%
  select(lopnr, method, biomarker, Risk) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = Risk
  ) %>%
  mutate(
    risk_diff = EKFC - `CKD-EPI`
  )

save(
  cohort_risk_diff_2y_long,
  file = here("Data", "analyses", "riskdiff_2y.rda")
)

## 5 years ---------------------------------------------------------------------

cohort_risk_diff_5y_long <- cohort_risk_5y_long |>
  filter(!str_detect(Equation, "R-LM")) %>%
  mutate(
    method = case_when(
      str_detect(Equation, "CKD-EPI") ~ "CKD-EPI",
      str_detect(Equation, "EKFC") ~ "EKFC",
      TRUE ~ NA_character_
    ),
    biomarker = case_when(
      str_detect(Equation, "cr-cys\\]") ~ "eGFR[cr-cys]",
      str_detect(Equation, "cr\\]") ~ "eGFR[cr]",
      str_detect(Equation, "cys\\]") ~ "eGFR[cys]",
      TRUE ~ NA_character_
    ),
    biomarker = factor(
      biomarker,
      levels = c("eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]")
    )
  ) %>%
  select(lopnr, method, biomarker, Risk) %>%
  tidyr::pivot_wider(
    names_from = method,
    values_from = Risk
  ) %>%
  mutate(
    risk_diff = EKFC - `CKD-EPI`
  )

save(
  cohort_risk_diff_5y_long,
  file = here("Data", "analyses", "riskdiff_5y.rda")
)

################################################################################
# eGFR distributions ###########################################################
################################################################################

cohort_egfr_long <- cohort |>
  select(lopnr, all_of(equations[!str_detect(equations, "revised")])) |>
  tidyr::pivot_longer(
    cols = all_of(equations[!str_detect(equations, "revised")]),
    names_to = "model",
    values_to = "eGFR"
  ) |>
  dplyr::mutate(
    method = dplyr::case_when(
      str_detect(model, "^ckd_epi") ~ "CKD-EPI",
      str_detect(model, "^ekfc") ~ "EKFC",
      TRUE ~ NA_character_
    ),
    biomarker = dplyr::case_when(
      str_detect(model, "cr_cys") ~ "eGFR[cr-cys]",
      str_detect(model, "_cr") & !str_detect(model, "cr_cys") ~ "eGFR[cr]",
      str_detect(model, "_cys") & !str_detect(model, "cr_cys") ~ "eGFR[cys]",
      TRUE ~ NA_character_
    ),
    biomarker = factor(
      biomarker,
      levels = c("eGFR[cr]", "eGFR[cys]", "eGFR[cr-cys]")
    )
  ) |>
  dplyr::select(lopnr, method, biomarker, eGFR)

save(
  cohort_egfr_long,
  file = here("Data", "analyses", "cohort_egfr_long.rda")
)
