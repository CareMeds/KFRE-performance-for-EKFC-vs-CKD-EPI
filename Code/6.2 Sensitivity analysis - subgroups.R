################################################################################
# FILENAME: 6.2 Sensitivity analysis - subgroups.R
# PROJECT: Performance and Clinical Utility of the KFRE Using EKFC Versus CKD-EPI eGFR
# PURPOSE: Validate KF predictions
# AUTHOR: Antoine Creon
# CREATED: 2026-03-12

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
library(riskRegression)
library(timeROC)
library(dcurves)

# load functions
source(here("Code", "Functions eGFR equations.R"))
source(here("Code", "Functions for analyses.R"))

# Load study population
load(here("Data", "clean", "cohort_predictions.RData"))

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
# Create subgroups #############################################################
################################################################################

cohort <- cohort |>
  mutate(
    age_group = cut(
      age,
      breaks = c(-Inf, 75, Inf),
      labels = c("<75", ">=75")
    ),
    ckd_group = cut(
      ckd_epi_2009_cr,
      breaks = c(-Inf, 30, Inf),
      labels = c("G4-5", "G3a-3b")
    )
  ) |>
  ungroup()


################################################################################
# Area under the curve #########################################################
################################################################################

## Sex -------------------------------------------------------------------------

AUCs_male <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(
    data = cohort |> filter(female == 0),
    horizon = .x,
    equation = .y
  )
) |>
  bind_rows() |>
  mutate(subgroup = "Male") %>%
  cbind(results, .)

AUCs_female <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(
    data = cohort |> filter(female == 1),
    horizon = .x,
    equation = .y
  )
) |>
  bind_rows() |>
  mutate(subgroup = "Female") |>
  bind_rows() %>%
  cbind(results, .)

## Age -------------------------------------------------------------------------

AUCs_younger <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(
    data = cohort |> filter(age_group == "<75"),
    horizon = .x,
    equation = .y
  )
) |>
  bind_rows() |>
  mutate(subgroup = "<75 years") |>
  bind_rows() %>%
  cbind(results, .)

AUCs_older <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(
    data = cohort |> filter(age_group == ">=75"),
    horizon = .x,
    equation = .y
  )
) |>
  bind_rows() |>
  mutate(subgroup = ">=75 years") |>
  bind_rows() %>%
  cbind(results, .)

## CKD stage -------------------------------------------------------------------

AUCs_g4_5 <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(
    data = cohort |> filter(ckd_group == "G4-5"),
    horizon = .x,
    equation = .y
  )
) |>
  bind_rows() |>
  mutate(subgroup = "G4-5") |>
  bind_rows() %>%
  cbind(results, .)

AUCs_g3a_3b <- map2(
  .x = results$horizons,
  .y = results$equations,
  ~ compute_AUC(
    data = cohort |> filter(ckd_group == "G3a-3b"),
    horizon = .x,
    equation = .y
  )
) |>
  bind_rows() |>
  mutate(subgroup = "G3a-3b") |>
  bind_rows() %>%
  cbind(results, .)

AUCs <- bind_rows(
  AUCs_male,
  AUCs_female,
  AUCs_younger,
  AUCs_older,
  AUCs_g4_5,
  AUCs_g3a_3b
)

# Save the results tables (long computation time)
save(AUCs, file = here("Data", "analyses", "subgroups_AUCs.rda"))

################################################################################
# Calibration plots ############################################################
################################################################################

## Sex -------------------------------------------------------------------------

calibration_plot_tbl_male <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(
    data = cohort |> filter(female == 0),
    equation = .x,
    horizon = .y
  )
) |>
  list_rbind() |>
  mutate(subgroup = "Male")

calibration_plot_tbl_female <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(
    data = cohort |> filter(female == 1),
    equation = .x,
    horizon = .y
  )
) |>
  list_rbind() |>
  mutate(subgroup = "Female")

## Age -------------------------------------------------------------------------

calibration_plot_tbl_younger <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(
    data = cohort |> filter(age_group == "<75"),
    equation = .x,
    horizon = .y
  )
) |>
  list_rbind() |>
  mutate(subgroup = "<75 years")

calibration_plot_tbl_older <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(
    data = cohort |> filter(age_group == ">=75"),
    equation = .x,
    horizon = .y
  )
) |>
  list_rbind() |>
  mutate(subgroup = ">=75 years")

## CKD stage -------------------------------------------------------------------

calibration_plot_tbl_g4_5 <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(
    data = cohort |> filter(ckd_group == "G4-5"),
    equation = .x,
    horizon = .y
  )
) |>
  list_rbind() |>
  mutate(subgroup = "G4-5")

calibration_plot_tbl_g3a_3b <- map2(
  .x = results$equations,
  .y = results$horizons,
  ~ calibration_plot(
    data = cohort |> filter(ckd_group == "G3a-3b"),
    equation = .x,
    horizon = .y
  )
) |>
  list_rbind() |>
  mutate(subgroup = "G3a-3b")

calibration_plot_tbl <- bind_rows(
  calibration_plot_tbl_male,
  calibration_plot_tbl_female,
  calibration_plot_tbl_younger,
  calibration_plot_tbl_older,
  calibration_plot_tbl_g4_5,
  calibration_plot_tbl_g3a_3b
)

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
  file = here("Data", "analyses", "subgroups_calibration_plots_tbl.rda")
)

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

# Set combinations
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

## 2 years =====================================================================

### Sex -------------------------------------------------------------------------

DCAs_2y_tbl_male <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 10, ckd_epi_2009_cr < 30) |>
        filter(female == 0),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.45
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = "Male")

DCAs_2y_tbl_female <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 10, ckd_epi_2009_cr < 30) |>
        filter(female == 1),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.45
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = "Female")


### Age -------------------------------------------------------------------------

DCAs_2y_tbl_younger <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 10, ckd_epi_2009_cr < 30) |>
        filter(age_group == "<75"),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.45
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = "<75")

DCAs_2y_tbl_older <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 10, ckd_epi_2009_cr < 30) |>
        filter(age_group == ">=75"),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.45
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = ">=75")

### CKD stage -------------------------------------------------------------------

# The analysis is already stratified by CKD stage

### Format results -------------------------------------------------------------

DCAs_2y_tbl <- bind_rows(
  DCAs_2y_tbl_male,
  DCAs_2y_tbl_female,
  DCAs_2y_tbl_younger,
  DCAs_2y_tbl_older
)

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

## 5 years =====================================================================

### Sex -------------------------------------------------------------------------

DCAs_5y_tbl_male <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 30, ckd_epi_2009_cr < 60) |>
        filter(female == 0),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.2
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = "Male")

DCAs_5y_tbl_female <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 30, ckd_epi_2009_cr < 60) |>
        filter(female == 1),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.2
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = "Female")

### Age -------------------------------------------------------------------------

DCAs_5y_tbl_younger <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 30, ckd_epi_2009_cr < 60) |>
        filter(age_group == "<75"),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.2
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = "<75")

DCAs_5y_tbl_older <- dca_combinations |>
  pmap(
    ~ compute_dca(
      data = cohort |>
        filter(ckd_epi_2009_cr >= 30, ckd_epi_2009_cr < 60) |>
        filter(age_group == ">=75"),
      equation1 = ..3,
      equation2 = ..4,
      horizon = ..1,
      max_threshold = 0.2
    )
  ) |>
  bind_rows() |>
  mutate(subgroup = ">=75")

### CKD stage ------------------------------------------------------------------

# The analysis is already stratified by CKD stage

### Format results -------------------------------------------------------------

DCAs_5y_tbl <- bind_rows(
  DCAs_5y_tbl_male,
  DCAs_5y_tbl_female,
  DCAs_5y_tbl_younger,
  DCAs_5y_tbl_older
)

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
  file = here("Data", "analyses", "subgroups_DCA.rda")
)
