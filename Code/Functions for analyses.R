################################################################################
# FILENAME: Functions for analyses.R
# PROJECT: KFRE Performance with non-CKD-EPI equations
# PURPOSE: Compute estimated glomerular filtration rate
# AUTHOR: Malou Magnani
# CREATED: 2025-02
# UPDATED:
# - 2025-10-22 (Antoine Creon)

# R VERSION: 4.5
################################################################################

################################################################################
# KFRE-predicted risk      #####################################################
################################################################################

## Compute the linear predictor (PI) -------------------------------------------
# It is the same for all time horizons
PI_KFRE <- function(age, male, egfr, alb) {
  PI <- -0.2201 *
    (age / 10 - 7.036) +
    0.2467 * (male - 0.5642) -
    0.5567 * (egfr / 5 - 7.222) +
    0.4510 * (log(alb) - 5.137)

  return(PI)
}

## Compute 2y KFRE risk (P) ----------------------------------------------------
KFRE_risk_2y <- function(PI) {
  risk_2y <- 1 - (0.9832**exp(PI)) # re-calibrated for Non-North America

  return(risk_2y)
}

## Compute 5y KFRE risk (P) ----------------------------------------------------
KFRE_risk_5y <- function(PI) {
  risk_5y <- 1 - (0.9365**exp(PI)) # re-calibrated for non-north america

  return(risk_5y)
}

################################################################################
# Area under the curve     #####################################################
################################################################################

compute_AUC <- function(data, horizon, equation) {
  equation_risk <- paste0("risk_", horizon, "y_", equation)

  roc_object <- timeROC(
    T = data[["time_to_event_full"]],
    delta = data[["outcome_full"]],
    marker = data[[equation_risk]], # predicted risk at 5 years
    cause = 1,
    weighting = "marginal",
    times = 365.25 * horizon,
    iid = TRUE
  )

  .auc <- roc_object$AUC_1[[2]]
  .auc_lower <- roc_object$AUC_1[[2]] -
    1.96 * roc_object$inference$vect_sd_1[[2]]
  .auc_upper <- roc_object$AUC_1[[2]] +
    1.96 * roc_object$inference$vect_sd_1[[2]]

  return(c("AUC" = .auc, AUC_upper = .auc_lower, AUC_lower = .auc_upper))
}


################################################################################
# Calibration intercept and slope ##############################################
################################################################################

## Compute calibration intercept & slope with 95% CIs --------------------------

cal_int_slope <- function(data, equation, horizon) {
  # extract outcome variables according to horizon
  pred_risks <- eval(parse(
    text = paste0("cohort$risk_", horizon, "y_", equation)
  ))
  data$time_to_event <- eval(parse(
    text = paste0("data$time_to_event_", horizon, "y")
  ))
  data$outcome <- eval(parse(text = paste0("data$outcome_", horizon, "y")))

  Score <- riskRegression::Score(
    # Give score list of predictions for all individuals
    # (This will be used for the predicted part of validation)
    list(pred_risks),
    # Define the model as subdistributional hazards model
    formula = Hist(time_to_event, outcome) ~ 1,
    # specify pseudovalues to be used
    cens.method = "pseudo",
    # define validation dataset with values for covariates,
    # this is to be used for observed probabilities
    data = data,
    # define prediction horizon
    times = horizon * 365.25,
    # define event of interest
    outcome = 1,
    conf.int = TRUE,
    # define validation methods
    plots = "calibration",
    summary = c("ipa"),
    metrics = "brier"
  )

  # Extract pseudo-values
  pseudos <- data.frame(Score$Calibration$plotframe)

  # add the cloglog risk estimates
  pseudos <- pseudos |>
    mutate(cll_pred = log(-log(1 - risk)))

  # Fit model for calibration intercept
  fit_cal_int <- geese(
    pseudovalue ~ offset(cll_pred),
    data = pseudos,
    id = riskRegression_ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence",
    jack = TRUE
  )

  # Fit model for calibration slope
  fit_cal_slope <- geepack::geese(
    pseudovalue ~ offset(cll_pred) + cll_pred,
    data = pseudos,
    id = riskRegression_ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence",
    jack = TRUE
  )

  # Extract results
  res_cal_int <- summary(fit_cal_int)$mean
  res_cal_slope <- summary(fit_cal_slope)$mean["cll_pred", ]

  # Return estimates with CIs
  return(tibble(
    cal_int = res_cal_int$estimate,
    cal_int_lower = res_cal_int$estimate - qnorm(0.975) * res_cal_int$san.se,
    cal_int_upper = res_cal_int$estimate + qnorm(0.975) * res_cal_int$san.se,
    cal_slope = 1 + res_cal_slope$estimate,
    cal_slope_lower = 1 +
      (res_cal_slope$estimate - qnorm(0.975) * res_cal_slope$san.se),
    cal_slope_upper = 1 +
      (res_cal_slope$estimate + qnorm(0.975) * res_cal_slope$san.se),
  ))
}

## Calibration plots -----------------------------------------------------------

calibration_plot <- function(data, equation, horizon) {
  .risk <- paste0("risk_", horizon, "y_", equation)

  .tte <- paste0("time_to_event_", horizon, "y")
  .outcome <- paste0("outcome_", horizon, "y")

  .cal_formula <- as.formula(paste0("Hist(", .tte, ",", .outcome, ") ~ 1"))

  Score <- riskRegression::Score(
    object = list(equation = data[[.risk]]),
    formula = .cal_formula,
    cens.method = "pseudo",
    data = data,
    times = horizon * 365.25,
    outcome = 1,
    conf.int = TRUE,
    plots = "calibration"
  )

  # Extract and smooth pseudo-values
  pseudos <- data.frame(Score$Calibration$plotframe) |>
    dplyr::arrange(risk)

  # LOESS model
  loess_model <- stats::loess(
    pseudovalue ~ risk,
    data = pseudos,
    degree = 1,
    span = 0.33
  )

  # add CI
  smooth_pseudos <- predict(loess_model, se = TRUE)

  # Store in a data frame
  temp_df <- data.frame(
    risk = pseudos$risk,
    observed = smooth_pseudos$fit,
    observed_upper = smooth_pseudos$fit + 1.96 * smooth_pseudos$se.fit,
    observed_lower = smooth_pseudos$fit - 1.96 * smooth_pseudos$se.fit,
    model = equation
  ) %>%
    mutate(horizon = horizon)

  temp_df
}


################################################################################
# Observed/estimated ratio #####################################################
################################################################################

oe_ratio <- function(
  equation,
  horizon,
  bootstrap = FALSE,
  B = 500,
  seed = 123
) {
  # Define variables
  pred_risks <- eval(parse(
    text = paste0("cohort$risk_", horizon, "y_", equation)
  ))

  CIF <- eval(parse(text = paste0("cin_", horizon, "y")))

  # Calculate original O/E ratio
  expected <- mean(pred_risks)
  observed <- cmprsk::timepoints(CIF, horizon * 365.25)
  oe <- observed$est[1] / expected

  # Initialize results with just the O/E ratio
  results <- tibble(
    O_E_Ratio = oe
  )

  # Add bootstrapping if requested
  if (bootstrap) {
    # Define bootstrap function
    boot_func <- function(data, indices, pred_risks, cif_obj, time_horizon) {
      # Resample risk predictions
      resampled_risks <- pred_risks[indices]

      # Calculate expected risk in bootstrap sample
      boot_expected <- mean(resampled_risks)

      # For observed, we need to recalculate CIF with resampled data
      # This depends on how your CIF object is created
      # This is a placeholder - you'll need to adjust based on your data structure
      boot_observed <- cmprsk::timepoints(cif_obj, time_horizon)

      # Calculate O/E ratio
      boot_oe <- boot_observed$est[1] / boot_expected

      return(as.numeric(boot_oe))
    }

    # Set seed for reproducibility
    set.seed(seed)

    # Create a dummy dataset for bootstrapping
    # We need this because boot() requires a data argument
    dummy_data <- data.frame(id = 1:length(pred_risks))

    # Run bootstrap
    boot_results <- boot::boot(
      data = dummy_data,
      statistic = boot_func,
      R = B,
      pred_risks = pred_risks,
      cif_obj = CIF,
      time_horizon = horizon * 365.25
    )

    # Calculate bootstrap confidence intervals
    oe_ci <- boot::boot.ci(boot_results, type = "perc")

    # Add only the CI to the results
    results$O_E_lower <- if (!is.null(oe_ci)) {
      oe_ci$percent[4]
    } else {
      NULL
    }

    results$O_E_upper <- if (!is.null(oe_ci)) {
      oe_ci$percent[5]
    } else {
      NULL
    }
  }

  return(results)
}

################################################################################
# Scaled Brier & Brier score, with CIs #########################################
################################################################################

brier_scores <- function(
  data,
  equation,
  horizon
) {
  # extract outcome variables according to horizon
  .pred_risks <- paste0("risk_", horizon, "y_", equation)
  .tte <- paste0("time_to_event_", horizon, "y")
  .outcome <- paste0("outcome_", horizon, "y")

  brier_formula <- as.formula(paste0("Hist(", .tte, ",", .outcome, ") ~ 1"))

  # Original Brier score calculation
  brier <- riskRegression::Score(
    list("with_egfr" = data[[.pred_risks]]),
    formula = brier_formula,
    data = data,
    times = horizon * 365.25,
    cause = 1,
    se.fit = 1L,
    metrics = "brier",
    summary = "ipa",
    contrast = FALSE
  )

  # Extract values
  brier_est <- brier$Brier$score[model == "with_egfr", Brier]
  brier_upper <- brier$Brier$score[model == "with_egfr", upper]
  brier_lower <- brier$Brier$score[model == "with_egfr", lower]

  IPA_est <- brier$IPA$score[model == "with_egfr", IPA]
  IPA_upper <- brier$IPA$score[model == "with_egfr", upper]
  IPA_lower <- brier$IPA$score[model == "with_egfr", lower]

  # Return results

  return(c(
    "brier_est" = brier_est,
    "brier_lower" = brier_lower,
    "brier_upper" = brier_upper,
    "IPA_est" = IPA_est,
    "IPA_lower" = IPA_lower,
    "IPA_upper" = IPA_upper
  ))
}


# Delta Scaled Brier with Bootstrap CI -----------------------------------------

delta_scaled_brier_boot <- function(
  data,
  equation_ref,
  equation_test,
  horizon,
  B = 1000,
  seed = 123
) {
  # Inner function for single calculation (used by boot)
  calc_delta_ipa <- function(
    data,
    indices,
    equation_ref,
    equation_test,
    horizon
  ) {
    # Resample data
    d <- data[indices, , drop = FALSE]

    # Extract variable names
    .pred_risks_ref <- paste0("risk_", horizon, "y_", equation_ref)
    .pred_risks_test <- paste0("risk_", horizon, "y_", equation_test)
    .tte <- paste0("time_to_event_", horizon, "y")
    .outcome <- paste0("outcome_", horizon, "y")

    brier_formula <- as.formula(paste0("Hist(", .tte, ",", .outcome, ") ~ 1"))

    # Reference IPA
    brier_ref <- riskRegression::Score(
      list("ref" = d[[.pred_risks_ref]]),
      formula = brier_formula,
      data = d,
      times = horizon * 365.25,
      cause = 1,
      se.fit = FALSE,
      metrics = "brier",
      summary = "ipa",
      contrast = FALSE
    )
    IPA_ref <- brier_ref$IPA$score[model == "ref", IPA]

    # Test IPA
    brier_test <- riskRegression::Score(
      list("test" = d[[.pred_risks_test]]),
      formula = brier_formula,
      data = d,
      times = horizon * 365.25,
      cause = 1,
      se.fit = FALSE,
      metrics = "brier",
      summary = "ipa",
      contrast = FALSE
    )
    IPA_test <- brier_test$IPA$score[model == "test", IPA]

    return(IPA_test - IPA_ref)
  }

  # Set seed
  set.seed(seed)

  # Run bootstrap
  boot_results <- boot::boot(
    data = data,
    statistic = calc_delta_ipa,
    R = B,
    equation_ref = equation_ref,
    equation_test = equation_test,
    horizon = horizon
  )

  # Get point estimate (original data)
  delta_est <- boot_results$t0

  # Get percentile CI
  boot_ci <- boot::boot.ci(boot_results, type = "perc", conf = 0.95)

  return(tibble(
    delta_scaled_brier = delta_est,
    lower = boot_ci$percent[4],
    upper = boot_ci$percent[5]
  ))
}


################################################################################
# Decision curve analysis ######################################################
################################################################################

compute_dca <- function(
  data,
  equation1,
  equation2,
  horizon,
  max_threshold = 0.2
) {
  risk_eq1 <- paste0("risk_", horizon, "y_", equation1)
  risk_eq2 <- paste0("risk_", horizon, "y_", equation2)
  tte <- paste0("time_to_event_", horizon, "y")
  outcome <- paste0("outcome_", horizon, "y_dca")

  dca_formula <- as.formula(
    paste0(
      "Surv(event=",
      outcome,
      ", time=",
      tte,
      ") ~",
      risk_eq1,
      "+",
      risk_eq2
    )
  )

  dca_tbl <- dca(
    dca_formula,
    data = data,
    time = horizon * 365.25,
    thresholds = seq(0, max_threshold, 0.01)
  ) %>%
    net_intervention_avoided() |>
    as_tibble() |>
    mutate(horizon = horizon, equations = paste0(equation1, "-", equation2))

  return(dca_tbl)
}
