################################################################################
# FILENAME: Functions eGFR equations.R
# PROJECT: Performance and Clinical Utility of the KFRE Using EKFC Versus CKD-EPI eGFR
# PURPOSE: Compute estimated glomerular filtration rate
# AUTHOR: Malou Magnani
# CREATED: 2025-02
# UPDATED:
# - 2025-10-22 (Antoine Creon)

# R VERSION: 4.5
################################################################################

################################################################################
# 1. Creatinine equations ----
################################################################################

# 1.1 CKD-EPIcr2009 ----
ckd_epi_2009_cr <- function(creatinine, age, female) {
  k <- ifelse(female == 1, 62, 80) # for umol/L
  alpha <- ifelse(female == 1, -0.329, -0.411)
  return(ifelse(
    female == 1,
    141 *
      (pmin(creatinine / k, 1)^alpha) *
      (pmax(creatinine / k, 1)^(-1.209)) *
      (0.9929^age) *
      1.018,
    141 *
      (pmin(creatinine / k, 1)^alpha) *
      (pmax(creatinine / k, 1)^(-1.209)) *
      (0.9929^age)
  ))
}

# 1.3 EKFCcr ----
ekfc_cr <- function(creatinine, age, female) {
  Q <- ifelse(
    age <= 25 & female == 0,
    exp(
      3.200 +
        0.259 * age -
        0.543 * log(age) -
        0.00763 * age^2 +
        0.0000790 * age^3
    ),
    ifelse(
      age <= 25 & female == 1,
      exp(
        3.080 +
          0.177 * age -
          0.223 * log(age) -
          0.00596 * age^2 +
          0.0000686 * age^3
      ),
      ifelse(age > 25 & female == 0, 80, ifelse(age > 25 & female == 1, 62, NA))
    )
  )

  return(ifelse(
    creatinine / Q < 1 & age <= 40,
    107.3 * (creatinine / Q)^-0.322,
    ifelse(
      creatinine / Q < 1 & age > 40,
      107.3 * (creatinine / Q)^-0.322 * 0.990^(age - 40),
      ifelse(
        creatinine / Q >= 1 & age <= 40,
        107.3 * (creatinine / Q)^-1.132,
        ifelse(
          creatinine / Q >= 1 & age > 40,
          107.3 * (creatinine / Q)^-1.132 * 0.990^(age - 40),
          NA
        )
      )
    )
  ))
}


################################################################################
#2. Cystatin C equations ----
################################################################################

#2.1 CKD-EPIcys2012 ----
ckd_epi_2012_cys <- function(cystatin, age, female) {
  return(ifelse(
    female == 1,
    133 *
      (pmin(cystatin / 0.8, 1)^(-0.499)) *
      (pmax(cystatin / 0.8, 1)^(-1.328)) *
      (0.9962^age) *
      0.932,
    133 *
      (pmin(cystatin / 0.8, 1)^(-0.499)) *
      (pmax(cystatin / 0.8, 1)^(-1.328)) *
      (0.9962^age)
  ))
}

#2.2 EKFCcys ----
ekfc_cys <- function(cystatin, age) {
  Q <- ifelse(age > 50, 0.83 + 0.005 * (age - 50), 0.83)

  return(ifelse(
    age >= 18 & age <= 40 & cystatin / Q < 1,
    107.3 * (cystatin / Q)^-0.322,
    ifelse(
      age >= 18 & age <= 40 & cystatin / Q >= 1,
      107.3 * (cystatin / Q)^-1.132,
      ifelse(
        age > 40 & age <= 50 & cystatin / Q < 1,
        107.3 * (cystatin / Q)^-0.322 * 0.990^(age - 40),
        ifelse(
          age > 40 & age <= 50 & cystatin / Q >= 1,
          107.3 * (cystatin / Q)^-1.132 * 0.990^(age - 40),
          ifelse(
            age > 50 & cystatin / Q < 1,
            107.3 * (cystatin / Q)^-0.322 * 0.990^(age - 40),
            ifelse(
              age > 50 & cystatin / Q >= 1,
              107.3 * (cystatin / Q)^-1.132 * 0.990^(age - 40),
              NA
            )
          )
        )
      )
    )
  ))
}

################################################################################
#3. Creatinine/cystatin C equations
################################################################################

#3.1 CKD-EPI2012crcys ----
ckd_epi_2012_cr_cys <- function(creatinine, cystatin, age, female) {
  k <- ifelse(female == 1, 62, 80)
  alpha <- ifelse(female == 1, -0.248, -0.207)
  return(ifelse(
    female == 1,
    135 *
      (pmin(creatinine / k, 1)^alpha) *
      (pmax(creatinine / k, 1)^(-0.601)) *
      (pmin(cystatin / 0.8, 1)^(-0.375)) *
      (pmax(cystatin / 0.8, 1)^(-0.711)) *
      (0.9952^age) *
      0.969,
    135 *
      (pmin(creatinine / k, 1)^alpha) *
      (pmax(creatinine / k, 1)^(-0.601)) *
      (pmin(cystatin / 0.8, 1)^(-0.375)) *
      (pmax(cystatin / 0.8, 1)^(-0.711)) *
      (0.9952^age)
  ))
}

#3.3 EKFCcrcys
ekfc_cr_cys <- function(creatinine, cystatin, age, female) {
  cr_value <- ekfc_cr(creatinine, age, female)
  cys_value <- ekfc_cys(cystatin, age)

  return(mean(c(cr_value, cys_value), na.rm = TRUE))
}
