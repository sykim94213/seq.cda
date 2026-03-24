
#' Causal decomposition under fixed intervention targets
#'
#' This function implements a sequential causal decomposition analysis for two
#' time-ordered intervention targets. It estimates the initial disparity, the
#' disparity remaining, and the disparity reduction under a hypothetical policy
#' that fixes both targets to specified values for all individuals, regardless
#' of group membership.

#' @param med_covariates Optional character vector of post-treatment mediator
#'   covariates to include in the outcome model. Use `NULL` if no mediator
#'   covariates are included.
#' @param base_covariates Character vector of baseline covariate names measured
#'   before the first intervention target.
#' @param group Character string giving the name of the group variable in
#'   `data`.
#' @param confounders1 Character vector of confounder names measured before the
#'   first intervention target.
#' @param confounders2 Character vector of confounder names measured after the
#'   first intervention target but before the second intervention target.
#' @param mediator1 Character string giving the name of the first intervention
#'   target in `data`.
#' @param mediator2 Character string giving the name of the second intervention
#'   target in `data`.
#' @param outcome Character string giving the name of the outcome variable in
#'   `data`.
#' @param data A data frame containing the variables used in the analysis.
#' @param comparison Character string giving the comparison group level of
#'   `group`.
#' @param reference Character string giving the reference group level of
#'   `group`.
#' @param mediator1_fixed Value to which the first intervention target is fixed
#'   under the hypothetical policy.
#' @param mediator2_fixed Value to which the second intervention target is fixed
#'   under the hypothetical policy.
#' @param cluster_id Character string giving the cluster identifier used for the
#'   clustered bootstrap.
#' @param n_bootstrap Integer giving the number of bootstrap resamples used to
#'   compute standard errors. Default is `1000`.
#' @param seed Integer random seed for reproducibility. Default is `123`.
#' @param verbose Logical; if `TRUE`, prints progress messages during estimation.
#'
#' @author
#'   Soojin Park, University of California, Riverside, \email{soojinp@@ucr.edu}.
#'
#' @references
#'    Park, S., Kim, S. Y., Zheng, X., & Lee, C. (2025).
#'    Causal decomposition analysis with synergistic interventions:
#'    A triply robust machine-learning approach to addressing multiple dimensions
#'    of social disparities. Psychological Methods. Advance online publication.
#'    https://doi.org/10.1037/met0000803
#'
#' @returns
#' An object of class `"cda_glm_fxfx"` containing:
#' \describe{
#'   \item{tau_glm}{Estimated initial disparity, computed as the bootstrap mean
#'   of the total disparity across replicates.}
#'   \item{delta_glm}{Estimated disparity reduction under the policy that fixes
#'   both intervention targets to `mediator1_fixed` and `mediator2_fixed`,
#'   computed as the bootstrap mean across replicates.}
#'   \item{zeta_glm}{Estimated disparity remaining under the fixed-target policy,
#'   computed as the bootstrap mean across replicates.}
#'   \item{tau_glm_se}{Bootstrap standard error for `tau_glm`.}
#'   \item{delta_glm_se}{Bootstrap standard error for `delta_glm`.}
#'   \item{zeta_glm_se}{Bootstrap standard error for `zeta_glm`.}
#'   \item{meta}{A list of metadata including the number of bootstrap
#'   replicates, random seed, group variable name, comparison group,
#'   reference group, and sample size.}
#' }
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' cda_glm_fxfx(med_covariates = NULL, base_covariates = c("C1", "C2"),
#'   group = "group", confounders1 = c("L1", "L2"), confounders2 = c("L3", "L4"),
#'   mediator1 = "M1", mediator2 = "M2", outcome = "Y", data = dat,
#'   comparison = "Group1", reference = "Group0",
#'   mediator1_fixed = 1, mediator2_fixed = 1,
#'   cluster_id = "cluster", n_bootstrap = 500, seed = 123, verbose = TRUE)
#' }

cda_glm_fxfx <- function(med_covariates = NULL, base_covariates, group, confounders1, confounders2,
                         mediator1, mediator2, outcome, data, comparison, reference, mediator1_fixed, mediator2_fixed,
                         cluster_id, n_bootstrap = 1000, seed = 123, verbose = TRUE) {
  # ---- Setup and reproducibility ----
  set.seed(seed)
  if (verbose) message("Starting cda_glm_fxfx()")

  # ---- Input validation ----
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  if (!is.character(group) || length(group) != 1) {
    stop("'group' must be a single column name.")
  }
  if (!is.character(cluster_id) || length(cluster_id) != 1) {
    stop("'cluster_id' must be a single column name.")
  }
  if (!is.numeric(n_bootstrap) || length(n_bootstrap) != 1 || is.na(n_bootstrap) || n_bootstrap < 1) {
    stop("'n_bootstrap' must be a single number >= 1.")
  }
  n_bootstrap <- as.integer(n_bootstrap)
  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
    stop("'seed' must be a single numeric value.")
  }
  needed_cols <- unique(c(
    base_covariates, group, confounders1, confounders2,
    mediator1, mediator2, outcome, cluster_id
  ))
  missing_cols <- setdiff(needed_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("The following required columns are missing in 'data': %s",
              paste(missing_cols, collapse = ", "))
    )
  }
  group_values <- unique(as.character(data[[group]]))
  if (!(comparison %in% group_values)) {
    stop(sprintf("'comparison' value '%s' is not present in data[['%s']].", comparison, group))
  }
  if (!(reference %in% group_values)) {
    stop(sprintf("'reference' value '%s' is not present in data[['%s']].", reference, group))
  }
  if (identical(comparison, reference)) {
    stop("'comparison' and 'reference' must be different values.")
  }
  m1_fixed_num <- suppressWarnings(as.integer(as.character(mediator1_fixed)))
  m2_fixed_num <- suppressWarnings(as.integer(as.character(mediator2_fixed)))
  if (is.na(m1_fixed_num) || !(m1_fixed_num %in% c(0L, 1L))) {
    stop("'mediator1_fixed' must be binary (0 or 1).")
  }
  if (is.na(m2_fixed_num) || !(m2_fixed_num %in% c(0L, 1L))) {
    stop("'mediator2_fixed' must be binary (0 or 1).")
  }
  mediator1_vals <- suppressWarnings(as.integer(as.character(stats::na.omit(data[[mediator1]]))))
  mediator2_vals <- suppressWarnings(as.integer(as.character(stats::na.omit(data[[mediator2]]))))
  if (length(mediator1_vals) == 0 || any(is.na(mediator1_vals)) || !all(mediator1_vals %in% c(0L, 1L))) {
    stop(sprintf("'%s' must be binary and coercible to 0/1.", mediator1))
  }
  if (length(mediator2_vals) == 0 || any(is.na(mediator2_vals)) || !all(mediator2_vals %in% c(0L, 1L))) {
    stop(sprintf("'%s' must be binary and coercible to 0/1.", mediator2))
  }

  # ---- Bootstrap replicate estimator ----
  disparity_calculation_glm <- function(data_boot) {
    # ---- (A) Covariate centering for marginal contrast interpretation ----
    # Centering over comparison/reference rows makes the group effect interpretable
    # as a contrast marginalized over baseline covariates.
    center_rows <- data_boot %>% dplyr::filter(!!rlang::sym(group) %in% c(comparison, reference))
    for (cov in base_covariates) {
      x <- center_rows[[cov]]
      if (is.factor(x)) x <- as.numeric(as.character(x))
      x_mean <- mean(x, na.rm = TRUE)
      # Coerce factor-coded numeric covariates before centering so linear-model
      # terms use numeric offsets rather than factor contrasts.
      if (is.factor(data_boot[[cov]])) {
        data_boot[[cov]] <- as.numeric(as.character(data_boot[[cov]]))
      }
      data_boot[[cov]] <- data_boot[[cov]] - x_mean
    }

    data_r1 <- data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison)
    data_r0 <- data_boot %>% dplyr::filter(!!rlang::sym(group) == reference)

    # ---- (B) Outcome model and nested imputation ----
    # Outcome regression supports nested imputation under fixed mediator policy.
    fit.y <- lm(y_form, data = data_boot)

    # Apply policy-fixed mediator values to create model-implied outcomes.
    dat_1 <- data_boot
    dat_1[[mediator1]] <- factor(mediator1_fixed, levels = levels(data_boot[[mediator1]]))
    dat_1[[mediator2]] <- factor(mediator2_fixed, levels = levels(data_boot[[mediator2]]))
    data_boot$muldm <- predict(fit.y, newdata = dat_1)

    # Nested imputation:
    # - Inner: E[Y | R, X, A=fixed, Z, M=fixed, C] from `muldm`
    # - Outer: E[Inner | R, X, A, C] from regression of `muldm`
    # - Marginalization over centered baseline covariates via model intercept.
    #
    # Comparison arm: outer model intentionally excludes confounders2 to match
    # the target decomposition definition used in this implementation.
    fit_final <- lm(
      as.formula(paste0(
        "muldm ~ ",
        paste(c(mediator1, base_covariates, confounders1), collapse = " + ")
      )),
      data = data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison)
    )
    data_boot$nu <- predict(fit_final, newdata = dat_1)
    b <- lm(
      as.formula(paste0("nu ~ ", paste(base_covariates, collapse = " + "))),
      data = data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison)
    )

    # Reference arm mirror of nested-imputation anchor for a comparable contrast.
    fit_final_ref <- lm(
      as.formula(paste0(
        "muldm ~ ",
        paste(c(mediator1, base_covariates, confounders1), collapse = " + ")
      )),
      data = data_boot %>% dplyr::filter(!!rlang::sym(group) == reference)
    )
    data_boot$nu_ref <- predict(fit_final_ref, newdata = dat_1)
    b_ref <- lm(
      as.formula(paste0("nu_ref ~ ", paste(base_covariates, collapse = " + "))),
      data = data_boot %>% dplyr::filter(!!rlang::sym(group) == reference)
    )

    # ---- (C) Weight-model fitting (comparison + reference) ----
    # Fit propensity models used to transport observed data to policy-fixed worlds.
    fit.m1 <- glm(m_form, data = data_r1, family = binomial(link = "logit"))
    fit.a1 <- glm(a_form, data = data_r1, family = binomial(link = "logit"))
    fit.m0_wgt <- glm(m_form, data = data_r0, family = binomial(link = "logit"))
    fit.a0_wgt <- glm(a_form, data = data_r0, family = binomial(link = "logit"))

    pm1 <- predict(fit.m1, data_boot, type = "response")
    pa1 <- predict(fit.a1, data_boot, type = "response")
    pm0 <- predict(fit.m0_wgt, data_boot, type = "response")
    pa0 <- predict(fit.a0_wgt, data_boot, type = "response")

    # ---- (D) Weight construction and percentile trimming ----
    # Build inverse-probability style weights for observed rows consistent with
    # fixed (mediator1, mediator2) intervention values.
    #
    # Assumption: mediator columns are binary and encoded as numeric-like values
    # (e.g., 0/1 or "0"/"1"). If factor levels are non-numeric labels, coercion
    # via as.integer(as.character(.)) can produce NA and invalidate weights.
    A_obs <- as.integer(as.character(data_boot[[mediator1]]))
    M_obs <- as.integer(as.character(data_boot[[mediator2]]))

    pA_fixed <- if (mediator1_fixed == 1) pa1 else (1 - pa1)
    pM_fixed <- if (mediator2_fixed == 1) pm1 else (1 - pm1)

    WM <- ifelse(
      data_boot[[group]] == comparison & A_obs == mediator1_fixed & M_obs == mediator2_fixed,
      1 / pM_fixed,
      0
    )

    WA <- ifelse(
      data_boot[[group]] == comparison & A_obs == mediator1_fixed,
      1 / pA_fixed,
      0
    )

    # Combined weight under fixed-both intervention.
    WAM <- WA * WM

    # Trimming stabilizes estimation when a few units receive extreme weights.
    wpos <- WAM[WAM > 0]
    if (length(wpos) == 0) {
      stop("No positive comparison-group intervention weights were generated.")
    }
    lower <- quantile(wpos, 0.01, na.rm = TRUE)
    upper <- quantile(wpos, 0.99, na.rm = TRUE)
    WAM[WAM > upper] <- upper
    WAM[WAM < lower & WAM > 0] <- lower
    data_boot$WAM <- WAM

    # Construct reference-group analog weights so comparison and reference terms
    # are computed under aligned intervention rules.
    pA_fixed_ref <- if (mediator1_fixed == 1) pa0 else (1 - pa0)
    pM_fixed_ref <- if (mediator2_fixed == 1) pm0 else (1 - pm0)
    WA_ref <- ifelse(
      data_boot[[group]] == reference & A_obs == mediator1_fixed,
      1 / pA_fixed_ref,
      0
    )
    WM_ref <- ifelse(
      data_boot[[group]] == reference & A_obs == mediator1_fixed & M_obs == mediator2_fixed,
      1 / pM_fixed_ref,
      0
    )
    WAM_ref <- WA_ref * WM_ref
    wpos_ref <- WAM_ref[WAM_ref > 0]
    if (length(wpos_ref) > 0) {
      lower_ref <- quantile(wpos_ref, 0.01, na.rm = TRUE)
      upper_ref <- quantile(wpos_ref, 0.99, na.rm = TRUE)
      WAM_ref[WAM_ref > upper_ref] <- upper_ref
      WAM_ref[WAM_ref < lower_ref & WAM_ref > 0] <- lower_ref
    }


    # ---- (E) Component estimators and triply robust contrast ----
    # Weighted observed outcome component.
    wmu1xdm <- lm(as.formula(paste0(outcome, " ~ ", paste(base_covariates, collapse = " + "))), weights = WAM, # WA * WM
                  data = data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison))

    # Weighted imputation component.
    wx_mu <- lm(as.formula(paste0("muldm", " ~ ", paste(base_covariates, collapse = " + "))), weights = WAM, # WA
                data = data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison))

    # Weighted nested-imputation component for triply robust assembly.
    wx_nu <- lm(as.formula(paste0("nu", " ~ ", paste(base_covariates, collapse = " + "))), weights = WAM, # WA
                data = data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison))
    data_boot$mu <- predict(fit.y, newdata = data_boot, type = "response")
    wxwm_mu <- lm(as.formula(paste0("mu", " ~ ", paste(base_covariates, collapse = " + "))), weights = WAM, # WA * WM
                  data = data_boot %>% dplyr::filter(!!rlang::sym(group) == comparison))

    # Reference-group mirror terms ensure the final contrast is on the same
    # estimand scale before subtraction.
    data_ref_tr <- data_boot %>% dplyr::filter(!!rlang::sym(group) == reference)
    WAM_ref_tr <- WAM_ref[data_boot[[group]] == reference]

    wx_mu_ref_tr <- lm(
      as.formula(paste0("muldm", " ~ ", paste(base_covariates, collapse = " + "))),
      weights = WAM_ref_tr,
      data = data_ref_tr
    )
    wx_nu_ref_tr <- lm(
      as.formula(paste0("nu_ref", " ~ ", paste(base_covariates, collapse = " + "))),
      weights = WAM_ref_tr,
      data = data_ref_tr
    )
    wmu1xdm_ref_tr <- lm(
      as.formula(paste0(outcome, " ~ ", paste(base_covariates, collapse = " + "))),
      weights = WAM_ref_tr,
      data = data_ref_tr
    )
    wxwm_mu_ref_tr <- lm(
      as.formula(paste0("mu", " ~ ", paste(base_covariates, collapse = " + "))),
      weights = WAM_ref_tr,
      data = data_ref_tr
    )
    term_comp <- b$coef[1] + (coef(wx_mu)[1] - coef(wx_nu)[1]) + (wmu1xdm$coef[1] - coef(wxwm_mu)[1])
    term_ref  <- b_ref$coef[1] + (coef(wx_mu_ref_tr)[1] - coef(wx_nu_ref_tr)[1]) +
      (wmu1xdm_ref_tr$coef[1] - coef(wxwm_mu_ref_tr)[1])
    zeta_tr <- term_comp - term_ref

    # ---- (F) Total disparity and residual ----
    # Total disparity with centered covariates equals the coefficient on the
    # binary group indicator in Y ~ R_ind_tau + base_covariates.
    tau_data <- data_boot %>% dplyr::filter(!!rlang::sym(group) %in% c(comparison, reference))
    tau_data$R_ind_tau <- ifelse(tau_data[[group]] == comparison, 1, 0)
    tau_fit <- lm(
      as.formula(paste0(outcome, " ~ R_ind_tau + ", paste(base_covariates, collapse = " + "))),
      data = tau_data
    )
    tau_glm <- coef(tau_fit)[["R_ind_tau"]]

    delta_tr <- tau_glm - zeta_tr

    # Return one replicate: (tau, delta, zeta)
    return(c(tau_glm, delta_tr, zeta_tr))
  }

  # ---- Formula construction ----
  if (verbose) message("Building model formulas...")
  # A ~ X + C
  a_rhs <- c(base_covariates, confounders1) %>%
    purrr::map(~ rlang::sym(.x)) %>%
    purrr::reduce(~ rlang::expr(!!.x + !!.y))
  a_form <- as.formula(rlang::expr(!!rlang::sym(mediator1) ~ !!a_rhs))

  # M ~ X + A + Z + C
  m_rhs <- c(base_covariates, confounders1, mediator1, confounders2) %>%
    purrr::map(~ rlang::sym(.x)) %>%
    purrr::reduce(~ rlang::expr(!!.x + !!.y))
  m_form <- as.formula(rlang::expr(!!rlang::sym(mediator2) ~ !!m_rhs))

  # Y ~ R + X + A + Z + M + C
  interaction_terms <- rlang::exprs(!!rlang::sym(group) * !!rlang::sym(mediator2), !!rlang::sym(group) * !!rlang::sym(mediator1))
  y_rhs <- c(base_covariates, group, confounders1, mediator1, confounders2, mediator2) %>%
    purrr::map(~ rlang::sym(.x)) %>%
    purrr::reduce(~ rlang::expr(!!.x + !!.y))
  y_rhs <- rlang::expr(!!y_rhs + !!(interaction_terms[[1]]) + !!(interaction_terms[[2]]))
  y_form <- as.formula(rlang::expr(!!rlang::sym(outcome) ~ !!y_rhs))

  # ---- Clustered bootstrap driver ----
  if (verbose) message(sprintf("Running clustered bootstrap (%d replicates)...", n_bootstrap))
  clusters <- unique(data[[cluster_id]])
  if (length(clusters) == 0) {
    stop(sprintf("No non-missing clusters found in '%s'.", cluster_id))
  }
  glm_results <- replicate(n_bootstrap, {
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    boot_data <- data %>% dplyr::filter(!!rlang::sym(cluster_id) %in% sampled_clusters)
    if (nrow(boot_data) == 0) {
      stop("A bootstrap replicate produced zero rows after cluster sampling.")
    }
    disparity_calculation_glm(boot_data)
  }, simplify = TRUE)

  # ---- Aggregate bootstrap means and standard errors ----
  # Means across replicates.
  tau_glm <- mean(glm_results[1, ])
  delta_glm <- mean(glm_results[2, ])
  zeta_glm <- mean(glm_results[3, ])

  # Bootstrap SDs as standard errors.
  tau_glm_se <- sd(glm_results[1, ])
  delta_glm_se <- sd(glm_results[2, ])
  zeta_glm_se <- sd(glm_results[3, ])

  # Return package-facing summary output:
  # - tau_glm, delta_glm, zeta_glm: bootstrap means of decomposition components
  # - *_se: bootstrap standard deviations (uncertainty estimates)
  out <- list(
    tau_glm = tau_glm, delta_glm = delta_glm, zeta_glm = zeta_glm,
    tau_glm_se = tau_glm_se, delta_glm_se = delta_glm_se, zeta_glm_se = zeta_glm_se,
    meta = list(
      n_bootstrap = n_bootstrap,
      seed = seed,
      group = group,
      comparison = comparison,
      reference = reference,
      n_obs = nrow(data)
    )
  )
  class(out) <- "cda_glm_fxfx"
  if (verbose) message("Finished cda_glm_fxfx().")
  return(out)
}
