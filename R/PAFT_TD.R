#' Fit a parametric AFT model with time-dependent covariates
#'
#' @description
#' This function estimates a parametric AFT (accelerated failure time) model that accommodates
#' both time-independent and time-dependent covariates using a log-likelihood-based optimization.
#'
#' @param formula A formula object with \code{Surv(start, stop, event)} on the left-hand side.
#' @param data A data frame in counting process format, containing start/stop times and covariates.
#' @param id A column name indicating the subject identifier.
#' @param dist The assumed distribution for the survival time (default is "lognormal").
#' @param init Optional list containing \code{beta} and \code{scale} as initial values for optimization.
#' @param control A list of optimization controls: \code{maxiter}, \code{tol}, \code{algorithm}.
#'
#' @return A list containing:
#' \item{coefficients}{Estimated regression coefficients and scale parameter}
#' \item{logLik}{Log-likelihood at the optimum}
#' \item{converged}{Logical indicator of convergence}
#' \item{status}{NLOPT status code}
#' \item{message}{Message returned from the optimizer}
#' \item{n.iter}{Number of iterations performed}
#' \item{time.dependent}{Logical indicating if TD covariates were used}
#' \item{call}{Matched function call}
#'
#' @export
#' @import dplyr survival nloptr
PAFT_TD <- function(formula, data, id_var = "ID", 
                    dist = "lognormal",                # flexible distribution
                    initi = FALSE, beta = NA, sigma_init = NA,
                    tol = 1.0e-5, maxiter = 2000,
                    algorithm = "NLOPT_LN_COBYLA") {

  # --- Distribution-specific log-likelihood component --- #
  get_loglik_function <- function(dist) {
    switch(dist,
           "lognormal" = function(z, psi, psi_d, delta, sigma) {
             f <- dnorm(z) * (psi_d / (sigma * psi))
             S <- 1 - pnorm(z)
             delta * f + (1 - delta) * S
           },
           stop("Unsupported distribution: ", dist))
  }
  loglik_fn <- get_loglik_function(dist)

  # --- Parse survival formula --- #
  mf <- model.frame(formula, data)
  surv_obj <- model.response(mf)
  if (!inherits(surv_obj, "Surv") || attr(surv_obj, "type") != "counting") {
    stop("Formula must be of the form Surv(start, stop, event) ~ covariates")
  }

  if (!(id_var %in% names(data))) {
    stop(paste0("The data must contain the ID variable '", id_var, "'"))
  }

  # --- Check if time-independent --- #
  id_count <- dplyr::count(data, !!rlang::sym(id_var))
  is_time_indep <- all(id_count$n == 1)

  if (is_time_indep) {
    data$start <- 0
    data$stop  <- surv_obj[, 2]
    data$delta <- surv_obj[, 3]
    formula_ti <- as.formula(
      paste("Surv(stop, delta) ~", paste(attr(terms(formula), "term.labels"), collapse = " + "))
    )
    fit <- survreg(formula_ti, data = data, dist = dist)
    est_params <- c(fit$coefficients, sigma = fit$scale)
    names(est_params)[length(est_params)] <- "sigma"
    return(list(EstBetas = est_params,
                method = "survreg",
                message = paste("Time-independent data fitted using survreg() with dist =", dist)))
  }

  # --- Prepare design matrix --- #
  start_time <- surv_obj[, 1]
  stop_time  <- surv_obj[, 2]
  event      <- surv_obj[, 3]
  X_mat_full <- model.matrix(formula, data)
  X_mat <- X_mat_full[, -1, drop = FALSE]
  covariate_names <- colnames(X_mat)

  data$start <- start_time
  data$stop  <- stop_time
  data$delta <- event
  for (j in seq_along(covariate_names)) {
    data[[covariate_names[j]]] <- X_mat[, j]
  }

  # --- Initial values --- #
  if (initi) {
    init_beta <- beta
    init_sigma <- sigma_init
  } else {
    data_last <- data %>% dplyr::group_by(.data[[id_var]]) %>% dplyr::arrange(stop) %>% dplyr::slice(n()) %>% dplyr::ungroup()
    X_last <- as.matrix(data_last[, covariate_names])
    fit_init <- survreg(Surv(stop, delta) ~ X_last, data = data_last, dist = dist)
    init_beta <- fit_init$coefficients
    init_sigma <- fit_init$scale
  }

  # --- Negative log-likelihood --- #
  neg_log_likelihood <- function(params) {
    intercept <- params[1]
    betas <- params[2:(length(covariate_names) + 1)]
    sigma <- params[length(params)]

    obs_time <- data %>% dplyr::group_by(.data[[id_var]]) %>% dplyr::slice(n()) %>% dplyr::ungroup() %>%
      dplyr::select(!!rlang::sym(id_var), stop) %>% dplyr::rename(obs_time = stop)
    dat_ext <- data %>% dplyr::left_join(obs_time, by = id_var) %>% dplyr::arrange(.data[[id_var]], start)

    Z_mat <- as.matrix(dat_ext[, covariate_names])
    lin_pred <- -Z_mat %*% betas
    dat_ext$psi_comp <- (dat_ext$stop - dat_ext$start) * exp(lin_pred)

    psi_df <- dat_ext %>% dplyr::group_by(.data[[id_var]]) %>%
      dplyr::summarise(psi = sum(psi_comp), .groups = "drop") %>% dplyr::arrange(.data[[id_var]])

    Z_last <- dat_ext %>% dplyr::group_by(.data[[id_var]]) %>% dplyr::slice(n()) %>%
      dplyr::ungroup() %>% dplyr::select(all_of(covariate_names)) %>% as.matrix()
    psi_d <- exp(-Z_last %*% betas)

    delta_df <- dat_ext %>% dplyr::group_by(.data[[id_var]]) %>% dplyr::slice(n()) %>%
      dplyr::ungroup() %>% dplyr::select(!!rlang::sym(id_var), delta)

    final_df <- psi_df %>% dplyr::mutate(psi_d = as.vector(psi_d)) %>%
      dplyr::left_join(delta_df, by = id_var) %>% dplyr::arrange(.data[[id_var]])

    log_psi <- log(final_df$psi)
    z <- (log_psi - intercept) / sigma
    L <- loglik_fn(z, final_df$psi, final_df$psi_d, final_df$delta, sigma)
    -sum(log(pmax(L, 1.0e-20)))
  }

  # --- Optimize using nloptr --- #
  init_params <- c(init_beta, init_sigma)
  bound_size <- max(abs(init_params)) * 10
  lower_bounds <- rep(-bound_size, length(init_params))
  upper_bounds <- rep(bound_size, length(init_params))

  res <- nloptr::nloptr(x0 = init_params,
                eval_f = neg_log_likelihood,
                lb = lower_bounds,
                ub = upper_bounds,
                opts = list(algorithm = algorithm,
                            xtol_rel = tol,
                            maxeval = maxiter))

  est_params <- res$solution
  names(est_params) <- c("Intercept", covariate_names, "sigma")

  return(list(EstBetas = est_params,
              method = "PAFT_TD",
              itr = res$iterations,
              status = res$status,
              message = res$message))
}
