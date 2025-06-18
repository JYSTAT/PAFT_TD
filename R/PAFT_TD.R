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
PAFT_TD <- function(formula, data, id, dist = "lognormal", init = NULL,
                    control = list(maxiter = 2000, tol = 1e-5, algorithm = "NLOPT_LN_COBYLA")) {
  
  call <- match.call()
  
  # ID variable name
  id_var <- deparse(substitute(id))
  
  # Extract formula variables
  all_vars <- all.vars(formula)
  time_var <- all_vars[2]
  event_var <- all_vars[3]
  covariates <- if (length(all_vars) > 3) all_vars[4:length(all_vars)] else NULL
  
  # Last row per subject
  data_last <- data %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[time_var]]) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(ID = all_of(id_var), Time = all_of(time_var), Event = all_of(event_var), all_of(covariates))
  
  # Time-dependency check
  is_time_dependent <- if (is.null(covariates)) FALSE else
    nrow(distinct(data[, c(id_var, covariates)])) != nrow(data_last)
  
  # Case 1: Time-independent
  if (!is_time_dependent) {
    aft_formula <- if (is.null(covariates)) Surv(Time, Event) ~ 1 else
      as.formula(paste("Surv(Time, Event) ~", paste(covariates, collapse = "+")))
    fit <- survreg(aft_formula, data = data_last, dist = dist)
    est <- c(fit$coefficients, Scale = fit$scale)
    names(est)[length(est)] <- "Scale"
    return(list(
      coefficients = round(est, 4),
      logLik = fit$loglik[2],
      converged = TRUE,
      time.dependent = FALSE,
      call = call
    ))
  }
  
  # Case 2: Time-dependent
  if (is.null(init)) {
    init_formula <- as.formula(paste("Surv(Time, Event) ~", paste(covariates, collapse = "+")))
    init_fit <- survreg(init_formula, data = data_last, dist = dist)
    init <- list(beta = init_fit$coefficients, scale = init_fit$scale)
  }
  
  # Negative log-likelihood
  neg_log_likelihood <- function(params) {
    dat <- data[, c(id_var, all_vars)]
    names(dat) <- c("ID", all_vars)
    
    obs_time <- data_last[, c("ID", "Time")]
    names(obs_time)[2] <- "ObsTime"
    dat <- merge(dat, obs_time, by = "ID")
    dat <- dat[order(dat$ID, dat$start), ]
    
    Z_mat <- as.matrix(dat[, covariates])
    beta_vec <- -params[2:(1 + length(covariates))]
    psi_vec <- (dat$stop - dat$start) * exp(Z_mat %*% beta_vec)
    psi <- dat %>%
      mutate(psi = psi_vec) %>%
      group_by(ID) %>%
      summarise(psi = sum(psi), .groups = "drop")
    
    Z_event <- as.matrix(data_last[, covariates])
    psi_event <- exp(Z_event %*% beta_vec)
    
    merged <- merge(psi, data_last[, c("ID", "Event")], by = "ID")
    sigma <- params[length(params)]
    logPsi <- log(merged$psi)
    logT <- (logPsi - params[1]) / sigma
    
    f_event <- dnorm(logT) * psi_event / (sigma * merged$psi)
    S_censor <- 1 - pnorm(logT)
    likelihood <- merged$Event * f_event + (1 - merged$Event) * S_censor
    likelihood <- pmax(likelihood, 1e-20)
    
    return(-sum(log(likelihood)))
  }
  
  # Optimization
  theta_init <- c(init$beta, init$scale)
  bound_val <- max(abs(theta_init), na.rm = TRUE) * 10
  lb <- c(rep(-bound_val, length(theta_init) - 1), 1e-5)
  ub <- rep(bound_val, length(theta_init))
  
  opt_result <- nloptr(
    x0 = theta_init,
    eval_f = neg_log_likelihood,
    lb = lb, ub = ub,
    opts = list(
      algorithm = control$algorithm,
      xtol_rel = control$tol,
      maxeval = control$maxiter
    )
  )
  
  names(opt_result$solution) <- c(paste0("Beta", 0:(length(opt_result$solution) - 2)), "Scale")
  
  return(list(
    coefficients = round(opt_result$solution, 4),
    logLik = -opt_result$objective,
    converged = opt_result$status > 0,
    status = opt_result$status,
    message = opt_result$message,
    n.iter = opt_result$iterations,
    time.dependent = TRUE,
    call = call
  ))
}
