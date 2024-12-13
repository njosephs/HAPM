# [ bsglm: Bayesian Inference and Sampling for GLMs ]

EPSILON <- .Machine$double.eps

# [ Helper functions ]
mat_mult <- function (A, b, mult_cond = is.vector(A))
  if (mult_cond) A * b else A %*% b

mat_add <- function (A, B, add_cond = is.vector(A)) {
  if (add_cond) diag(B) <- diag(B) + A else B <- B + A
  B
}

mat_quad <- function (A, b, quad_cond = is.vector(A))
  if (quad_cond) sum(A * b ^ 2) else crossprod(b, crossprod(A, b))

chol_solve <- function (C, b)
  backsolve(C, backsolve(C, b, transpose = TRUE))


rinvchisq <- function (n, nu, nu_tau2) 1 / rgamma(n, nu / 2, nu_tau2 / 2)

mcmc_init_array <- function (ns, nc, params)
  array(dim = c(ns, nc, length(params)),
        dimnames = list(iterations = NULL,
                        chains = paste0("chain:", seq_len(nc)),
                        parameters = params))


# [ Main functions ]

#' Control parameters for bsglm routines
#'
#' Specify control parameter for bsglm routines: tolerance and maximum number
#' of iterations for fitting, number of chains and samples for MCMC sampling,
#' and an execution trace flag, similar to `glm.control`.
#'
#' @param epsilon convergence tolerance
#' @param maxit maximum number of iterations
#' @param nsamples number of MCMC samples
#' @param nchains number of MCMC chains
#' @param trace log execution trace?
#' @return List with control parameters.
#' @export
control <- function (epsilon = 1e-8, maxit = 25, nsamples = 1000,
                           nchains = 1, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("invalid epsilon")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("invalid max number of iterations")
  if (!is.numeric(nsamples) || nsamples <= 0)
    stop("invalid number of samples")
  if (!is.numeric(nchains) || nchains <= 0)
    stop("invalid number of chains")
  list(epsilon = epsilon, maxit = maxit, nsamples = nsamples,
       nchains = nchains, trace = trace)
}


#' Fit a Bayesian linear model
#'
#' Find the maximum a posteriori (MAP) estimator of Bayesian linear model with
#' coefficients \eqn{\beta}{beta} and dispersion \eqn{\phi}{\phi}.
#' The coefficients have a normal prior with mean
#' \eqn{\beta_0}{beta0} and precision matrix \eqn{\Omega}{Omega},
#' \eqn{\beta \sim N(\beta_0, \Omega^{-1})}{beta ~ N(beta0, Omega^{-1})}, and
#' the dispersion has an inverse scaled chi-square prior with degrees of
#' freedom \eqn{\nu}{nu} and scale \eqn{\tau^2}{tau2},
#' \eqn{1/\phi \sim Ga(\nu/2, \nu\tau^2/2)}{1/phi ~ Ga(nu / 2, nu * tau2 / 2)}.
#'
#' This routine is a `lm` method and should be passed to a `lm` call, e.g.,
#' `lm(formula, data, ..., method = bsglm::lm_fitter(...))`.
#'
#' @param prior_coef hyper-prior parameters for the coefficients, specified as
#' a list with entries `mean`, a vector, and `precision`, either a matrix or a
#' vector storing the diagonal entries
#' @param prior_disp hyper-prior parameters for the dispersion, specified as a
#' list with entries `df`, the number of degrees of freedom, and `scale`, the
#' variance prior scale
#' @return a `lm` object.
#' @export
lm_fitter <- function (prior_coef = NULL, prior_disp = NULL) {
  function (x, y, weights = NULL, start = NULL,
            etastart = NULL, mustart = NULL, offset = NULL,
            family = gaussian(), control = list(),
            intercept = TRUE, singular.ok = TRUE) {
    # initialize
    control <- do.call("control", control)
    y <- drop(y)
    if (!is.null(offset)) y <- y - offset
    nvars <- ncol(x); nobs <- nrow(x)
    xnames <- colnames(x)
    if (is.null(xnames)) xnames <- paste0("x", 1L:nvars)
    if (is.null(weights)) weights <- rep.int(1, nobs)
    deviance <- Inf
    phi <- 1
    if (is.null(prior_coef))
      prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
    if (is.null(prior_disp))
      prior_disp <- list(df = 0, scale = 0)
    Omega <- prior_coef$precision
    beta0 <- drop(prior_coef$mean)
    Omega_beta0 <- mat_mult(Omega, beta0)
    nu <- prior_disp$df
    nu_tau2 <- nu * prior_disp$scale

    # [ fit ]
    for (iter in 1:control$maxit) {
      W <- weights / phi
      z <- crossprod(x, y * W) + Omega_beta0
      C <- chol(mat_add(Omega, crossprod(x, x * W)))
      coef <- drop(chol_solve(C, z))
      fitted_values <- drop(x %*% coef)
      deviance_new <- sum(weights * (y - fitted_values) ^ 2)
      phi <- (nu_tau2 + deviance_new) / (nu + nobs) # E-step: E[1 / sigma^2]

      rel_error <- abs((deviance_new - deviance) / deviance)
      if (!is.infinite(deviance) &&
          ((rel_error < control$epsilon) || (deviance_new < EPSILON)))
        break
      deviance <- deviance_new
    }

    wtdmu <- if (intercept)
      (sum(weights * y) / phi + Omega_beta0[1]) /
        (sum(weights) / phi + Omega[1])
    else
      linkinv(offset + Omega_beta0[1])
    null_deviance <- sum(weights * (y - wtdmu) ^ 2)
    residuals <- sqrt(weights) * (y - fitted_values)
    if (!is.null(offset)) fitted_values <- fitted_values + offset
    names(coef) <- xnames
    n_ok <- nobs - sum(weights == 0)
    df_null <- n_ok - as.integer(intercept)
    df_residual <- n_ok - nvars # rank = nvars
    aic <- nobs * (log(2 * pi * phi) + 1) - sum(log(weights)) + 2 * nvars

    list(coefficients = coef, dispersion = phi, deviance = deviance,
         prior.coef = prior_coef, prior.disp = prior_disp,
         df.null = df_null, df.residual = df_residual,
         null.deviance = null_deviance, aic = aic,
         qr = structure(list(qr = C, pivot = 1L:nvars), class = "qr"),
         residuals = residuals, fitted.values = fitted_values,
         linear.predictors = fitted_values, iter = iter,
         weights = weights, y = y, rank = nvars)
  }
}

#' Sample from a Bayesian linear model
#'
#' Generate MCMC Gibbs samples for a Bayesian linear model with coefficients
#' \eqn{\beta}{beta} and dispersion \eqn{\phi}{\phi}.
#' The coefficients have a normal prior with mean
#' \eqn{\beta_0}{beta0} and precision matrix \eqn{\Omega}{Omega},
#' \eqn{\beta \sim N(\beta_0, \Omega^{-1})}{beta ~ N(beta0, Omega^{-1})}, and
#' the dispersion has an inverse scaled chi-square prior with degrees of
#' freedom \eqn{\nu}{nu} and scale \eqn{\tau^2}{tau2},
#' \eqn{1/\phi \sim Ga(\nu/2, \nu\tau^2/2)}{1/phi ~ Ga(nu / 2, nu * tau2 / 2)}.
#'
#' This routine is a `lm` method and should be passed to a `lm` call, e.g.,
#' `lm(formula, data, ..., method = bsglm::lm_sampler(...))`.
#'
#' @param prior_coef hyper-prior parameters for the coefficients, specified as
#' a list with entries `mean`, a vector, and `precision`, either a matrix or a
#' vector storing the diagonal entries
#' @param prior_disp hyper-prior parameters for the dispersion, specified as a
#' list with entries `df`, the number of degrees of freedom, and `scale`, the
#' variance prior scale
#' @return a `lm` object.
#' @export
lm_sampler <- function (prior_coef = NULL, prior_disp = NULL) {
  function (x, y, weights = NULL, start = NULL,
            etastart = NULL, mustart = NULL, offset = NULL,
            family = gaussian(), control = list(),
            intercept = TRUE, singular.ok = TRUE) {
    # initialize
    control <- do.call("control", control)
    y <- drop(y)
    if (!is.null(offset)) y <- y - offset
    nvars <- ncol(x); nobs <- nrow(x)
    xnames <- colnames(x)
    if (is.null(xnames)) xnames <- paste0("x", 1L:nvars)
    if (is.null(weights)) weights <- rep.int(1, nobs)
    if (is.null(prior_coef))
      prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
    if (is.null(prior_disp))
      prior_disp <- list(df = 0, scale = 0)
    Omega <- prior_coef$precision
    beta0 <- drop(prior_coef$mean)
    Omega_beta0 <- mat_mult(Omega, beta0)
    nu <- prior_disp$df
    nu_tau2 <- nu * prior_disp$scale
    phi <- 1


    # [ sample ]
    sims <- mcmc_init_array(control$nsamples, control$nchains,
                            c(xnames, "dispersion", "__lp"))
    for (chain in 1:control$nchains) {
      for (iter in 1:control$nsamples) {
        # sample coef
        W <- weights / phi
        z <- crossprod(x, y * W) + Omega_beta0
        C <- chol(mat_add(Omega, crossprod(x, x * W)))
        coef <- backsolve(C, z, transpose = TRUE)
        coef <- drop(backsolve(C, coef + rnorm(nvars)))
        fitted_values <- drop(x %*% coef)
        deviance <- sum(weights * (y - fitted_values) ^ 2)
        # sample phi
        phi <- rinvchisq(1, nu + nobs, nu_tau2 + deviance)
        # compute log posterior
        lp <- -.5 * (mat_quad(Omega, coef - beta0) + deviance / phi)

        sims[iter, chain, 1:nvars] <- coef
        sims[iter, chain, nvars + 1] <- phi
        sims[iter, chain, nvars + 2] <- lp
      }
    }

    s_mean <- apply(sims, 3, mean)
    coef <- s_mean[1L:nvars]; phi <- s_mean[nvars + 1]
    # intercept only:
    wtdmu <- if (intercept)
      (sum(weights * y) / phi + Omega_beta0[1]) /
        (sum(weights) / phi + Omega[1])
    else
      linkinv(offset + Omega_beta0[1])
    null_deviance <- sum(weights * (y - wtdmu) ^ 2)
    W <- weights / phi
    C <- chol(mat_add(Omega, crossprod(x, x * W)))
    fitted_values <- drop(x %*% coef)
    residuals <- sqrt(weights) * (y - fitted_values)
    if (!is.null(offset)) fitted_values <- fitted.values + offset
    deviance <- sum(residuals ^ 2)
    names(coef) <- xnames
    n_ok <- nobs - sum(weights == 0)
    df_null <- n_ok - as.integer(intercept)
    df_residual <- n_ok - nvars # rank = nvars
    aic <- nobs * (log(2 * pi * phi) + 1) - sum(log(weights)) + 2 * nvars

    list(coefficients = coef, dispersion = phi, deviance = deviance,
         prior.coef = prior_coef, prior.disp = prior_disp, samples = sims,
         df.null = df_null, df.residual = df_residual,
         null.deviance = null_deviance, aic = aic,
         qr = structure(list(qr = C, pivot = 1L:nvars), class = "qr"),
         residuals = residuals, fitted.values = fitted_values,
         linear.predictors = fitted_values, iter = iter,
         weights = weights, y = y, rank = nvars)
  }
}



#' Fit a Bayesian generalized linear model
#'
#' Find the maximum a posteriori (MAP) estimator of Bayesian generalized
#' linear model with coefficients \eqn{\beta}{beta} and dispersion
#' \eqn{\phi}{\phi}. The coefficients have a normal prior with mean
#' \eqn{\beta_0}{beta0} and precision matrix \eqn{\Omega}{Omega},
#' \eqn{\beta \sim N(\beta_0, \Omega^{-1})}{beta ~ N(beta0, Omega^{-1})}, and
#' the dispersion has an inverse scaled chi-square prior with degrees of
#' freedom \eqn{\nu}{nu} and scale \eqn{\tau^2}{tau2},
#' \eqn{1/\phi \sim Ga(\nu/2, \nu\tau^2/2)}{1/phi ~ Ga(nu / 2, nu * tau2 / 2)}.
#'
#' This routine is a `glm` method and should be passed to a `glm` call, e.g.,
#' `glm(formula, data, family, ..., method = bsglm::fitter(...))`.
#'
#' @param dispersion value of dispersion if assumed deterministic or initial
#' value for fitting
#' @param adjust_dispersion is dispersion random?
#' @param prior_coef hyper-prior parameters for the coefficients, specified as
#' a list with entries `mean`, a vector, and `precision`, either a matrix or a
#' vector storing the diagonal entries
#' @param prior_disp hyper-prior parameters for the dispersion, specified as a
#' list with entries `df`, the number of degrees of freedom, and `scale`, the
#' variance prior scale
#' @return a `glm` object.
#' @export
fitter <- function (dispersion = 1, adjust_dispersion = FALSE,
                    prior_coef = NULL, prior_disp = NULL) {
  function (x, y, weights = NULL, start = NULL,
            etastart = NULL, mustart = NULL, offset = NULL,
            family = gaussian(), control = list(),
            intercept = TRUE, singular.ok = TRUE) {
    control <- do.call("control", control)
    x <- as.matrix(x)
    nobs <- nrow(x)
    nvars <- ncol(x)
    xnames <- colnames(x)
    if (is.null(xnames)) xnames <- paste0("x", 1L:nvars)
    if (nvars == 0L) stop("empty models are not allowed")
    if (is.null(weights)) weights <- rep.int(1, nobs)
    if (is.null(offset)) offset <- rep.int(0, nobs)

    # cache family functions
    if (is.function(family)) family <- family() # canonical link
    link_is_canonical <- (!is.null(family$linkcan)) && family$linkcan
    variance <- family$variance
    linkinv <- family$linkinv
    mu_eta <- family$mu.eta
    dev_resids <- family$dev.resids
    if (!is.function(variance) || !is.function(linkinv) ||
        !is.function(mu_eta) || !is.function(dev_resids))
      stop("'family' does not seem to be valid")

    # [ initialize ]
    if (!is.null(etastart)) {
      eta <- etastart
    } else {
      if (!is.null(start)) {
        if (length(start) != nvars)
          stop("inconsistent size of 'start'")
        eta <- drop(if (nvars == 1L) x * start else x %*% start)
      } else {
        # set mustart
        if (is.null(mustart)) {
          eval(family$initialize)
        } else {
          mukeep <- mustart
          eval(family$initialize)
          mustart <- mukeep
        }
        eta <- family$linkfun(mustart)
      }
    }
    eta <- eta + offset
    mu <- linkinv(eta)
    pdeviance <- Inf
    phi_new <- phi <- dispersion
    if (is.null(prior_coef))
      prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
    if (is.null(prior_disp))
      prior_disp <- list(df = 0, scale = 0)
    Omega <- prior_coef$precision
    beta0 <- mat_mult(Omega, drop(prior_coef$mean))
    nu <- prior_disp$df
    nu_tau2 <- nu * prior_disp$scale

    # [ iterate ]
    conv <- FALSE
    for (iter in 1L:control$maxit) {
      varmu <- variance(mu)
      W <- weights * varmu / phi
      residuals <- weights * (y - mu) / phi
      if (!link_is_canonical) { # adjust weights?
        thetaeta <- mu_eta(eta) / varmu
        W <- W * thetaeta ^ 2
        residuals <- residuals * thetaeta
      }
      z <- crossprod(x, W * (eta - offset) + residuals) + beta0
      C <- chol(mat_add(Omega, crossprod(x, x * W)))
      beta_new <- drop(chol_solve(C, z))
      eta <- drop(mat_mult(x, beta_new, nvars == 1L)) + offset
      mu <- linkinv(eta)
      deviance <- sum(dev_resids(y, mu, weights))
      if (adjust_dispersion)
        phi_new <- (nu_tau2 + deviance) / (nu + nobs) # E-step: E[1 / phi]

      beta_s <- mat_mult(Omega, beta_new)
      pdeviance_new <- phi_new * sum(beta_new * beta_s) + deviance
      if (control$trace)
        message("[", iter, "] Post. deviance = ", pdeviance_new)
      rel_error <- abs((pdeviance_new - pdeviance) / pdeviance)
      if (!is.infinite(pdeviance) && # converged?
          ((rel_error < control$epsilon) || (pdeviance_new < EPSILON))) {
        conv <- TRUE
        break
      }
      beta <- beta_new; phi <- phi_new; pdeviance <- pdeviance_new
    }
    residuals <- residuals / W # (y - mu) / mueta FIXME: check

    names(beta) <- xnames
    wtdmu <- if (intercept)
      (sum(weights * y) / phi + beta0[1]) / (sum(weights) / phi + Omega[1])
    else
      linkinv(offset + beta0[1])
    null_deviance <- sum(dev_resids(y, wtdmu, weights))
    n_ok <- nobs - sum(weights == 0)
    df_null <- n_ok - as.integer(intercept)
    df_residual <- n_ok - nvars # rank = nvars
    aic <- family$aic(y, n_ok, mu, weights, pdeviance) + 2 * nvars

    list(coefficients = beta, dispersion = phi, deviance = pdeviance,
         prior.coef = prior_coef, prior.disp = prior_disp,
         df.null = df_null, df.residual = df_residual,
         null.deviance = null_deviance, aic = aic, converged = conv,
         # extra
         qr = structure(list(qr = C, pivot = 1L:nvars), class = "qr"),
         residuals = residuals, fitted.values = mu,
         family = family, linear.predictors = eta, iter = iter,
         prior.weights = weights, weights = W, y = y, rank = nvars)
  }
}




# Extend `family` with d.V/d.mu and d(d.mu/d.eta)/d.eta for Hamiltonian MC
extend_family <- function (fam) {
  fam$variance.deriv <- switch(fam$family,
    gaussian = function (mu) rep(0, length(mu)),
    binomial = function (mu) 1 - 2 * mu,
    poisson = function (mu) rep(1, length(mu)),
    Gamma = function (mu) 2 * mu,
    inverse.gaussian = function (mu) 3 * mu ^ 2,
    stop(gettextf("%s family not recognized", sQuote(fam$family))))
  fam$mu.eta.deriv <- switch(fam$link,
    logit = function (eta) {
      me <- fam$mu.eta(eta)
      me - 2 * me ^ 2 * (1 + exp(eta))
    },
    probit = function (eta) -eta * dnorm(eta),
    cauchit = function (eta) -2 * eta / (pi * (1 + eta ^ 2) ^ 2),
    cloglog = function (eta) {
      ee <- exp(eta)
      exp(eta - ee) * (1 - ee)
    },
    identity = function (eta) rep.int(0, length(eta)),
    log = function (eta) pmax(exp(eta), EPSILON),
    sqrt = function (eta) rep.int(2, length(eta)),
    `1/mu^2` = function (eta) .75 / (eta ^ 2.5),
    inverse = function (eta) 2 / (eta ^ 3),
    stop(gettextf("%s link not recognized", sQuote(fam$link))))
  fam
}

mala_correction <- function (x, xw, C, delta) {
  p <- ncol(x)
  A <- xw %*% backsolve(C, diag(p))
  tg <- apply(x, 2, function (xj) sum((A ^ 2) * delta * xj))
  ug <- numeric(p)
  for (j in 1:p) {
    ug <- ug + crossprod(xw, xw * delta * x[, j]) %*%
      chol_solve(C, 1:p == j)
  }
  tg - 2 * ug
}

# TODO: graft into new `bsglm_sampler_hmc` and use `nsteps`
# `pm` is momentum
hmc_genleapfrog <- function (beta, pm, max_steps = 10, tol = 1e-6) {
  bc <- beta; pc <- pm
  for (hmc_step in 1:nsteps) {
    # [ First GLF step: p(t + eps / 2) ]
    # p(t + eps / 2) = p(t) - eps / 2 * dH/dbeta(beta(t), p(t + eps/2))
    eta <- drop(mat_mult(x, bc, nvars == 1L)) + offset
    mu <- linkinv(eta); mueta <- mu_eta(eta); var <- family$variance(mu)
    W <- weights * mueta ^ 2 / var / phi
    xw <- x * sqrt(W)
    C <- chol(mat_add(Omega, crossprod(xw)))
    # sg = dl/dbeta (score)
    sg <- crossprod(x, (y - mu) / mueta * W) - mat_mult(Omega, bc - beta0)
    # tg = [tr(G(beta)^{-1} * dG/dbeta_i]_i
    mueta_d <- mu_eta_deriv(eta); var_d <- variance_deriv(mu)
    delta <- 2 * mueta_d / mueta - mu * var_d / var
    A <- xw %*% backsolve(C, diag(ncol(x)))
    tg <- apply(x, 2, function (xj) sum((A ^ 2) * delta * xj))
    # solve via fixed point
    u <- pc - step_size / 2 * (-sg + .5 * tg)
    for (fi_step in 1:max_steps) {
      # ug = [pc' * G(beta)^{-1} * dG/dbeta_i * G(beta)^{-1} * pc]_i
      pg <- xw %*% chol_solve(C, pc)
      ug <- apply(x, 2, function (xj) sum((pg ^ 2) * delta * xj))
      pn <- u + step_size / 2 * ug # fixed point iteration
      if (sum(abs(pn - pc)) < tol) break
      pc <- pn
    }

    # [ Second GLF step: beta(t + eps) ]
    # beta(t + eps) = beta(t) + eps / 2 * (dH/dp(beta(t), p(t + eps/2)) +
    #                                      dH/dp(beta(t + eps), p(t + eps/2)))
    u <- bc + step_size / 2 * chol_solve(C, pc)
    for (fi_step in 1:max_steps) {
      eta <- drop(mat_mult(x, bc, nvars == 1L)) + offset
      mu <- linkinv(eta); mueta <- mu_eta(eta); var <- family$variance(mu)
      W <- weights * mueta ^ 2 / var / phi
      xw <- x * sqrt(W)
      C <- chol(mat_add(Omega, crossprod(xw)))
      bn <- u + step_size / 2 * chol_solve(C, pc) # fixed point iteration
      if (sum(abs(bn - bc)) < tol) break
      bc <- bn
    }
    # [ Third GLF step: p(t + eps) ]
    # p(t + eps) = p(t + eps / 2) - eps / 2 *
    #                   dH/dbeta(beta(t + eps), p(t + eps / 2))
    sg <- crossprod(x, (y - mu) / mueta * W) - mat_mult(Omega, bc - beta0)
    mueta_d <- mu_eta_deriv(eta); var_d <- variance_deriv(mu)
    delta <- 2 * mueta_d / mueta - mu * var_d / var
    A <- xw %*% backsolve(C, diag(ncol(x)))
    tg <- apply(x, 2, function (xj) sum((A ^ 2) * delta * xj))
    pg <- xw %*% chol_solve(C, pc)
    ug <- apply(x, 2, function (xj) sum((pg ^ 2) * delta * xj))
    pc <- pc - step_size / 2 * (-sg + .5 * tg - ug)
  }
  list(beta = bc, pm = pc)
}


#' Sample from a Bayesian generalized linear model
#'
#' Generate Hamiltonian MCMC samples for a Bayesian generalized
#' linear model with coefficients \eqn{\beta}{beta} and dispersion
#' \eqn{\phi}{\phi}. The coefficients have a normal prior with mean
#' \eqn{\beta_0}{beta0} and precision matrix \eqn{\Omega}{Omega},
#' \eqn{\beta \sim N(\beta_0, \Omega^{-1})}{beta ~ N(beta0, Omega^{-1})}, and
#' the dispersion has an inverse scaled chi-square prior with degrees of
#' freedom \eqn{\nu}{nu} and scale \eqn{\tau^2}{tau2},
#' \eqn{1/\phi \sim Ga(\nu/2, \nu\tau^2/2)}{1/phi ~ Ga(nu / 2, nu * tau2 / 2)}.
#'
#' This routine is a `glm` method and should be passed to a `glm` call, e.g.,
#' `glm(formula, data, family, ..., method = bsglm::sampler(...))`.
#'
#' @param dispersion value of dispersion if assumed deterministic or initial
#' value for fitting
#' @param adjust_dispersion is dispersion random?
#' @param prior_coef hyper-prior parameters for the coefficients, specified as
#' a list with entries `mean`, a vector, and `precision`, either a matrix or a
#' vector storing the diagonal entries
#' @param prior_disp hyper-prior parameters for the dispersion, specified as a
#' list with entries `df`, the number of degrees of freedom, and `scale`, the
#' variance prior scale
#' @param step_size Step size in leap frog solver
#' @param nsteps Number of steps in leap frog solver
#' @param simplified Simplified version? With `nsteps = 1` this is MALA
#' @return a `glm` object.
#' @export
# NOTE: check accepting rate with `mean(abs(diff(g$samples$coef[,1])) > 1e-10)`
sampler <- function (dispersion = 1, adjust_dispersion = FALSE,
                     prior_coef = NULL, prior_disp = NULL,
                     step_size = 1, nsteps = 1, simplified = TRUE) {
  function (x, y, weights = NULL, start = NULL,
            etastart = NULL, mustart = NULL, offset = NULL,
            family = gaussian(), control = list(),
            intercept = TRUE, singular.ok = TRUE) {
    control <- do.call("control", control)
    x <- as.matrix(x)
    nobs <- nrow(x)
    nvars <- ncol(x)
    xnames <- colnames(x)
    if (is.null(xnames)) xnames <- paste0("x", 1L:nvars)
    if (nvars == 0L) stop("empty models are not allowed")
    if (is.null(weights)) weights <- rep.int(1, nobs)
    if (is.null(offset)) offset <- rep.int(0, nobs)
    if (is.null(start)) start <- prior_coef$mean
    ss2 <- step_size ^ 2 / 2

    # cache family functions
    if (is.function(family)) family <- family() # canonical link
    variance <- family$variance
    linkinv <- family$linkinv
    mu_eta <- family$mu.eta
    dev_resids <- family$dev.resids
    if (!is.function(variance) || !is.function(linkinv) ||
        !is.function(mu_eta) || !is.function(dev_resids))
      stop("'family' does not seem to be valid")
    if (!simplified) {
      family <- extend_family(family)
      mu_eta_deriv <- family$mu.eta.deriv
      variance_deriv <- family$variance.deriv
    }

    # [ initialize ]
    if (is.null(prior_coef))
      prior_coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
    if (is.null(prior_disp))
      prior_disp <- list(df = 0, scale = 0)
    Omega <- prior_coef$precision
    beta0 <- drop(prior_coef$mean)
    nu <- prior_disp$df
    nu_tau2 <- nu * prior_disp$scale

    phi <- dispersion
    beta <- drop(if (is.null(start)) prior_coef$mean else start)
    eta <- drop(mat_mult(x, beta, nvars == 1L)) + offset
    mu <- linkinv(eta); mueta <- mu_eta(eta); var <- variance(mu)
    deviance <- sum(dev_resids(y, mu, weights))
    W <- weights * mueta ^ 2 / var / phi
    xw <- x * sqrt(W)
    C <- chol(mat_add(Omega, crossprod(xw)))
    beta_mu <- crossprod(x, (y - mu) / mueta * W) -
      mat_mult(Omega, beta - beta0) # score
    if (!simplified) {
      mueta_d <- mu_eta_deriv(eta); var_d <- variance_deriv(mu)
      delta <- 2 * mueta_d / mueta - mu * var_d / var
      beta_mu <- beta_mu + mala_correction(x, xw, C, delta)
    }
    beta_mu <- drop(beta + ss2 * chol_solve(C, beta_mu))

    sims <- mcmc_init_array(control$nsamples, control$nchains,
                            c(xnames, "dispersion", "__lp"))
    # [ iterate ]
    for (chain in 1:control$nchains) {
      for (iter in 1:control$nsamples) {
        beta_c <- drop(beta_mu + step_size * backsolve(C, rnorm(nvars)))
        eta_c <- drop(mat_mult(x, beta_c, nvars == 1L)) + offset
        mu_c <- linkinv(eta_c); mueta_c <- mu_eta(eta_c); var_c <- variance(mu_c)
        deviance_c <- sum(dev_resids(y, mu_c, weights))
        W_c <- weights * mueta_c ^ 2 / var_c / phi
        xw_c <- x * sqrt(W_c)
        C_c <- chol(mat_add(Omega, crossprod(xw_c)))
        beta_mu_c <- crossprod(x, (y - mu_c) / mueta_c * W_c) -
          mat_mult(Omega, beta_c - beta0) # score
        if (!simplified) {
          mueta_d_c <- mu_eta_deriv(eta_c); var_d_c <- variance_deriv(mu_c)
          delta_c <- 2 * mueta_d_c / mueta_c - mu_c * var_d_c / var_c
          beta_mu_c <- beta_mu_c + mala_correction(x, xw_c, C_c, delta_c)
        }
        beta_mu_c <- drop(beta_c + ss2 * chol_solve(C_c, beta_mu_c))

        if (is.vector(Omega)) {
          lc <- -.5 * sum(Omega * (beta_c - beta0) ^ 2)
          lb <- -.5 * sum(Omega * (beta - beta0) ^ 2)
        }
        else {
          u <- beta_c - beta0
          lc <- -.5 * crossprod(u, crossprod(Omega, u))
          u <- beta - beta0
          lb <- -.5 * crossprod(u, crossprod(Omega, u))
        }
        lqc <- sum(log(diag(C))) -
          .5 * crossprod(C %*% (beta_c - beta_mu) / step_size)
        lc <- lc - .5 * deviance_c / phi
        lcq <- lc - lqc
        lqb <- sum(log(diag(C_c))) -
          .5 * crossprod(C_c %*% (beta - beta_mu_c) / step_size)
        lb <- lb - .5 * deviance / phi
        lbq <- lb - lqb
        lp <- lb
        if (lcq >= lbq || log(runif(1)) <= lcq - lbq) { # accept?
          beta <- beta_c; deviance <- deviance_c
          mu <- mu_c; mueta <- mueta_c; var <- var_c
          if (!simplified) delta <- delta_c
          # dispersion-dependent:
          beta_mu <- beta_mu_c; C <- C_c
          lp <- lc
        }

        if (adjust_dispersion) {
          # sample
          phi <- rinvchisq(1, nu + nobs, nu_tau2 + deviance)
          # re-cache beta_mu and C
          W <- weights * mueta ^ 2 / var / phi
          xw <- x * sqrt(W)
          C <- chol(mat_add(Omega, crossprod(xw)))
          beta_mu <- crossprod(x, (y - mu) / mueta * W) -
            mat_mult(Omega, beta - beta0) # score
          if (!simplified)
            beta_mu <- beta_mu + mala_correction(x, xw, C, delta)
          beta_mu <- drop(beta + ss2 * chol_solve(C, beta_mu))
        }

        sims[iter, chain, 1:nvars] <- beta
        sims[iter, chain, nvars + 1] <- phi
        sims[iter, chain, nvars + 2] <- lp
      }
    }

    # summarize
    phi <- mean(sims[, , nvars + 1])
    beta <- apply(sims, 3, mean)[1L:nvars]; names(beta) <- xnames
    eta <- drop(if (nvars == 1L) x * beta else x %*% beta) + offset
    mu <- linkinv(eta); mueta <- mu_eta(eta); var <- variance(mu)
    residuals <- (y - mu) / mueta
    deviance <- sum(dev_resids(y, mu, weights))
    W <- weights * mueta ^ 2 / var / phi
    C <- chol(mat_add(Omega, crossprod(x, x * W)))
    wtdmu <- if (intercept)
      (sum(weights * y) / phi + beta0[1]) / (sum(weights) / phi + Omega[1])
    else
      linkinv(offset + beta0[1])
    null_deviance <- sum(dev_resids(y, wtdmu, weights))
    n_ok <- nobs - sum(weights == 0)
    df_null <- n_ok - as.integer(intercept)
    df_residual <- n_ok - nvars # rank = nvars
    aic <- family$aic(y, n_ok, mu, weights, deviance) + 2 * nvars

    list(coefficients = beta, dispersion = phi, deviance = deviance,
         prior.coef = prior_coef, prior.disp = prior_disp,
         samples = sims,
         df.null = df_null, df.residual = df_residual,
         null.deviance = null_deviance, aic = aic,
         converged = TRUE, # since we're sampling
         # extra
         qr = structure(list(qr = C, pivot = 1L:nvars), class = "qr"),
         residuals = residuals, fitted.values = mu,
         family = family, linear.predictors = eta, iter = iter,
         prior.weights = weights, weights = W, y = y, rank = nvars)
  }
}

# TODO:
# - semi-conjugate option
# - move `mala_correction` and parts of leap frog to C

