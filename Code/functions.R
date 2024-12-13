# Iteratively reweighted least squares for canonical links
bglm.fit <- function (formula, family, weights = 1,
                      dispersion = 1, adjust.dispersion = F,
                      prior.coef = list(mean = 0, precision = 0),
                      prior.disp = list(nu = 0, nu.tau2 = 0),
                      offset = 0, maxit = 25, tol = 1e-8, int = 1) {
  y <- model.response(model.frame(formula))
  X <- model.matrix(formula) # design matrix
  
  if (int == 0) X <- X[,-1]
  
  fam <- family
  if (class(family) == 'function') fam <- family() # canonical link
  nobs <- length(y)
  
  etastart = 0
  eval(fam$initialize) # set mustart
  offset = as.numeric(offset) #change offset to numeric
  eta <- fam$linkfun(mustart) + offset
  
  mu <- fam$linkinv(eta)
  logpost <- Inf
  phi.new <- phi <- dispersion
  S.inv0 <- prior.coef$precision
  if (is.matrix(S.inv0))
    beta0 <- S.inv0 %*% prior.coef$mean
  else
    beta0 <- S.inv0 * prior.coef$mean
  
  for (iter in 1:maxit) {
    mu.eta <- fam$mu.eta(eta) # 1 / g'(mu) [ = V(mu) for canonical link]
    W <- weights * mu.eta ^ 2 / fam$variance(mu) / phi
    z <- crossprod(X, ((eta - offset) + (y - mu) / mu.eta) * W) + beta0
    V.inv <- crossprod(X, X * W)
    if (is.vector(S.inv0))
      diag(V.inv) <- diag(V.inv) + S.inv0
    else
      V.inv <- V.inv + S.inv0
    C <- chol(V.inv)
    beta.new <- backsolve(C, backsolve(C, z, transpose=TRUE))
    eta <- drop(X %*% beta.new) + offset
    mu <- fam$linkinv(eta)
    deviance <- sum(fam$dev.resids(y, mu, weights))
    if (adjust.dispersion)
      phi.new <- (prior.disp$nu.tau2 + deviance) / (prior.disp$nu + nobs + 2)
    
    if (is.matrix(S.inv0))
      logpost.new <- crossprod(beta.new, S.inv0 %*% beta.new)
    else
      logpost.new <- crossprod(beta.new, S.inv0 * beta.new)
    logpost.new <- phi.new * logpost.new + deviance
    rel.error <- abs(logpost.new - logpost) / abs(logpost)
    if (!is.infinite(logpost) && (rel.error < tol)) # converged?
      break
    beta <- beta.new; phi <- phi.new; logpost <- logpost.new
  }
  
  list(coef = drop(beta), dispersion = phi, var = chol2inv(C), C = C)
}

# Bayesian Linear Models
# y | beta, sigma2 ~ N(X * beta, sigma2 * I_n)
# beta ~ N(beta0, S.inv0^(-1))
# sigma2 ~ Inv-scaled-Chisq(nu, tau2)

# [ Auxiliary ]
rinvchisq <- function (n, nu, nu.tau2) 1 / rgamma(n, nu / 2, nu.tau2 / 2)

mcmc.init.array <- function (ns, nc, params)
  array(dim=c(ns, nc, length(params)),
        dimnames=list(iterations=NULL,
                      chains=paste0("chain:", 1L:nc),
                      parameters=params))

# [ Simple interface ]

bslm.fit <- function (y, x, prior.coef = NULL, prior.disp = NULL,
                      maxit = 25, epsilon = 1e-8) {
  nvars <- ncol(x); nobs <- nrow(x)
  dn <- colnames(x); if (is.null(dn)) dn <- paste0("x", 1L:nvars)
  if (is.null(prior.coef))
    prior.coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
  if (is.null(prior.disp))
    prior.disp <- list(df = 0, scale = 0)
  S.inv0 <- prior.coef$precision
  beta0 <- prior.coef$mean
  beta0 <- if (is.vector(S.inv0)) S.inv0 * beta0 else S.inv0 %*% beta0
  nu <- prior.disp$df
  nu.tau2 <- nu * prior.disp$scale
  
  rss <- sum((y - mean(y)) ^ 2)
  sigma2 <- (nu.tau2 + rss) / (nu + nobs)
  for (iter in 1:maxit) {
    z <- crossprod(x, y) / sigma2 + beta0
    V.inv <- crossprod(x) / sigma2
    if (is.vector(S.inv0))
      diag(V.inv) <- diag(V.inv) + S.inv0
    else
      V.inv <- V.inv + S.inv0
    C <- chol(V.inv)
    coef <- backsolve(C, z, transpose=TRUE)
    coef <- drop(backsolve(C, coef))
    rss.new <- sum((y - drop(x %*% coef)) ^ 2)
    sigma2 <- (nu.tau2 + rss.new) / (nu + nobs)
    
    rel.error <- abs((rss.new - rss) / rss)
    if (!is.infinite(rss.new) && (rel.error < epsilon))
      break
    rss <- rss.new
  }
  names(coef) <- dn
  
  list(coef=coef, sigma2=sigma2, C=C)
}


bslm.sample <- function (y, x, prior.coef = NULL, prior.disp = NULL,
                         nsamples = 1000, nchains = 1) {
  nvars <- ncol(x); nobs <- nrow(x)
  dn <- colnames(x); if (is.null(dn)) dn <- paste0("x", 1L:nvars)
  if (is.null(prior.coef))
    prior.coef <- list(mean = rep(0, nvars), precision = rep(0, nvars))
  if (is.null(prior.disp))
    prior.disp <- list(df = 0, scale = 0)
  S.inv0 <- prior.coef$precision
  beta0 <- prior.coef$mean
  beta0 <- if (is.vector(S.inv0)) S.inv0 * beta0 else S.inv0 %*% beta0
  nu <- prior.disp$df
  nu.tau2 <- nu * prior.disp$scale
  
  rss <- sum((y - mean(y)) ^ 2)
  sigma2 <- (nu.tau2 + rss) / (nu + nobs)
  sims <- mcmc.init.array(nsamples, nchains, c(dn, "sigma2"))
  for (chain in 1:nchains) {
    for (iter in 1:nsamples) {
      z <- crossprod(x, y) / sigma2 + beta0
      V.inv <- crossprod(x) / sigma2
      if (is.vector(S.inv0))
        diag(V.inv) <- diag(V.inv) + S.inv0
      else
        V.inv <- V.inv + S.inv0
      C <- chol(V.inv)
      coef <- backsolve(C, z, transpose=TRUE)
      coef <- drop(backsolve(C, coef + rnorm(nvars)))
      rss <- sum((y - drop(x %*% coef)) ^ 2)
      sigma2 <- rinvchisq(1, nu + nobs, nu.tau2 + rss)
      
      sims[iter, chain, 1:nvars] <- coef
      sims[iter, chain, nvars + 1] <- sigma2
    }
  }
  
  sims
}

elbow_finder <- function(x_values, y_values) {
  # https://stackoverflow.com/a/42808962
  
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)
  
  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}

jaccard_dist <- function(x, y) length(intersect(x, y)) / length(union(x, y))