### function for generating data

gen_data <- function(seed, n_individuals, n_clusters, frailty_theta, treatment_effect, lambda, p) {
  # set seed
  set.seed(seed)
  # cluster id
  grpid = rep(1:n_clusters, each = n_individuals)
  # treatment effect
  trt = rbinom(n_individuals * n_clusters, size = 1, prob = 0.5)
  # Gamma frailty
  fr = rgamma(n_clusters, shape = 1 / frailty_theta, scale = frailty_theta)
  frvec = rep(fr, each = n_individuals)
  # survival times
  s = (-log(runif(n_individuals * n_clusters)) / (lambda * exp(trt * treatment_effect) * frvec)) ^ (1 / p)
  # censoring time
  c = (-log(runif(n_individuals * n_clusters)) / lambda) ^ (1 / p)
  # survival time
  t = pmin(s, c)
  # event indicator variable
  d = as.numeric(s <= c)
  return(data.frame(grpid, trt, frvec, t, d))
}
# gen_data(seed = 1234, n_individuals = 100, n_clusters = 5, frailty_theta = 1, treatment_effect = 0.5, lambda = 3, p = 1.5) %>% View()

gen_data_normal <- function(seed, n_individuals, n_clusters, frailty_sigma, treatment_effect, lambda, p) {
  # set seed
  set.seed(seed)
  # cluster id
  grpid = rep(1:n_clusters, each = n_individuals)
  # treatment status
  trt = rbinom(n_individuals * n_clusters, size = 1, prob = 0.5)
  # random effect + treatment effect
  reff = rnorm(n_clusters, mean = 0, sd = frailty_sigma)
  reffvec = rep(reff, each = n_individuals)
  # survival times
  s = (-log(runif(n_individuals * n_clusters)) / (lambda * exp(trt * (treatment_effect + reffvec)))) ^ (1 / p)
  # censoring time
  c = (-log(runif(n_individuals * n_clusters)) / lambda) ^ (1 / p)
  # survival time
  t = pmin(s, c)
  # event indicator variable
  d = as.numeric(s <= c)
  return(data.frame(grpid, trt, reffvec, t, d))
}
# gen_data_normal(seed = 1234, n_individuals = 100, n_clusters = 5, frailty_sigma = 0.50, treatment_effect = 0.5, lambda = 3, p = 1.5) %>% View()

### functions for model estimation

# function for estimating a model with a (potentially shared) Gamma frailty term
sim_an_vs_gq <- function(seed, n_individuals, n_clusters, frailty_theta, treatment_effect, lambda, p, ngl = 35) {
  # packages
  if (!requireNamespace("pacman")) install.packages("pacman")
  pacman::p_load("pracma")

  if (ngl < 2) stop("number of knots 'ngl' needs to be >= 2")
  # generate data
  # seed = 123; n_individuals = 100; n_clusters = 5; frailty_theta = 0.25; treatment_effect = 0.50; lambda = 0.50; p = 1; ngl = 35
  df = gen_data(seed = seed,
                n_individuals = n_individuals,
                n_clusters = n_clusters,
                frailty_theta = frailty_theta,
                treatment_effect = treatment_effect,
                lambda = lambda,
                p = p)

  # starting parameters
  sr = survival::survreg(survival::Surv(t, d) ~ trt, dist = "weibull", data = df)
  # we need to change survreg's parametrisation
  # 1/sr$scale # p, shape
  # exp(sr$coef[1]) ^ (-1/sr$scale) # lambda, scale
  start = c(log(1 / sr$scale), # p
            log(exp(sr$coef[1]) ^ (-1 / sr$scale)), # lambda
            log(1), # theta
            -sr$coef[2] * (1 / sr$scale)) # beta1, coefficient of trt, AFT --> PH
  names(start) <- c("p", "lambda", "theta", "trt")

  # analytical likelihood
  if (n_clusters == 1) {
    mll = function(pars){
      p = exp(pars[1])
      lambda = exp(pars[2])
      theta = exp(pars[3])
      beta1 = pars[4]
      log_hi = log(p) + log(lambda) + (p - 1) * log(df$t) + df$trt * beta1
      log_Si = -lambda * df$t ^ p * exp(df$trt * beta1)
      ll = sum(df$d * (log_hi) - (theta ^ (-1) + df$d) * log(1 - theta * log_Si)) + sum(log(df$t))
      return(-ll)
    }
  } else {
    mll = function(pars) {
      p = exp(pars[1])
      lambda = exp(pars[2])
      theta = exp(pars[3])
      beta1 = pars[4]
      lli = vapply(1:n_clusters,
                   FUN = function(i) {
                     log_hi = log(p) + log(lambda) + (p - 1) * log(df$t[df$grpid == i]) + df$trt[df$grpid == i] * beta1
                     log_Si = -lambda * df$t[df$grpid == i] ^ p * exp(df$trt[df$grpid == i] * beta1)
                     Di = sum(df$d[df$grpid == i])
                     lli = sum(df$d[df$grpid == i] * (log_hi)) - (theta ^ (-1) + Di) * log(1 - theta * sum(log_Si)) + Di * log(theta) + lgamma(theta ^ (-1) + Di) - lgamma(1 / theta) + sum(log(df$t[df$grpid == i]))
                     return(lli)
                   },
                   FUN.VALUE = numeric(1))
      ll = sum(lli)
      return(-ll)
    }
  }

  # Gauss-Laguerre knots and weights
  gl_rule = gaussLaguerre(ngl)

  # quadrature likelihood
  mll_quad = function(pars) {
    p = exp(pars[1])
    lambda = exp(pars[2])
    theta = exp(pars[3])
    beta1 = pars[4]

    lli = vapply(1:n_clusters,
                 FUN = function(i) {
                   log_hi = log(p) + log(lambda) + (p - 1) * log(df$t[df$grpid == i]) + (df$trt[df$grpid == i] * beta1)
                   log_Si = -lambda * df$t[df$grpid == i] ^ p * exp(df$trt[df$grpid == i] * beta1)
                   Di = sum(df$d[df$grpid == i])
                   intgrdq = function(alpha) exp(alpha + Di * log(alpha) + alpha * sum(log_Si) + (1 / theta - 1) * log(alpha) - alpha / theta)
                   vintgrdq = Vectorize(intgrdq)
                   int = sum(gl_rule$w * vintgrdq(gl_rule$x))
                   ll = sum(df$d[df$grpid == i] * (log_hi)) - lgamma(1 / theta) - (1 / theta) * log(theta) + log(int) + sum(log(df$t[df$grpid == i]))
                   return(ll)},
                 FUN.VALUE = numeric(1))
    ll = sum(lli)
    return(-ll)
  }

  # mll_quad(start)

  # error object to return in case of optim errors
  error_obj = list(
    par = c(p = NA, lambda = NA, theta = NA, trt = NA),
    value = NA,
    counts = c(`function` = NA, `gradient` = NA),
    convergence = -99,
    message = NULL,
    hessian = matrix(NA, nrow = 4, ncol = 4, dimnames = list(c("p", "lambda", "theta", "trt"), c("p", "lambda", "theta", "trt")))
  )

  # optimisation routines
  o_an = tryCatch(optim(par = start, fn = mll, method = "BFGS", hessian = TRUE),
                  error = function(cond) {
                    return(error_obj)
                  })
  if (o_an$convergence != 0) {
    o_an$se = c(p = NA, lambda = NA, theta = NA, trt = NA)
  } else {
    o_an$se = sqrt(diag(solve(o_an$hessian)))
  }

  o_gq = tryCatch(optim(par = start, fn = mll_quad, method = "BFGS", hessian = TRUE),
                  error = function(cond) {
                    return(error_obj)
                  })
  if (o_gq$convergence != 0) {
    o_gq$se = c(p = NA, lambda = NA, theta = NA, trt = NA)
  } else {
    o_gq$se = sqrt(diag(solve(o_gq$hessian)))
  }

  output = data.frame(seed = seed, n_individuals = n_individuals, n_clusters = n_clusters, frailty_theta = frailty_theta, treatment_effect = treatment_effect, lambda = lambda, p = p, ngl = ngl, AF_p = o_an$par[1], AF_p_se = o_an$se[1], AF_lambda = o_an$par[2], AF_lambda_se = o_an$se[2], AF_theta = o_an$par[3], AF_theta_se = o_an$se[3], AF_trt = o_an$par[4], AF_trt_se = o_an$se[4], AF_value = -o_an$value, AF_convergence = o_an$convergence, GQ_p = o_gq$par[1], GQ_p_se = o_gq$se[1], GQ_lambda = o_gq$par[2], GQ_lambda_se = o_gq$se[2], GQ_theta = o_gq$par[3], GQ_theta_se = o_gq$se[3], GQ_trt = o_gq$par[4], GQ_trt_se = o_gq$se[4], GQ_value = -o_gq$value, GQ_convergence = o_gq$convergence)
  rownames(output) <- NULL

  # return results
  return(output)
}

# function for estimating a model with a (potentially shared) Normal frailty term
sim_normal_gq <- function(seed, n_individuals, n_clusters, frailty_sigma, treatment_effect, lambda, p, ngh = 35) {
  # packages
  if (!requireNamespace("pacman")) install.packages("pacman")
  pacman::p_load("fastGHQuad")

  if (ngh < 2) stop("number of knots 'ngh' needs to be >= 2")
  # generate data
  # seed = 352486; n_individuals = 1000; n_clusters = 15; frailty_sigma = 1; treatment_effect = 0.5; lambda = 3; p = 1.5; ngh = 35
  df = gen_data_normal(seed = seed,
                       n_individuals = n_individuals,
                       n_clusters = n_clusters,
                       frailty_sigma = frailty_sigma,
                       treatment_effect = treatment_effect,
                       lambda = lambda,
                       p = p)

  # starting parameters
  sr = survival::survreg(survival::Surv(t, d) ~ trt, dist = "weibull", data = df)
  # we need to change survreg's parametrisation
  # 1/sr$scale # p, shape
  # exp(sr$coef[1]) ^ (-1/sr$scale) # lambda, scale
  start = c(log(1 / sr$scale), # p
            log(exp(sr$coef[1]) ^ (-1 / sr$scale)), # lambda
            log(1), # sigma
            -sr$coef[2] * (1 / sr$scale)) # beta1, coefficient of trt, AFT --> PH
  names(start) <- c("p", "lambda", "sigma", "trt")

  # Gauss-Hermite knots and weights, including adjustment:
  gh_rule = gaussHermiteData(ngh)

  # quadrature likelihood
  mll = function(pars) {
    p = exp(pars[1])
    lambda = exp(pars[2])
    sigma = exp(pars[3])
    beta1 = pars[4]

    lli = vapply(1:n_clusters,
                 FUN = function(i) {
                   intgrd = function(bi) {
                     log_hi = log(p) + log(lambda) + (p - 1) * log(df$t[df$grpid == i]) + (df$trt[df$grpid == i] * (beta1 + bi))
                     log_Si = -lambda * df$t[df$grpid == i] ^ p * exp(df$trt[df$grpid == i] * (beta1 + bi))
                     # ret = exp(sum(df$d[df$grpid == i] * (log_hi + log(df$t[df$grpid == i]))) + sum(log_Si) - (bi ^ 2) / (2 * sigma ^ 2))
                     ret = exp(sum(df$d[df$grpid == i] * (log_hi)) + sum(log_Si))
                     return(ret)}
                   vintgrd = Vectorize(intgrd)
                   int = sum(gh_rule$w / sqrt(pi) * vintgrd(gh_rule$x * sqrt(2) * sigma))
                   # ll = -(1/2) * log(2 * (sigma ^ 2) * pi) + log(int)
                   ll = log(int) + sum(log(df$t[df$grpid == i]))
                   return(ll)},
                 FUN.VALUE = numeric(1))
    ll = sum(lli)
    return(-ll)
  }

  # mll(start)

  # error object to return in case of optim errors
  error_obj = list(
    par = c(p = NA, lambda = NA, sigma = NA, trt = NA),
    value = NA,
    counts = c(`function` = NA, `gradient` = NA),
    convergence = -99,
    message = NULL,
    hessian = matrix(NA, nrow = 4, ncol = 4, dimnames = list(c("p", "lambda", "sigma", "trt"), c("p", "lambda", "sigma", "trt")))
  )

  # optimisation routines
  o_gq = tryCatch(
    optim(par = start, fn = mll, method = "BFGS", hessian = TRUE),
    error = function(cond) return(error_obj))

  if (o_gq$convergence != 0) {
    o_gq$se = c(p = NA, lambda = NA, sigma = NA, trt = NA)
  } else {
    o_gq$se = sqrt(diag(solve(o_gq$hessian)))
  }

  output = data.frame(seed = seed, n_individuals = n_individuals, n_clusters = n_clusters, frailty_sigma = frailty_sigma, treatment_effect = treatment_effect, lambda = lambda, p = p, ngh = ngh, GQ_p = o_gq$par[1], GQ_p_se = o_gq$se[1], GQ_lambda = o_gq$par[2], GQ_lambda_se = o_gq$se[2], GQ_sigma = o_gq$par[3], GQ_sigma_se = o_gq$se[3], GQ_trt = o_gq$par[4], GQ_trt_se = o_gq$se[4], GQ_value = -o_gq$value, GQ_convergence = o_gq$convergence)
  rownames(output) <- NULL

  # return results
  return(output)
}
