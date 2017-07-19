### functions for generating data

make_data <- function(n_individuals, n_clusters, fv_dist, fv, treatment_effect, lambda, p, maxt = 5, scenario) {
  # cluster id
  grpid = rep(1:n_clusters, each = n_individuals)
  # get RNG seed
  seed = .Random.seed
  # treatment effect
  trt = rbinom(n_individuals * n_clusters, size = 1, prob = 0.5)
  # frailty
  if (fv_dist == "Gamma") {
    fr = rgamma(n_clusters, shape = 1 / fv, scale = fv)
  } else {
    fr = exp(rnorm(n_clusters, mean = 0, sd = sqrt(fv)))
  }
  frvec = rep(fr, each = n_individuals)
  # draw from uniform
  uu = runif(n_individuals * n_clusters)
  # survival times
  s = (-log(uu) / (lambda * exp(trt * treatment_effect) * frvec)) ^ (1 / p)
  # applying administrative censoring at `maxt` years
  t = pmin(s, maxt)
  # event indicator variable
  d = as.numeric(s <= maxt)
  out = data.frame(grpid, trt, frvec, t, d, n_individuals, n_clusters, fv_dist, fv, treatment_effect, lambda, p, scenario)
  attr(out, "seed") = seed
  return(out)
}

### define DGMs
if (!requireNamespace("pacman")) install.packages("pacman")
pacman::p_load("dplyr", "tidyr")

### data for an_vs_gq comparison
dgms <- crossing(
  n_individuals = c(2, 30, 100, 500),
  n_clusters = c(15, 50, 1000),
  treatment_effect = c(-0.50, 0.00, 0.50),
  fv = c(0.25, 0.75, 1.25),
  fv_dist = "Gamma",
  lambda = 0.5,
  p = 0.6) %>%
  filter((n_individuals * n_clusters) <= 10000 & (n_individuals * n_clusters) > 300) %>%
  arrange(p, lambda, fv_dist, fv, treatment_effect, n_individuals, n_clusters)

### generate B datasets for each scenario
B = 1000
set.seed(294635876)
lapply(1:nrow(dgms),
       function(d) {
         cat(paste0(d, "\n"))
         out = replicate(B,
                         make_data(n_individuals = dgms[d,]$n_individuals,
                                   n_clusters = dgms[d,]$n_clusters,
                                   fv = dgms[d,]$fv,
                                   fv_dist = dgms[d,]$fv_dist,
                                   treatment_effect = dgms[d,]$treatment_effect,
                                   lambda = dgms[d,]$lambda,
                                   p = dgms[d,]$p,
                                   scenario = d),
                         simplify = FALSE)
         saveRDS(object = out, file = paste0("Data/simdata_an_vs_gq_", d, ".RDS"))
       }) %>%
  invisible()

### data for normal_gq comparison
dgms <- crossing(
  n_individuals = c(2, 30, 100, 500),
  n_clusters = c(15, 50, 1000),
  treatment_effect = c(-0.50, 0.00, 0.50),
  fv = c(0.25, 0.75, 1.25),
  fv_dist = "Log-Normal",
  lambda = 0.5,
  p = 0.6) %>%
  filter((n_individuals * n_clusters) <= 10000 & (n_individuals * n_clusters) > 300) %>%
  arrange(p, lambda, fv_dist, fv, treatment_effect, n_individuals, n_clusters)

### generate B datasets for each scenario
B = 1000
set.seed(1853497892)
lapply(1:nrow(dgms),
       function(d) {
         cat(paste0(d, "\n"))
         out = replicate(B,
                         make_data(n_individuals = dgms[d,]$n_individuals,
                                   n_clusters = dgms[d,]$n_clusters,
                                   fv = dgms[d,]$fv,
                                   fv_dist = dgms[d,]$fv_dist,
                                   treatment_effect = dgms[d,]$treatment_effect,
                                   lambda = dgms[d,]$lambda,
                                   p = dgms[d,]$p,
                                   scenario = d),
                         simplify = FALSE)
         saveRDS(object = out, file = paste0("Data/simdata_normal_gq_", d, ".RDS"))
       }) %>%
  invisible()
