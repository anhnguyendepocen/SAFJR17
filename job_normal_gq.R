# job script for submission to HPC
rm(list = ls())
gc()
source("functions.R")
if (!requireNamespace("pacman")) install.packages("pacman")
pacman::p_load("foreach", "doSNOW", "parallel", "dplyr")

# define scenarios
scenarios = expand.grid(
  n_individuals = c(25, 50, 100, 250, 500, 1000),
  n_clusters = c(25, 50, 100, 150),
  frailty_sigma = c(0.25, 0.50, 1.00),
  treatment_effect = c(-0.50, 0.00, 0.50),
  lambda = 3,
  p = 1.5) %>%
  filter(!(n_individuals %in% c(500, 1000) & n_clusters == 100) & !(n_individuals %in% c(250, 500, 1000) & n_clusters == 150) & !(n_individuals == 1000 & n_clusters == 50))

# get parameter from array id
TID = commandArgs(trailingOnly = T)

# setup parallel cluster
cl <- makeSOCKcluster(24)
registerDoSNOW(cl)

# generate seeds for each simulation
set.seed(TID)
sds <- round(runif(1000, 0, 1e7))

# run 1000 simulations for the current scenario
s = foreach(siid = sds, .combine = bind_rows) %dopar% {
  with(scenarios[TID,],
       sim_normal_gq(seed = siid,
                     n_individuals = n_individuals,
                     n_clusters = n_clusters,
                     frailty_sigma = frailty_sigma,
                     treatment_effect = treatment_effect,
                     lambda = lambda,
                     p = p))
}

# stop the parallel cluster
stopCluster(cl)

# save results
saveRDS(s, paste0("Results/s_normal_gq_", TID, ".RDS"))
