# job script for submission to HPC
rm(list = ls())
gc()
source("functions.R")
library(foreach)
library(doSNOW)
library(parallel)
library(dplyr)

# define scenarios
scenarios = expand.grid(
  n_individuals = c(25, 50, 250, 500, 1000),
  n_clusters = c(10, 25, 100, 200),
  frailty_sigma = c(0.25, 0.50, 0.75),
  treatment_effect = c(-0.50, 0.00, 0.50),
  ngl = c(35, 75, 105),
  lambda = 0.5,
  p = 1) %>%
  filter((n_individuals %in% c(25, 50) & n_clusters %in% c(100, 200)) | (n_individuals %in% c(250, 500, 1000) & n_clusters %in% c(10, 25)))

# get parameter from array id
TID = commandArgs(trailingOnly = T)

# setup parallel cluster
cl <- makeCluster(10, type = "SOCK")
registerDoSNOW(cl)

# generate seeds for the simulation
set.seed(TID)

# run 1000 simulations for the current scenario
s = foreach(seeds = round(runif(1000, 0, 1e6)), .combine = rbind) %dopar% {
  with(scenarios[TID,],
       sim_normal_gq(seed = seeds,
                     n_individuals = n_individuals,
                     n_clusters = n_clusters,
                     frailty_sigma = frailty_sigma,
                     treatment_effect = treatment_effect,
                     ngh = ngh,
                     lambda = lambda,
                     p = p))
}

# stop the parallel cluster
stopCluster(cl)

# save results
saveRDS(s, paste0("Results/s_normal_gq_", TID, ".RDS"))
