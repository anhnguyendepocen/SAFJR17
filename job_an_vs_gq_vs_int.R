# job script for submission to HPC
rm(list = ls())
gc()
source("functions.R")
library(foreach)
library(doSNOW)
library(parallel)

# define scenarios
scenarios = expand.grid(
  n_individuals = c(25, 50, 100, 250, 500, 1000),
  n_clusters = c(15, 30, 100, 200),
  frailty_theta = c(0.25, 0.50, 1.00),
  treatment_effect = c(-0.50, 0.00, 0.50),
  lambda = 0.5,
  p = 1)

scenarios = scenarios[!(scenarios$n_individuals %in% c(500, 1000) & scenarios$n_clusters == 100) & !(scenarios$n_individuals %in% c(250, 500, 1000) & scenarios$n_clusters == 200),]

# get parameter from array id
TID = commandArgs(trailingOnly = T)

# setup parallel cluster
cl <- makeCluster(20, type = "SOCK")
registerDoSNOW(cl)

# generate seeds for the simulation
set.seed(TID)

# run 1000 simulations for the current scenario
s = foreach(seeds = round(runif(1000, 0, 1e9)), .combine = rbind) %dopar% {
  with(scenarios[TID,],
       sim_an_vs_gq_vs_int(seed = seeds,
                           n_individuals = n_individuals,
                           n_clusters = n_clusters,
                           frailty_theta = frailty_theta,
                           treatment_effect = treatment_effect,
                           lambda = lambda,
                           p = p))
}

# stop the parallel cluster
stopCluster(cl)

# save results
saveRDS(s, paste0("Results/s_an_vs_gq_vs_int_", TID, ".RDS"))
