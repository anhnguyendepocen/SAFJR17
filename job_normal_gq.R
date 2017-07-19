cat("Starting time:", paste(Sys.time()))

# Job script for submission to HPC
rm(list = ls())
gc()
source("functions.R")

if (!requireNamespace("pacman")) install.packages("pacman")
pacman::p_load("doSNOW", "snow", "Rmpi", "dplyr", "foreach")

# Get array ID
TID = Sys.getenv("PBS_ARRAYID")

# Load data
data = readRDS(paste0("Data/simdata_normal_gq_", TID, ".RDS"))

# Setup MPI cluster
cl <- getMPIcluster()
registerDoSNOW(cl)

# Run replications in batches
s = foreach(i = 1:length(data), .combine = rbind) %dopar% {
  out = sim_normal_gq_vs_int(data[[i]])
  return(out)}

# Save the results
saveRDS(s, paste0("Results/res_normal_gq_", TID, ".RDS"))

# Write OK to signal everything went well
cat("TID", TID, "OK. Finishing time:", paste(Sys.time()))

# Finally, stop the parallel cluster
stopCluster(cl)
