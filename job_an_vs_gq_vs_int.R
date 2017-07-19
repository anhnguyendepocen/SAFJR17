cat("Starting time:", paste(Sys.time()))

# Job script for submission to HPC
rm(list = ls())
gc()
source("functions.R")

if (!requireNamespace("pacman")) install.packages("pacman")
pacman::p_load("doSNOW", "snow", "Rmpi", "dplyr", "foreach")

# Get number of procs requested
npr = as.numeric(Sys.getenv("PBS_NUM_NODES")) * as.numeric(Sys.getenv("PBS_NUM_PPN"))

# Get array ID
TID = Sys.getenv("PBS_ARRAYID")

# Load data
data = readRDS(paste0("Data/simdata_an_vs_gq_", TID, ".RDS"))

# Setup MPI cluster
cl <- getMPIcluster()
registerDoSNOW(cl)

# Export data and functions to the nodes
clusterExport(cl, "data")
clusterExport(cl, "sim_an_vs_gq_vs_int")

# Split replications in batches
indx = split(1:length(data), ceiling(seq_along(1:length(data)) / (npr - 1)))

# Start timer
tstart = Sys.time()
tstart

# Run replications in batches
s_list = lapply(indx,
                function(xx) {
                  cat("Batch:", xx, "\n")
                  # Run sim_an_vs_gq_vs_int on each simulated dataset
                  ss = parLapply(cl, xx, function(xx) {
                    out = sim_an_vs_gq_vs_int(data[[xx]], xx)
                    return(out)
                  })
                  # Bind rows here
                  ss = bind_rows(ss)
                  return(ss)
                })

# Stop timer
tstop = Sys.time()
tstop

# Bind rows and calculate timediff
s = bind_rows(s_list) %>%
  mutate(timediff = difftime(tstop, tstart, units = "hours") %>% as.numeric())

# Save the results
saveRDS(s, paste0("Results/res_an_vs_gq_", TID, ".RDS"))

# Write OK to signal everything went well
cat("TID", TID, "OK. Finishing time:", paste(Sys.time()))

# Finally, stop the parallel cluster
stopCluster(cl)
