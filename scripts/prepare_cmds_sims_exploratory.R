##################################################################
# 
# Prepare commands to run BiSSE simulations.
# 
# spec0 - speciation rate of diploid lineages
# ext0  - extinction rate of diploid lineages
# spec1 - speciation rate of polyploid lineages
# ext1  - extinction rate of polyploid lineages
# q01   - diploid-to-polyploid transition rate
# q10   - polyploid-to-diploid transition rate
# 
# RF means relative fraction
# 
##################################################################

spec0 <- c(1)
rf_ext0 <- c(0.10, 0.50, 0.90) # RF of spec0
rf_spec1 <- c(0.25, 0.50, 1.0, 1.50, 1.75) # RF of spec0
rf_ext1 <- c(0.10, 0.50, 0.90) # RF of spec1
rf_q01 <- c(1/8, 1/4, 1/2) # RF of spec0
rf_q10 <- c(0.25, 0.50, 1.0, 1.50, 1.75) # RF of q01


N <- length(spec0) * length(rf_ext0)
N <- N * length(rf_spec1) * length(rf_ext1)
N <- N * length(rf_q01) * length(rf_q10)

df <- matrix(nrow = N, ncol = 6)
colnames(df) <- c('spec0', 'ext0', 'spec1', 'ext1', 'q01', 'q10')

row_index <- 0
for(i in 1:length(rf_ext0)){
  ext0 <- rf_ext0[i] * spec0
  for(j in 1:length(rf_spec1)){
    spec1 <- rf_spec1[j] * spec0
    for(k in 1:length(rf_ext1)){
      ext1 <- rf_ext1[k] * spec1
      for(l in 1:length(rf_q01)){
        q01 <- rf_q01[l] * spec0
        for(m in 1:length(rf_q10)){
          q10 <- rf_q10[m] * q01
          pars <- c(spec0, ext0, spec1, ext1, q01, q10)
          row_index <- row_index + 1
          df[row_index, ] <- pars
        }
      }
    }
  }
}

stopifnot(N == dim(df)[1])


# Get the unique values of the rate parameters
print("spec0")
unique(df[, "spec0"])
print("ext0")
unique(df[, "ext0"])
print("spec1")
unique(df[, "spec1"])
print("ext1")
unique(df[, "ext1"])
print("q01")
unique(df[, "q01"])
print("q10")
unique(df[, "q10"])


# Print commands
Rscript_exe <- "Rscript"
Rscript_opt <- "--no-save"
script_file <- "simulate_histories.R"
cmd_file <- "sims_cmds_exploratory.txt"


n_taxa <- 10^3
n_reps <- 10^3

cmds_list <- c()

for(i in 1:nrow(df)){
  out_file <- paste0("sim", "_", i, ".csv")
  
  cmd <- paste(Rscript_exe,
               Rscript_opt,
               script_file,
               paste0("--spec0", "=", df[i, "spec0"]),
               paste0("--spec1", "=", df[i, "spec1"]),
               paste0("--ext0", "=", df[i, "ext0"]),
               paste0("--ext1", "=", df[i, "ext1"]),
               paste0("--q01", "=", df[i, "q01"]),
               paste0("--q10", "=", df[i, "q10"]),
               paste0("--ntaxa", "=", n_taxa),
               paste0("--nreps", "=", n_reps),
               paste0("--out_file", "=", out_file),
               sep=" "
               )
  
  cmds_list <- c(cmds_list, cmd)
}


file_conn <- file(cmd_file)
writeLines(cmds_list, file_conn)
close(file_conn)
