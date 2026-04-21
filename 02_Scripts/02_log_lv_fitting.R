# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Logistic Lotka-Volterra fitting 

library(rstan)
library(readr)
library(dplyr)

# Reference --------------------------------------------------------
# Stan & Lotka-Volterra 
# https://canada1.discourse-cdn.com/flex030/uploads/mc_stan/original/2X/b/b6acad3d9f2cd8e7b925c90809f9fb4e4b005a4f.pdf 
 

# Data.frame --------------------------------------------------------
indv_gr <- read_tsv("01_RawData/individual_strains_growth_curves_filtered.tsv")
indv_gr <- as.data.frame(indv_gr)

# data.frame with aditional column for replica
indv_gro <- read_tsv("01_RawData/modified_individual_strain_growth_curves_ord_replica - individual_strains_growth_curves_filtered.tsv")
indv_gro <- as.data.frame(indv_gro)

# Stan function --------------------------------------------------

# Corrections 
# Theta, x_i & x_r are no needed 
# vectors instead of real[] objects 

# Function testing ------------------------------------------------------------
writeLines(loglv, con = "loglv.stan")

# time series 
time <- seq(0, 18, by = 2)

# OD
rep1_CH23 <- subset(indv_gr, Cepa == "CH23" &
                      ((rep == 1 & temp == 30))) %>%
             arrange(temp) %>%
             pull(OD_real)

rep2_CH23 <- subset(indv_gr, Cepa == "CH23" &
                      ((rep == 2 & temp == 30))) %>%
             arrange(temp) %>%
             pull(OD_real)

rep3_CH23 <- subset(indv_gr, Cepa == "CH23" &
                      ((rep == 4 & temp == 30))) %>%
             arrange(temp) %>%
             pull(OD_real)

# OD [ column are the different replicas]
OD_CH23 <- cbind(
  # replica 1
  rep1_CH23,
  # replica 2
  rep2_CH23,
  # replica 3
  rep3_CH23
)

# initial values 
in_val <- c(rep1_CH23[1], rep2_CH23[1], rep3_CH23[1])
N_val <- length(time)

log_CH23_df <- list(
  
    N = 9,
    ts = time[-1],
    y0 = in_val,
    y = OD_CH23[-1,]
  
)

log_CH23fit <- stan(model_code = loglv, 
                    data = log_CH23_df, 
                    save_dso = FALSE, 
                    iter = 2000,  # iterations
                    chains = 4,   # n. chains 
                    init = stan_initial_standarized ) # the initial values i want stan to start from 
print(log_CH23fit)

# CH23 - 30° results stan function --------------------------------------------
# r = 0.30 
# k = 0.8

# Table interpretation notes --------------------------------------------------
# r & k are the ones that i'm interested in 
# r - fitness 
# k - carrying capacity 
# mean 
# se_mean 
# sd - standard deviation 
# % -> quantiles, confidence intervals ?? these provide an estimate of range where 
# the true value is likely to be 
# n_eff -> number of effective samples / if it is too low the aproximations 
# may lack precision, we want this number to be specially high 
# Rhat - Do the mcmc chains get to the same value? - 1.01 is a good indicator
# if this values goes above from 1.01 or 1.05 then that means the value of interest (r, k or the initial one)
# are not realiable values and we cannot use them in the model 
# Values close to 1.0 indicate that the chains have successfully converged to the same distribution

# Pipeline -------------------------------------------------------------------

pruebafunc <- stan_ccfunct(df = indv_gro, temp_col = "temp", replica_col = "ord_replica", strain_col = "Cepa", interest_col = "OD_real", 
                           time_series = seq(0, 18, by = 2), time_alternative = c(0, 10, 12, 14, 16, 18), niterations = 3000, nchains = 4)
rk_valslist <- saveRDS(pruebafunc, file = "03_Output/rk_valslist")

