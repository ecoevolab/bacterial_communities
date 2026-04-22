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

pruebafunc <- stan_ccfunct(df = indv_gro, temp_col = "temp", replica_col = "rep", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           time_series = seq(0, 18, by = 2), time_alternative = c(0, 10, 12, 14, 16, 18), 
                           niterations = 3000, nchains = 4)
rk_valslist <- saveRDS(pruebafunc, file = "03_Output/rk_valslist")















stan_ccfunct <- function (df, temp_col, replica_col, strain_col, sample_byh, interest_col, time_series, time_alternative, niterations, nchains){
  
  # assigning objects to specific values in the data.frame 
  
  spps <- unique(df[[strain_col]]) 
  ntemps <- sort(unique(df[[temp_col]]), decreasing = FALSE)

  ntemps_numeric <- as.numeric(ntemps) 
  ntemps_character <- as.character(ntemps)
  
  vector_freplica <- list() 
  p <- 1
  
  
  for (m in 1:length(spps)) {
    for (o in 1:length(ntemps_numeric)) {
      
      df_complete <- df[df[[strain_col]] == spps[m] & df[[temp_col]] == ntemps_numeric[o], ]
      
      timepoints <- length(unique(df_complete[[sample_byh]]))
      ntotalc <- nrow(df_complete) / timepoints
      
      if (nrow(df_complete) > 0){
        
        df_filtered <- df_complete %>% 
          arrange(.data[[replica_col]], .data[[sample_byh]]) %>%
          pull(.data[[interest_col]])
        
        df_matrix <- matrix(df_filtered, ncol = ntotalc)
        colnames(df_matrix) <- rep(ntemps_numeric[o], times = ntotalc)
        
        vector_freplica[[p]] = df_matrix
        
        names(vector_freplica)[p] <- paste0(spps[m], "_T", ntemps_numeric[o])
        p <- p + 1
      }
    }
  }
  
  
  # Generating initial values and every data.frame for the stan input 
  init_v <- list()
  m_data <- list()
  stan_input <- list()
  ts_vector <- list()
  
  for (q in seq_along(vector_freplica)) {
    
    # extract the first row of every data.frame to get the initial values for stan
    init_v <- as.numeric(vector_freplica[[q]][1, ])
    
    # Extract the rest of the data 
    m_data[[q]] <- vector_freplica[[q]][-1, , drop = FALSE]
    
    # if - to identify specific variations in the data.frame for the time_series 
    
    if (nrow(vector_freplica[[q]]) == length(time_alternative)) {
      
      ts_vector <- time_alternative[-1] 
      
    } else {
      
      ts_vector <- time_series[-1]
      
    }
    
    # That way it doesn't matter if we have another time_series alternative 
    
    # Create the stan input 
    stan_input[[names(vector_freplica)[q]]] <- list(
      N = nrow(m_data[[q]]), 
      S = ncol(m_data[[q]]),
      ts = ts_vector, 
      y0 = init_v,
      y  = m_data[[q]]
    )
  }
  
  stan_output <- list()
  
  for (r in seq_along(stan_input)){
    
    current_S <- stan_input[[r]]$S
    init_fun <- function() {
      list(
        r = 0.3,
        k = 0.8, 
        z0 = as.array(rep(0.001, times = current_S)), # Dinámico según S
        sigma = 0.1
      )
    }
    
    stan_output[[names(stan_input)[r]]] <- stan(
                             model_code = loglv, # here it is the stan function i created earlier
                             data = stan_input[[r]], # stan_input 
                             save_dso = FALSE, 
                             iter = niterations,  # iterations
                             chains = nchains,   # n. chains 
                             init = init_fun) # this function is available in the "Functions" script
    
  }
  return(stan_output)
}

# we're going to use the data.frame without CH29 because the sampling is different from the other spps 

indv_withoutCH29 <- indv_gro %>% filter(Cepa != "CH29")
pruebafunc <- stan_ccfunct(df = indv_withoutCH29, temp_col = "temp", replica_col = "rep", 
                           strain_col = "Cepa", sample_byh = "hr", interest_col = "OD_real", 
                           time_series = seq(0, 18, by = 2), time_alternative = c(0, 10, 12, 14, 16, 18), 
                           niterations = 3000, nchains = 4)

rk_valslist <- saveRDS(pruebafunc, file = "03_Output/rk_by_temp_and_strain")



