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
View(indv_gr)


# Stan function --------------------------------------------------

loglv <- '

functions {
 real[] lvfnc(real t,       // time
               real[] z,                     // x values [for each replica]
               real[] theta,                 // parameters [r & k]
               real[] x_r,              // real data 
               int[] x_i){              // integer data
              
    real dzdt[3];
    real r = theta[1];
    real k = theta[2];
    
      for (j in 1:3) {
        dzdt[j] = r * z[j] * (1 - z[j] / k);
    }

  return dzdt;
}
}

               
  data {
    int <lower=1> N;  // total num measurements  
    real ts[N];      // measurement times > 0 
    real <lower=0> y0[3];    // initial measured population [3 for each temp measure]
    real <lower=0> y[N, 3]; // measured population at measurement times 
   }
  
  parameters {
  real<lower=0> r;
  real<lower=0> k;
  array[3] real<lower=0> z0;
  real<lower=0> sigma;
}

  
  transformed parameters {
     real theta[2];
     real z[N,3];

     theta[1] = r;
     theta[2] = k;

         z = integrate_ode_rk45(lvfnc,
                                z0,
                                0.0,
                                ts,
                                theta,
                                rep_array(0.0,0),
                                rep_array(0,0)
                                );

                           }
  
  
 model {
 
  // priors
  r ~ normal(1, 0.5);
  k ~ normal(1, 0.5);
  sigma ~ lognormal(0, 0.5); 


  // initial state
   for (j in 1:3) {
    z0[j] ~ normal(y0[j], 0.8);
   }
  
  
  // likelihood (lognormal)
  
  for (m in 1:3) {
    y0[m] ~ normal(z0[m], sigma);

    y[ , m] ~ normal(z[, m], sigma);
  }
  
}

'

writeLines(loglv, con = "loglv.stan")

# time series 
time <- seq(0, 18, by = 2)

# od
indv_gr

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

OD_CH23 <- cbind(
  # replica 1
  rep1_CH23,
  # replica 2
  rep2_CH23,
  # replica 3
  rep3_CH23
)

OD_CH23[OD_CH23 <= 0] <- 1e-6


in_val <- c(rep1_CH23[1], rep2_CH23[1], rep3_CH23[1])
N_val <- length(time)

log_CH23_df <- list(
  
    N = N_val,
    ts = time,
    y0 = in_val,
    y = OD_CH23
  
)

CH23_init <- function() {
  list(
    r = 0.3,
    k = 0.8,
    z0 = as.numeric(y0),
    sigma = 0.1
  )
}

log_CH23fit <- stan(model_code = loglv, 
                    data = log_CH23_df,
                    save_dso = FALSE, 
                    iter = 20,
                    chains = 4)
print(exfit)
plot(exfit)

diff()

