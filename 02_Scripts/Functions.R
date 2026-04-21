# Author: Mayra Beatriz Mendoza Velázquez 
# Title: Functions 

# Stan functions ---------------------------------------------------------------

# 1. Initial values standarized 
stan_initial_standarized <- function(){
  list(
    r = 0.3,
    k = 0.8, 
    z0 = as.array(c(0.001, 0.001, 0.001)),
    sigma = 0.1
  )
}

# Reconstruction methods -------------------------------------------------------