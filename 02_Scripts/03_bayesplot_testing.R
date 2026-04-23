# Author: Mayra Beatriz Mendoza Velazquez 
# Title: Rstan 

library (bayesplot)
library(ggplot2)

# To do list: 
# * bayesplot - pairs 
# * bayesplot - trace 
# * distribution graph 
# r x iter graph 

# Load the data (stan outputs) ---------------------------------------------------
stan_out_allbut <- readRDS("03_Output/rk_by_temp_and_strain_withoutCH29")
CH29_30_10 <- readRDS("03_Output/CH29_T30_10times")
CH29_30_6 <- readRDS("03_Output/CH29_T30_6times")
CH29_37_10 <- readRDS("03_Output/CH29_T37_10times")
CH29_37_6 <- readRDS("03_Output/CH29_T37_6times")
CH29_42_10 <- readRDS("03_Output/CH29_T42_10times")
CH29_42_6 <- readRDS("03_Output/CH29_T42_6times")

# Merge two lists (for posterior manipulation) -------------------------------------
CH29_allist <- list( CH29_T30_10t = CH29_30_10,
                     CH29_T30_6t = CH29_30_6, 
                     CH29_T37_10t = CH29_37_10,
                     CH29_T37_6t = CH29_37_6,
                     CH29_T42_10t = CH29_42_10,
                     CH29_T42_6t = CH29_42_6)

stan_toutputs <- append(stan_out_allbut, CH29_allist)

# r (& k) * iter graph -------------------------------------------------------------------
# Ref: https://mc-stan.org/rstan/reference/stanfit-method-plot.html


# mcmc_trace ----------------------------------------------------------------
pdf("03_Output/Traceplots_CCStrains.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]
  
  pst_cp <- as.array(stan_toutputs[[i]])
  
  color_scheme_set("mix-brightblue-gray")
  
  # adding title to each plot 
  f <- mcmc_trace(pst_cp, pars = c("r", "k", "sigma")) +
    labs(title = paste("Traceplot:", strain),
         subtitle = "Sampling (1500)")
  
  print(f)
}

dev.off()


# plot()----------------------------------------------
pdf("03_Output/Traceplots_plotfun.pdf", width = 15, height = 7)

for (i in 1:length(stan_toutputs)){
  # Strain for title 
  strain <- names(stan_toutputs)[i]

  # adding title to each plot 
  x <- plot(stan_toutputs[[i]], plotfun = "trace", pars = c("r", "k", "sigma"), inc_warmup = TRUE) +
    labs(title = paste("Traceplot:", strain),
         subtitle = "Shaded area: Warmup / Unshaded area: Sampling")
  
  print(x)
}

dev.off()
