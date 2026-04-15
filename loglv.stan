

functions {
 real[] lvfnc(real t,       // time
               real[] z,                     // x values [for each replica]
               real[] theta,                 // parameters [r & k]
               real[] x_r,              // real data 
               int[] x_i){              // integer data
              
    real dzdt[3]; 
    
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


