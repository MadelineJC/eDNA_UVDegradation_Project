functions{
  vector exp_growth(real t, // Time
                    vector y, // State
                    real r){ // r parameter
    vector[1] dndt;
    dndt[1] = -r*y[1]; 
    return dndt;
  }
}
data {
  int<lower=1> n_times; // Number of time steps minus one
  int<lower=0> y; // Observed data, minus initial state
  int<lower=0> y0_obs; // Observed initial state
  real t0; // First time step
  real ts; // Time steps
}
parameters {
  vector[1] y0; // Initial state
  real r; // Per-capita rate of population change
}
transformed parameters{
    // ODE solver
  vector[1] lambda[n_times] = ode_rk45(exp_growth, y0, t0, ts, r);
}
model {
  // Priors
  r ~ normal(0, 1);
  y0 ~ exponential(1);
  
  // Likelihood
  y0_obs ~ poisson(y0[1]);
  for (t in 1:n_times) {
    y[t] ~ poisson(lambda[t]);
  }
}
generated quantities{
  array[n_times+1] int n_sim;
  n_sim[1] = poisson_rng(y0[1]);
  for(t in 1:n_times){
    n_sim[t+1] = poisson_rng(lambda[t,1]);
  }
}
