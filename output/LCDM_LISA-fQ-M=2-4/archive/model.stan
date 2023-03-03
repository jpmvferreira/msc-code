//// LCDM_GW.stan
// stan model file to test ΛCDM, using the luminosity distance

// user defined functions
functions{
  // 1/E(z)
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    // real x: function argument (integration variable)
    // real xc: complement of the function argument to avoid precision loss (not used explicitly)
    // real[] theta: array of the model parameters (e.g.: theta = {mu, sigma})
    // real[] x_r: data values used to evaluate the integral (can be null)
    // real[] x_i: integer data values used to evaluate the integral (can be null)
    real Omega_m = theta[2];
    return ( Omega_m*(1.0+x)^3 + (1.0-Omega_m) )^(-0.5);
  }
}

// declare variables that will hold the data required to the model
// must match the columns for each input .csv file!
// if two file have the same columns then they will stack as one big array.
// the order of the files in the CLI must be the same as the order on which we define each data variable here
// unless we're stacking several measurements with the same column names, those can show up anywhere
data {
  int N1;
  array[N1] real redshift;
  array[N1] real luminosity_distance;
  array[N1] real error;
}

// define constants and transform the data
// only evaluated in the beginning of each chain
transformed data {
  // create null data values to give to integrate_1d because it's required
  array[0] real x_r;
  array[0] int x_i;
}

// declare the model parameters
// these will be sampled and optimized
parameters {
  real<lower=0> h;
  real<lower=0, upper=1> Omega_m;
}

// allows new variables to be defined in terms of data and/or parameters that may be used later
// will be evaluated on each leapfrog step
transformed parameters {
  // compute the luminosity distance
  array[N1] real dL;
  for (i in 1:N1) {
    dL[i] = (1.0 + redshift[i]) * (2.9979/h) * integrate_1d(integrand, 0, redshift[i], {h, Omega_m}, x_r, x_i);
  }
}

// likelihood and priors
// will be evaluated on each leapfrog step
model {
  // priors
  h ~ normal(0.7, 10);
  Omega_m ~ normal(0.284, 10);

  // likelihood for the GW
  luminosity_distance ~ normal(dL, error);
}

// allows derived quantities based on parameters, data, and rng
// is executed once per iteration
generated quantities {
  // output the natural logarithm of the likelihood to use in model selection criteria
  vector[N1] log_lik;
  for (n in 1:N1) {
    log_lik[n] = normal_lpdf(luminosity_distance[n] | dL[n], error[n]);
  }
}
