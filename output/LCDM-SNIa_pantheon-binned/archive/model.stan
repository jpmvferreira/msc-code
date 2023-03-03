// Description:
// - implementation of ΛCDM to be constrained using SNIa

// user defined functions
functions{
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    // real x: function argument (integration variable)
    // real xc: complement of the function argument to avoid precision loss (not used explicitly)
    // real[] theta: array of the model parameters (e.g.: theta = {mu, sigma})
    // real[] x_r: data values used to evaluate the integral (can be null)
    // real[] x_i: integer data values used to evaluate the integral (can be null)
    real Omega_m = theta[1];
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
  array[N1] real zcmb;
  array[N1] real mb;
  array[N1] real dmb;
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
  real<lower=0,upper=1> Omega_m;
}

// allows new variables to be defined in terms of data and/or parameters that may be used later
// will be evaluated on each leapfrog step
transformed parameters {
  real Delta;
  array[N1] real a;
  array[N1] real b;
  array[N1] real c;

  for (i in 1:N1) {
    Delta = mb[i] - 5.0*log10((1.0+zcmb[i]) * integrate_1d(integrand, 0, zcmb[i], {Omega_m}, x_r, x_i));
    a[i] = (Delta/dmb[i])^2;
    b[i] = Delta/dmb[i]^2;
    c[i] = dmb[i]^(-2);
  }

  real A = sum(a);
  real B = sum(b);
  real C = sum(c);
}

// likelihood and priors
// will be evaluated on each leapfrog step
model {
  // priors
  Omega_m ~ normal(0.284, 10);

  // likelihood for the SNIa
  // increment the log posterior probability
  // already removed factors which do not depend on the parameters
  target += -0.5*(A - B^2/C);
}

// block for generated quantites
// allows derived quantities based on parameters, data, and rng
// is executed once per iteration
generated quantities {
  // output the natural logarithm of the likelihood to use in model selection criteria
  vector[N1] log_lik;
  for (i in 1:N1) {
    log_lik[i] = - 0.5*(a[i] + 2*(-B/C)*b[i] + c[i]*(B/C)^2);
  }
}
