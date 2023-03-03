// Description
// - an implementation of the f(Q) model found in arXiv:2104.15123, neglecting radiation, to be constrained using SNIa

/* block for user defined functions */
functions {
  // The arguments to this function are:
  // - `time`: the time to evaluate the DAE system
  // - `state`: the state of the DAE system at the time specified
  // - `state_derivative`: the time derivatives of the state of the DAE system at the time specified
  // - `...`: sequence of arguments passed unmodified from the DAE solve function call.
  vector residual(real z, vector state, vector deriv, real lambda, real Omega_m) {
    real E = state[1];
    real lhs = (E^2 - 2*lambda)*exp(lambda/E^2);  // left hand side of the equation for E
    real rhs = Omega_m*(1+z)^3;  // right hand side of the equation for E

    // compute the residuals
    vector[2] res;
    res[1] = rhs - lhs;       // E satisfies the algebraic relation, meaning that res[1] should be 0
    res[2] = 1/E - deriv[2];  // the derivative of S₂ is in fact 1/E, this should also be 0

    return res;
  }
}

/* block for datasets used
- must match the columns for each input .csv file, and the number of events must be N1, N2, ...
- if two file have the same columns then they will stack as one big array.
- the order of the files in the CLI must be the same as the order on which we define each data variable here, unless we're stacking several measurements with the same column names, those can show up anywhere */
data {
  int N1;
  array[N1] real zcmb;
  array[N1] real mb;
  array[N1] real dmb;
}

/* block for transformed data
- only evaluated in the beginning of each chain */
transformed data {
}

/* block for model parameters */
parameters {
  real<lower=0,upper=1> Omega_m;
}

/* block for transformed parameters
- allows derived quantities from parameters and data
- is evaluated on each leapfrog step */
transformed parameters {
  // quantity derived from Omega_m
  real lambda;
  lambda = 0.5 + lambert_w0( -Omega_m/(2*exp(0.5)) );

  // compute dE/dz(z = 0), knowing that E(0) = 1
  real deriv_0;
  deriv_0 = 1/(2*exp(lambda)) * 3*Omega_m/(1 - lambda + 2*lambda^2);

  // call the differential algebraic equation solver
  // `residual`: DAE residual function,
  // `initial_state`: initial state, type vector
  // `initial_state_derivative`: time derivative of the initial state, type vector
  // `initial_time`: initial time, type real
  // `times`: solution times, type array[] real
  // `...`: sequence of arguments that will be passed through unmodified to the DAE residual function.
  array[N1] vector[2] S;
  S = dae(residual, [1, 0.0]', [deriv_0, 1]', 0.0, zcmb, lambda, Omega_m);

  // compute Δ and use it to compute A, B and C for the marginalized likelihood for the SNIa
  real Delta;
  array[N1] real a;
  array[N1] real b;
  array[N1] real c;
  
  for (i in 1:N1) {
    Delta = mb[i] - 5.0*log10((1.0+zcmb[i]) * S[i,2]);
    a[i] = (Delta/dmb[i])^2;
    b[i] = Delta/dmb[i]^2;
    c[i] = dmb[i]^(-2);
  }

  real A = sum(a);
  real B = sum(b);
  real C = sum(c);
}

/* block for model definition */
model {
  // priors
  Omega_m ~ normal(0.284, 10);

  // increment the log likelihood (aka target) by a given ammount
  target += -0.5*(A - B^2/C);
}

/* block for generated quantites
- allows derived quantities based on parameters, data, and rng
- is executed once per iteration */
generated quantities {
  // output the natural logarithm of the likelihood to use in model selection criteria
  vector[N1] log_lik;
  for (i in 1:N1) {
    log_lik[i] = - 0.5*(a[i] + 2*(-B/C)*b[i] + c[i]*(B/C)^2);
  }
}
