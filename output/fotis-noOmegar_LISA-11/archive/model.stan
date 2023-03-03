/* fotis.stan
- an implementation of the f(Q) model found in arXiv:2104.15123, neglecting radiation, to be constrained using standard sirens
*/

/* block for user defined functions */
functions {
  // The arguments to this function are:
  // - `time`: the time to evaluate the DAE system
  // - `state`: the state of the DAE system at the time specified
  // - `state_derivative`: the time derivatives of the state of the DAE system at the time specified
  // - `...`: sequence of arguments passed unmodified from the DAE solve function call.
  vector residual(real z, vector state, vector deriv, real lambda, real Omega_m) {
    real E = state[1];
    real lhs = (E^2 - 2*lambda)*exp(lambda*E^-2);  // left hand side of the equation for E
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
  array[N1] real redshift;
  array[N1] real luminosity_distance;
  array[N1] real error;
}

/* block for transformed data
- only evaluated in the beginning of each chain */
transformed data {
}

/* block for model parameters */
parameters {
  real<lower=0> h;
  real<upper=1.21306131942526684721> Omega_m;  // constrain such that Ωm ≤ 2e^(-0.5)
}

/* block for transformed parameters
- allows derived quantities from parameters and data
- is evaluated on each leapfrog step */
transformed parameters {
  array[N1] real dGW;   // luminosity distance for gravitational waves
  real correction;      // correction from modification to gravity
  real E;               // E(z)
  real k = 2.9979;      // c/H0 = (k Gpc)/h

  // quantity derived from Omega_m
  real lambda;
  lambda = 0.5 + lambert_w0( -Omega_m/(2*exp(0.5)) );

  // compute dE/dz(z = 0), knowing that E(0) = 1
  real deriv_0;
  deriv_0 = 1/(2*exp(lambda)) * 4*Omega_m/(1 - lambda + 2*lambda^2);

  // first entry is S₁(z) and second is S₂(z), order is also verified below
  array[N1] vector[2] S;

  // call the differential algebraic equation solver
  // `residual`: DAE residual function,
  // `initial_state`: initial state, type vector
  // `initial_state_derivative`: time derivative of the initial state, type vector
  // `initial_time`: initial time, type real
  // `times`: solution times, type array[] real
  // `...`: sequence of arguments that will be passed through unmodified to the DAE residual function.
  S = dae(residual, [1, 0.0]', [deriv_0, 1]', 0.0, redshift, lambda, Omega_m);

  // compute theoretical luminosity distance for gravitational waves for each measured redshift
  for (i in 1:N1) {
    E = S[i,1];
    correction = ((1-lambda)/(1-lambda/E^2))^0.5 * exp(lambda/2 * (E^2 - 1)/E^2);
    dGW[i] = correction * (1.0+redshift[i]) * (k/h) * S[i,2];
  }
}

/* block for model definition */
model {
  // priors
  h ~ normal(0.68, 10);
  Omega_m ~ normal(0.353, 10);

  // likelihood
  luminosity_distance ~ normal(dGW, error);
}

/* block for generated quantites
- allows derived quantities based on parameters, data, and rng
- is executed once per iteration */
generated quantities {
  // output the natural logarithm of the likelihood to use in model selection criteria
  vector[N1] log_lik;
  for (n in 1:N1) {
    log_lik[n] = normal_lpdf(luminosity_distance[n] | dGW[n], error[n]);
  }
}
