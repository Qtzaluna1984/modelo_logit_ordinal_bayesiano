functions {
  real int_power(real x, int n) {
    real out = 1;
    for (k in 1:n) {
      out *= x;
    }
    return out;
  }

  real ar1_transition_sd(real rho, real sigma_state, int delta) {
    real rho_delta = int_power(rho, delta);
    return sigma_state * sqrt((1 - square(rho_delta)) / (1 - square(rho)));
  }
}

data {
  int<lower=1> N_train;
  int<lower=1> N_test;
  int<lower=1> I;
  int<lower=1> T_train;
  int<lower=1> T_test;
  int<lower=1> K;
  int<lower=2> J;

  matrix[N_train, K] X_train;
  array[N_train] int<lower=1, upper=J> y_train;
  array[N_train] int<lower=1, upper=I> cell_train;
  array[N_train] int<lower=1, upper=T_train> time_train;
  array[N_train] int<lower=1, upper=12> month_train;
  array[T_train] int<lower=1> train_delta;

  matrix[N_test, K] X_test;
  array[N_test] int<lower=1, upper=J> y_test;
  array[N_test] int<lower=1, upper=I> cell_test;
  array[N_test] int<lower=1, upper=T_test> time_test;
  array[N_test] int<lower=1, upper=12> month_test;
  array[T_test] int<lower=1> test_delta;
}

parameters {
  vector[K] beta;
  ordered[J - 1] cutpoints;

  real<lower=-0.99, upper=0.99> rho;
  real<lower=0> sigma_state;
  matrix[I, T_train] state_raw;

  real<lower=0> sigma_cell;
  vector[I] cell_raw;

  real<lower=0> sigma_month;
  vector[12] month_raw;
}

transformed parameters {
  matrix[I, T_train] state;
  vector[I] cell_effect;
  vector[12] month_effect;

  cell_effect = sigma_cell * (cell_raw - mean(cell_raw));
  month_effect = sigma_month * (month_raw - mean(month_raw));

  for (i in 1:I) {
    state[i, 1] = sigma_state * state_raw[i, 1] / sqrt(1 - square(rho));

    for (t in 2:T_train) {
      real rho_delta = int_power(rho, train_delta[t]);
      real transition_sd = ar1_transition_sd(rho, sigma_state, train_delta[t]);
      state[i, t] = rho_delta * state[i, t - 1] + transition_sd * state_raw[i, t];
    }
  }
}

model {
  vector[N_train] eta_train;

  beta ~ normal(0, 0.7);
  cutpoints ~ normal(0, 3);

  rho ~ normal(0, 0.5);
  sigma_state ~ exponential(1);
  to_vector(state_raw) ~ std_normal();

  sigma_cell ~ exponential(1);
  cell_raw ~ std_normal();

  sigma_month ~ exponential(1);
  month_raw ~ std_normal();

  for (n in 1:N_train) {
    eta_train[n] =
      X_train[n] * beta +
      cell_effect[cell_train[n]] +
      month_effect[month_train[n]] +
      state[cell_train[n], time_train[n]];
  }

  y_train ~ ordered_logistic(eta_train, cutpoints);
}

generated quantities {
  vector[N_train] log_lik_train;
  array[N_train] int y_pred_train;
  vector[N_test] log_lik_test;
  array[N_test] int y_pred_test;

  for (n in 1:N_train) {
    real eta =
      X_train[n] * beta +
      cell_effect[cell_train[n]] +
      month_effect[month_train[n]] +
      state[cell_train[n], time_train[n]];

    log_lik_train[n] = ordered_logistic_lpmf(y_train[n] | eta, cutpoints);
    y_pred_train[n] = ordered_logistic_rng(eta, cutpoints);
  }

  {
    matrix[I, T_test] state_test;

    for (i in 1:I) {
      real previous_state = state[i, T_train];

      for (t in 1:T_test) {
        real rho_delta = int_power(rho, test_delta[t]);
        real transition_sd = ar1_transition_sd(rho, sigma_state, test_delta[t]);
        state_test[i, t] = rho_delta * previous_state + transition_sd * normal_rng(0, 1);
        previous_state = state_test[i, t];
      }
    }

    for (n in 1:N_test) {
      real eta =
        X_test[n] * beta +
        cell_effect[cell_test[n]] +
        month_effect[month_test[n]] +
        state_test[cell_test[n], time_test[n]];

      log_lik_test[n] = ordered_logistic_lpmf(y_test[n] | eta, cutpoints);
      y_pred_test[n] = ordered_logistic_rng(eta, cutpoints);
    }
  }
}
