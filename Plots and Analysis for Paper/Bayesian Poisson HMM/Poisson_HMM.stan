data {
  int<lower=0> N;
  int<lower=0> K;
  int y[N];
}

parameters {
  simplex[K] theta[K];
  // real mu[K];
  positive_ordered[K] mu;
}

model {
  // priors
  //target+= normal_lpdf(mu[1] | 1, 1);
  //target+= normal_lpdf(mu[2] | 4, 1);
  //target+= normal_lpdf(mu[3] | 12, 1);
  // forward algorithm
  {
  real acc[K];
  real gamma[N, K];
  for (k in 1:K)
    gamma[1, k] = poisson_lpmf(y[1] | mu[k]);
  for (t in 2:N) {
    for (k in 1:K) {
      for (j in 1:K)
        acc[j] = gamma[t-1, j] + log(theta[j, k]) + poisson_lpmf(y[t] | mu[k]);
      gamma[t, k] = log_sum_exp(acc);
    }
  }
  target += log_sum_exp(gamma[N]);
  }
}

generated quantities {
  int<lower=1,upper=K> z_star[N];
  real log_p_z_star;
  {
    int back_ptr[N, K];
    real best_logp[N, K];
    for (k in 1:K)
      best_logp[1, k] = poisson_lpmf(y[1] | mu[k]);
    for (t in 2:N) {
      for (k in 1:K) {
        best_logp[t, k] = negative_infinity();
        for (j in 1:K) {
          real logp;
          logp = best_logp[t-1, j] + log(theta[j, k]) + poisson_lpmf(y[t] | mu[k]);
          if (logp > best_logp[t, k]) {
            back_ptr[t, k] = j;
            best_logp[t, k] = logp;
          }
        }
      }
    }
    log_p_z_star = max(best_logp[N]);
    for (k in 1:K)
      if (best_logp[N, k] == log_p_z_star)
        z_star[N] = k;
    for (t in 1:(N - 1))
      z_star[N - t] = back_ptr[N - t + 1, z_star[N - t + 1]];
  }
}