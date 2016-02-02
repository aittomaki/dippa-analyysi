data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of microRNAs in model
  vector[N] P; // protein expression
  matrix[N,J] M; // microRNA expression
  vector[N] G; // mRNA expression
}
parameters {
  real alpha; // intercept
  vector[J] beta; // microRNA coefs
  real gamma; // mRNA coef
  real<lower=0> sigma;
}
model {
  P ~ normal(alpha + M * beta + gamma * G, sigma);
}
