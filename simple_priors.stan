data {
  int<lower=0> N; // number of samples
  int<lower=0> J; // number of microRNAs in model
  vector[N] P; // protein expression
  matrix[N,J] M; // microRNA expression
  vector[N] G; // mRNA expression
}
parameters {
  vector[J] beta; // microRNA coefs
  real gamma; // mRNA coef
  real alpha; // intercept
  real<lower=0> sigma; // error term of regression
}
model {
  P ~ normal(M * beta + gamma * G + alpha, sigma);
}
