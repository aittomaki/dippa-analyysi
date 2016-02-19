
/** simple_priors.stan
 * Simplest possible Stan regression model
 * for master's thesis
 */

data {

	int<lower=0> n; // number of samples
	int<lower=0> d; // number of microRNAs in model
	vector[n] P; // protein expression
	matrix[n,d] M; // microRNA expression
	vector[n] G; // mRNA expression

}

parameters {

	// weights for microRNA and mRNA
	vector[d] beta;
	real gamma;

	// intercept and noise std
	real alpha;
	real<lower=0> sigma;

}

model {

	P ~ normal(M*beta + G*gamma + alpha, sigma);

}
