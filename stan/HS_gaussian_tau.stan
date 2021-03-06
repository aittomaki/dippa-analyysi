
/** shrinkage_prior.stan
 * Stan model with hierarchical shrinkage priors
 * for master's thesis.
 * Modified from lg_t.stan by Juho Piironen.
 */

functions {
	// square root of a vector (elementwise)
	vector sqrt_vec(vector x) {
		vector[dims(x)[1]] res;

		for (m in 1:dims(x)[1]){
			res[m] <- sqrt(x[m]);
		}
		return res;
	}
}

data {
	int<lower=0> n; // number of samples
	int<lower=0> d; // number of microRNA
	vector[n] P;	// protein expr output
	matrix[n,d] M;	// microRNA expr input
	vector[n] G;    // gene expr iput
	real<lower=1> nu; // degrees of freedom for the half t-priors
	real<lower=0> pn;// assumed number of meaningful variables
}

parameters {

	// intercept and noise std
	real w0;
	real<lower=0> sigma;

	// gene expression weight
	real wg;

	// auxiliary variables for the variance parameters
	vector[d] z;
	//real<lower=0> r1_global;
	//real<lower=0> r2_global;
    real<lower=0> tau;
	vector<lower=0>[d] r1_local;
	vector<lower=0>[d] r2_local;
}

transformed parameters {

	// global and local variance parameters, and the input weights
	//real<lower=0> tau;
	vector<lower=0>[d] lambda;
	vector[d] w;

	//tau <- r1_global * sqrt(r2_global); // this is a student's T!
	lambda <- r1_local .* sqrt_vec(r2_local); // this is a student's T!
	w <- z .* lambda*tau;
}

model {

	// observation model
	P ~ normal(M*w + G*wg + w0, sigma);

	// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
	z ~ normal(0, 1);
	r1_local ~ normal(0.0, 1.0);
	r2_local ~ inv_gamma(0.5*nu, 0.5*nu);

	// tight half cauchy for tau
	//r1_global ~ normal(0.0, pn/n);
	//r2_global ~ inv_gamma(0.5, 0.5);
    tau ~ normal(0.0, pn/n);

	// weakly informative prior for the intercept and gene weight
	w0 ~ normal(0,5);
	wg ~ normal(0,5);

	// using uniform prior on the noise variance
}
