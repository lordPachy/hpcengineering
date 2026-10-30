
data {
    int<lower=0> J;
    int<lower=0> N;
    array[N] real y;
    array[N] int<lower=1, upper=J> g;

    real<lower=0> s0;
    real<lower=0> taumax;
    real mu0;
}

parameters {
    real mu;
    real<lower=0> tau;
    
    array[J] real theta;
    real<lower=0> sigma;
}


model {
    mu ~ normal(mu0, s0);
    tau ~ uniform(0, taumax);
    sigma ~ uniform(0, 100);
    theta ~ normal(mu, tau);
    for (i in 1:N) {
        y ~ normal(theta[g[i]], sigma);
    }
}

