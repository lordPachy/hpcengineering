
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
    
    vector[J] theta_raw;
    real<lower=0> sigma;
}

transformed parameters {
    vector[J] theta;
    theta = mu + theta_raw * tau;
}


model {
    mu ~ normal(mu0, s0);
    tau ~ uniform(0, taumax);
    sigma ~ uniform(0, 100);
    theta_raw ~ normal(0, 1);
    for (i in 1:N) {
        y ~ normal(theta[g[i]], sigma);
    }
}

