
    data {
        int<lower=0> dim;
        int<lower=0> N;
        array[N] vector[dim] x;
        
        real<lower=0> nu0;
        cov_matrix[dim] psi0;
    }
    
    parameters {
        cov_matrix[dim] sigma;
    }
    
    transformed parameters {
        matrix[dim, dim] sigma_chol = cholesky_decompose(sigma);
    }
    
    model {
        vector[dim] mu = rep_vector(0, dim);
        sigma ~ inv_wishart(nu0, psi0);
        x ~ multi_normal_cholesky(mu, sigma_chol);
    }

