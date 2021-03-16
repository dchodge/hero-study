data {
    int N_gp;      // number of sensitivity groups in gp
    real x_gp[N_gp];     // sensitivity groups
    int N;          // number of samples in data
    int n[N];       // total number of samples
    int y[N];      // number of positive samples
}
transformed data {
    real delta = 1e-9;
}
parameters {
    real<lower = 0> alpha;
    real<lower = 0> rho;
    vector[N_gp] eta;
}
transformed parameters {
    vector[N_gp] f;
    vector[N_gp] p;
    matrix[N_gp, N_gp] K;
    matrix[N_gp, N_gp] L_K;

    profile("gp_calc") { 
        K = gp_exp_quad_cov(x_gp, alpha, rho) + diag_matrix(rep_vector(delta, N_gp));
        L_K = cholesky_decompose(K);
        f = L_K * eta;
    }    
    p = inv_logit(f);
}
model {
    profile("prior") { 
        N
        alpha ~ normal(0, 1);
        eta ~ std_normal();    
    }

    y ~ binomial(n, p);
}
generated quantities {
    vector[N_gp] p_post;
    p_post = p;
}
