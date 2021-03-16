data {
    int N_gp;      // number of sensitivity groups in gp
    real x_gp[N_gp];     // sensitivity groups
    int N;          // number of sensitivtiy groups in data
    int D;          // number of age groups in data
    int n[N, D];       // total number of samples
    int y[N, D];      // number of positive samples
}
transformed data {
    real delta = 1e-9;
}
parameters {
    vector<lower = 0>[D] alpha;
    real<lower = 0> rho;
    cholesky_factor_corr[D] L_Omega;
    matrix[N_gp, D] eta;
}
transformed parameters {
    matrix[N_gp, D] f;
    matrix[N_gp, D] p;
    matrix[N_gp, N_gp] K;
    matrix[N_gp, N_gp] L_K;

    profile("gp_calc") { 
        K = cov_exp_quad(x_gp, 1.0, rho) + diag_matrix(rep_vector(delta, N_gp));
        L_K = cholesky_decompose(K);
        f = L_K * eta * diag_pre_multiply(alpha, L_Omega)';
    }    
    p = inv_logit(f);
    
}
model {
    profile("prior") { 
        rho ~ inv_gamma(5, 5);
        alpha ~ normal(0, 1);
        to_vector(eta) ~ std_normal();    
        L_Omega ~ lkj_corr_cholesky(3);
    }
    for (i in 1:N) {
       for (j in 1:D){
            y[i, j] ~ binomial(n[i, j], p[i, j]);
       }
    }
}
generated quantities {
    matrix[N_gp, D] p_post;
    for (i in 1:N_gp) {
       for (j in 1:D) {
            p_post[i, j] = p[i, j];
       }
    }
}
