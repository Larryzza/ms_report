
functions {
  real switch_eta(real eta, real week_index) {
    return(eta + (1 - eta) / (1 + exp(week_index - 6)));
  }
  
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      
      real V = y[1];
      real V3 = y[2];
      real V41 = y[3];
      real V42 = y[4];
      real S = y[5];
      real E = y[6];
      real I = y[7];
      real R = y[8];
      
      real N = x_i[1];
      real beta = theta[1];
      real gamma = theta[2];
      real a = theta[3];
      //real alpha1 = 0.0001;
      //real alpha2 = 0.0001;
      //real alpha3 = 0.0001;
      //real alpha4 = 0.0001;
      real eta = theta[4];
      real week_index = theta[5];
      real epsilon1 = theta[6];
      real epsilon2 = theta[7];
      real epsilon3 = theta[8];
      real epsilon4 = theta[9];
      
      real forcing_function = switch_eta(eta, week_index); // switch function
      real beta_eff = beta * forcing_function;
      
      real dV_dt = -0.0001 * V - (1 - epsilon1) * beta_eff * I * V / N;
      real dV3_dt = -0.0001 * V3 - (1 - epsilon2) * beta_eff * I * V3 / N;
      
      real dV41_dt = -0.0001 * V41 - (1 - epsilon3) * beta_eff * I * V41 / N;
      real dV42_dt = -0.0001 * V42 - (1 - epsilon4) * beta_eff * I * V42 / N;
      
      real dS_dt = 0.0001 * V + 0.0001 * V3 +
                   0.0001 * V41 + 0.0001 * V42 - 
                   beta_eff * I * S / N;
      
      real dE_dt = (1 - epsilon1) * beta_eff * I * V / N +
                   (1 - epsilon2) * beta_eff * I * V3 / N +
                   (1 - epsilon3) * beta_eff * I * V41 / N + 
                   (1 - epsilon4) * beta_eff * I * V42 / N +
                   beta_eff * I * S / N - a * E;
      
      real dI_dt =  a * E - gamma * I;
      
      real dR_dt =  gamma * I;
      
      return {dV_dt, dV3_dt, dV41_dt, dV42_dt, 
              dS_dt, dE_dt, dI_dt, dR_dt};
  }
}

data {
  int<lower=1> n_days;
  int<lower=1> n_weeks;
  real y0[8];
  real ts[7];
  int N;
  int cases[n_days];
  real vac_num[n_weeks, 5];
}

transformed data {
  real x_r[0];
  int x_i[1] = { N };
}

parameters {
  //real<lower=0> gamma;
  real<lower=0> beta;
  //real<lower=0> a;
  //real<lower=0> phi_inv;
  real<lower=0, upper=1> eta;
  real<lower=0, upper=1> epsilon1;
  real<lower=0, upper=1> epsilon2;
  real<lower=0, upper=1> epsilon3;
  real<lower=0, upper=1> epsilon4;
}

transformed parameters{
  vector[n_days] y_out;
  real temp[7,8];
  real gamma = 1./5;
  real phi = 1/exp(-0.6);
  real theta[9] = {beta, gamma, 1./3, eta, 1.0,
                   epsilon1, epsilon2, epsilon3, epsilon4};

  temp = integrate_ode_rk45(sir, y0, 0.0, ts, theta, x_r, x_i);
  y_out[1:7] = col(to_matrix(temp), 7);
  temp[7, 1] = temp[7, 1] + vac_num[1, 1];                                 
  temp[7, 2] = temp[7, 2] + vac_num[1, 2];
  temp[7, 3] = temp[7, 3] + vac_num[1, 3];                                 
  temp[7, 4] = temp[7, 4] + vac_num[1, 4];
  temp[7, 5] = temp[7, 5] - vac_num[1, 5];
  
  for(n in 2:n_weeks){
    theta[5] = n;
    temp = integrate_ode_rk45(sir, temp[7], 0.0, ts, theta, x_r, x_i);
    y_out[(7*n-6):(7*n)] = col(to_matrix(temp), 7);
    temp[7, 1] = temp[7, 1] + vac_num[n, 1];                                 
    temp[7, 2] = temp[7, 2] + vac_num[n, 2];
    temp[7, 3] = temp[7, 3] + vac_num[n, 3];                                 
    temp[7, 4] = temp[7, 4] + vac_num[n, 4];
    temp[7, 5] = temp[7, 5] - vac_num[n, 5];
  }
}

model {
  //priors
  beta ~ normal(2, 1);
  //gamma ~ normal(0.5, 0.3);
  //a ~ normal(0.3, 0.3) T[0,1];
  //phi_inv ~ exponential(5);
  eta ~ beta(2, 4);
  epsilon1 ~ normal(0.3, 0.2);
  epsilon2 ~ normal(0.6, 0.2) T[epsilon1,1];
  epsilon3 ~ normal(0.7, 0.2) T[epsilon2,1];
  epsilon4 ~ normal(0.7, 0.2) T[epsilon2,1];
  //sampling distribution
  for(i in 22:93){
    target += neg_binomial_2_lpmf(cases[i] | y_out[i], phi);
    //target += poisson_lpmf(cases[i] | y_out[i]);
  } 
}

generated quantities {
  real R0 = beta / gamma;
  real pred_cases[n_days];
  for(i in 1:n_days){
    pred_cases[i] = neg_binomial_2_rng(y_out[i], phi);
    //pred_cases[i] = poisson_rng(y_out[i]);
  } 
}

