//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  real switch_eta(real eta, real week_index, real xi) {
    return(eta + (1 - eta) / (1 + exp(xi * (week_index - 5))));
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
      real alpha1 = 0.0001;
      real alpha2 = 0.0001;
      real alpha3 = 0.0001;
      real alpha4 = 0.0001;
      real epsilon1 = 0.8;
      real epsilon2 = 1;
      real epsilon3 = 1;
      real epsilon4 = 1;
      real eta = theta[4];
      real xi = theta[5];
      real week_index = theta[6];
      
      real forcing_function = switch_eta(eta, week_index, xi); // switch function
      real beta_eff = beta * forcing_function;
      
      real dV_dt = -alpha1 * V - (1 - epsilon1) * beta_eff * I * V / N;
      real dV3_dt = -alpha2 * V3 - (1 - epsilon2) * beta_eff * I * V3 / N;
      
      real dV41_dt = -alpha3 * V41 - (1 - epsilon3) * beta_eff * I * V41 / N;
      real dV42_dt = -alpha4 * V42 - (1 - epsilon4) * beta_eff * I * V42 / N;
      
      real dS_dt = alpha1 * V + alpha2 * V3 +
                   alpha3 * V41 + alpha4 * V42 - 
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
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> a;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> eta;
  real<lower=0,upper=1> xi_raw;
}

transformed parameters{
  vector[n_days] y_out;
  real temp[7,8];
  real phi = 1. / phi_inv;
  real xi = xi_raw + 0.5;
  real theta[6] = {beta, gamma, a, eta, xi,1.0};
  
  temp = integrate_ode_rk45(sir, y0, 0.0, ts, theta, x_r, x_i);
  y_out[1:7] = col(to_matrix(temp), 7);
  temp[7, 1] = temp[7, 1] + vac_num[1, 1];                                 
  temp[7, 2] = temp[7, 2] + vac_num[1, 2];
  temp[7, 3] = temp[7, 3] + vac_num[1, 3];                                 
  temp[7, 4] = temp[7, 4] + vac_num[1, 4];
  temp[7, 5] = temp[7, 5] - vac_num[1, 5];
  
  for(n in 2:n_weeks){
    theta[6] = n;
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
  beta ~ normal(3, 2);
  gamma ~ normal(0.3, 0.5);
  a ~ normal(0.3, 0.5);
  phi_inv ~ exponential(5);
  eta ~ beta(2.5, 4);
  xi_raw ~ beta(1, 1);
  //sampling distribution
  cases ~ neg_binomial_2(y_out, phi);
}

generated quantities {
  real R0 = beta / gamma;
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(y_out, phi);
}

