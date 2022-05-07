
functions {

  real switch_eta(real eta, real week_index, real k, real day_shift) {
    return(eta + (1 - eta) / (1 + k * exp(week_index - day_shift)));
    //return(eta + (1 - eta) / (1 + exp(week_index - 7)));
  }
  
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      
      real V = y[1];
      real V31 = y[2];
      real V32 = y[3];
      real V41 = y[4];
      real V42 = y[5];
      real S = y[6];
      real E = y[7];
      real EV = y[8];
      real EV1 = y[9];
      real EV2 = y[10];
      real I = y[11];
      real IV = y[12];
      real IV1 = y[13];
      real IV2 = y[14];
      real R = y[15];
      
      real N = x_i[1];
      real beta = theta[1];
      real gamma = theta[2];
      real sigma = theta[3];
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
      real epsilon5 = theta[10];
      real k = theta[11];
      real day_shift = theta[12];
      
      real c = switch_eta(eta, week_index, k, day_shift); // switch function with k
      //real c = switch_eta(eta, week_index); // switch function without k
      
      real dV_dt = -0.0001 * V - (1 - epsilon1) * beta * c * I * V / N;
      real dV31_dt = -0.0001 * V31 - (1 - epsilon2) * beta * c * I * V31 / N;
      real dV32_dt = -0.0001 * V32 - (1 - epsilon3) * beta * c * I * V32 / N;
      real dV41_dt = -0.0001 * V41 - (1 - epsilon4) * beta * c * I * V41 / N;
      real dV42_dt = -0.0001 * V42 - (1 - epsilon5) * beta * c * I * V42 / N;
      
      real dS_dt = 0.0001 * V + 
                   0.0001 * V31 + 0.0001 * V32 +
                   0.0001 * V41 + 0.0001 * V42 - 
                   beta * c * I * S / N;
      
      real dE_dt = beta * c * I * S / N - sigma * E;
      real dEV_dt = (1 - epsilon1) * beta * c * I * V / N  - sigma * EV;
      real dEV1_dt = (1 - epsilon2) * beta * c * I * V31 / N +
                     (1 - epsilon4) * beta * c * I * V41 / N -  
                     sigma * EV1;
      real dEV2_dt = (1 - epsilon3) * beta * c * I * V32 / N +
                     (1 - epsilon5) * beta * c * I * V42 / N -  
                     sigma * EV2;
      
      real dI_dt =  sigma * E - gamma * I;
      real dIV_dt =  sigma * EV - gamma * IV;
      real dIV1_dt =  sigma * EV1 - gamma * IV1;
      real dIV2_dt =  sigma * EV2 - gamma * IV2;
      
      real dR_dt =  gamma * I;
      
      return {dV_dt, dV31_dt, dV32_dt, dV41_dt, dV42_dt, dS_dt, 
              dE_dt, dEV_dt, dEV1_dt, dEV2_dt, 
              dI_dt, dIV_dt, dIV1_dt, dIV2_dt, dR_dt};
  }
  
  real[] col_sums(matrix X) {
     real s[cols(X)] ;
     int col_index[4] = {11,12,13,14};
     for (j in 1:4) s[j] = sum(col(X, col_index[j])) ;
	   return s ;
  }
}

data {
  int<lower=1> n_weeks;
  real y0[15];
  real ts[7];
  int N;
  int cases[n_weeks, 4];
  real vac_num[n_weeks, 6];
}

transformed data {
  real x_r[0];
  int x_i[1] = { N };
}

parameters {
  //real<lower=0> gamma;
  real<lower=0> beta;
  //real<lower=0> sigma;
  //real<lower=0> phi_inv;
  real<lower=0, upper=1> eta;
  real<lower=0, upper=1> epsilon1;
  real<lower=0, upper=1> epsilon2;
  real<lower=0, upper=1> epsilon3;
  real<lower=0, upper=1> epsilon4;
  real<lower=0, upper=1> epsilon5;
  real<lower=0> k;
  real day_shift;
}

transformed parameters{
  real y_out[n_weeks, 4];
  real temp[7,15];
  real gamma = 1./5;
  real phi = 10;
  real theta[12] = {beta, gamma, 1./3, eta, 1.0,
                   epsilon1, epsilon2, epsilon3, epsilon4, epsilon5, 
                   k, day_shift};

  temp = integrate_ode_rk45(sir, y0, 0.0, ts, theta, x_r, x_i);
  y_out[1] = col_sums(to_matrix(temp));
  temp[7, 1] = temp[7, 1] + vac_num[1, 1];                                 
  temp[7, 2] = temp[7, 2] + vac_num[1, 2];
  temp[7, 3] = temp[7, 3] + vac_num[1, 3];                                 
  temp[7, 4] = temp[7, 4] + vac_num[1, 4];
  temp[7, 5] = temp[7, 5] + vac_num[1, 5];
  temp[7, 6] = temp[7, 6] - vac_num[1, 6];
  
  for(n in 2:n_weeks){
    theta[5] = n;
    temp = integrate_ode_rk45(sir, temp[7], 0.0, ts, theta, x_r, x_i);
    y_out[n] = col_sums(to_matrix(temp));
    temp[7, 1] = temp[7, 1] + vac_num[n, 1];                                 
    temp[7, 2] = temp[7, 2] + vac_num[n, 2];
    temp[7, 3] = temp[7, 3] + vac_num[n, 3];                                 
    temp[7, 4] = temp[7, 4] + vac_num[n, 4];
    temp[7, 5] = temp[7, 5] + vac_num[n, 5];
    temp[7, 6] = temp[7, 6] - vac_num[n, 6];
  }
}

model {
  //priors
  beta ~ normal(2, 1) T[0,];
  k ~ normal(1, 2) T[0,];
  day_shift ~ normal(6, 2);
  //gamma ~ normal(0.5, 0.3);
  //sigma ~ normal(0.3, 0.3) T[0,1];
  //phi_inv ~ exponential(5);
  eta ~ beta(2, 4);
  epsilon1 ~ normal(0.3, 0.2) T[0,1];
  epsilon2 ~ normal(0.6, 0.2) T[epsilon1,1];
  epsilon3 ~ normal(0.6, 0.2) T[epsilon1,1];
  epsilon4 ~ normal(0.7, 0.2) T[epsilon2,1];
  epsilon5 ~ normal(0.7, 0.2) T[epsilon3,1];
  //sampling distribution
  for(i in 4:n_weeks){
    for(j in 1:4){
      target += neg_binomial_2_lpmf(cases[i,j] | y_out[i,j], phi);
    //target += poisson_lpmf(cases[i] | y_out[i]);
    }
  } 
}

generated quantities {
  real R0 = beta / gamma;
  real pred_cases[n_weeks, 4];
  for(i in 1:n_weeks){
    for(j in 1:4){
      pred_cases[i,j] = neg_binomial_2_rng(y_out[i,j], phi);
    //pred_cases[i] = poisson_rng(y_out[i]);
    }
  } 
}

