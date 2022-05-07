
functions {

  real switch_eta(real eta, real week_index, real k, real day_shift) {
    return(eta + (1 - eta) / (1 + k * exp(week_index - day_shift)));
    //return(eta + (1 - eta) / (1 + exp(week_index - 7)));
  }
  
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      
      real V2 = y[1];
      real V3 = y[2];
      real V4 = y[3];
      real S = y[4];
      real E = y[5];
      real EV2 = y[6];
      real EV3 = y[7];
      real I = y[8];
      real IV2 = y[9];
      real IV3 = y[10];
      //real R = y[11];
      
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
      real k = theta[9];
      real day_shift = theta[10];
      
      real c = switch_eta(eta, week_index, k, day_shift); // switch function with k
      //real c = switch_eta(eta, week_index); // switch function without k
      
      real dV2_dt = -0.0001 * V2 - (1 - epsilon1) * beta * c * I * V2 / N;
      real dV3_dt = -0.0001 * V3 - (1 - epsilon2) * beta * c * I * V3 / N;
      real dV4_dt = -0.0001 * V4 - (1 - epsilon3) * beta * c * I * V4 / N;
      
      real dS_dt = 0.0001 * V2 + 0.0001 * V3 +  0.0001 * V4 - 
                   beta * c * I * S / N;
      
      real dE_dt = beta * c * I * S / N - sigma * E;
      real dEV2_dt = (1 - epsilon1) * beta * c * I * V2 / N  - sigma * EV2;
      real dEV3_dt = (1 - epsilon2) * beta * c * I * V3 / N + 
                     (1 - epsilon3) * beta * c * I * V4 / N
                     - sigma * EV3;
                     
      real dI_dt =  sigma * E - gamma * I;
      real dIV2_dt =  sigma * EV2 - gamma * IV2;
      real dIV3_dt =  sigma * EV3 - gamma * IV3;
      
      //real dR_dt =  gamma * I + gamma * IV2 + gamma * IV3;
      
      return {dV2_dt, dV3_dt, dV4_dt, dS_dt, 
              dE_dt, dEV2_dt, dEV3_dt,  
              dI_dt, dIV2_dt, dIV3_dt};
  }
  
  real[] col_sums(matrix X) {
     real s[cols(X)] ;
     int col_index[3] = {8, 9, 10};
     for (j in 1:3) s[j] = sum(col(X, col_index[j])) ;
	   return s ;
  }
}

data {
  int<lower=1> n_weeks;
  real y0[10];
  real ts[7];
  int N;
  int cases[n_weeks, 3];
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
  //real<lower=0, upper=1> epsilon4;
  //real<lower=0, upper=1> epsilon5;
  real<lower=0> k;
  real day_shift;
}

transformed parameters{
  real y_out[n_weeks, 3];
  real temp[7,10];
  real gamma = 1./5;
  real phi = 10;
  real theta[10] = {beta, gamma, 1./3, eta, 1.0,
                   epsilon1, epsilon2, epsilon3, 
                   k, day_shift};

  temp = integrate_ode_rk45(sir, y0, 0.0, ts, theta, x_r, x_i);
  y_out[1] = col_sums(to_matrix(temp));
  temp[7, 1] = temp[7, 1] + vac_num[1, 1];                                 
  temp[7, 2] = temp[7, 2] + vac_num[1, 2] + vac_num[1, 3];
  temp[7, 3] = temp[7, 3] + vac_num[1, 4] + vac_num[1, 5];                                 
  temp[7, 4] = temp[7, 4] - vac_num[1, 6];
  
  for(n in 2:n_weeks){
    theta[5] = n;
    temp = integrate_ode_rk45(sir, temp[7], 0.0, ts, theta, x_r, x_i);
    y_out[n] = col_sums(to_matrix(temp));
    temp[7, 1] = temp[7, 1] + vac_num[n, 1];                                 
    temp[7, 2] = temp[7, 2] + vac_num[n, 2] + vac_num[n, 3];
    temp[7, 3] = temp[7, 3] + vac_num[n, 4] + vac_num[n, 5];                                 
    temp[7, 4] = temp[7, 4] - vac_num[n, 6];
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
  epsilon2 ~ normal(0.5, 0.2) T[epsilon1,1];
  epsilon3 ~ normal(0.6, 0.2) T[epsilon2,1];
  //sampling distribution
  for(i in 4:n_weeks){
    for(j in 1:3){
      target += neg_binomial_2_lpmf(cases[i,j] | y_out[i,j], phi);
    //target += poisson_lpmf(cases[i] | y_out[i]);
    }
  } 
}

generated quantities {
  real R0 = beta / gamma;
  real pred_cases[n_weeks, 3];
  for(i in 1:n_weeks){
    for(j in 1:3){
      pred_cases[i,j] = neg_binomial_2_rng(y_out[i,j], phi);
    //pred_cases[i] = poisson_rng(y_out[i]);
    }
  } 
}

