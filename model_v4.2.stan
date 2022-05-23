
functions {

  real switch_eta(real eta, real week_index, real k, real week_shift) {
    return(eta + (1 - eta) / (1 + exp(k * (week_index - 6 - week_shift))));
  }
  
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      
      real V2 = y[1];
      real V31 = y[2];
      real V41 = y[3];
      real V32 = y[4];
      real V42 = y[5];
      real S = y[6];
      real S2 = y[7];
      real S31 = y[8];
      real S32 = y[9];
      real E = y[10];
      real EV2 = y[11];
      real EV31 = y[12];
      real EV32 = y[13];
      real I = y[14];
      real IV2 = y[15];
      real IV31 = y[16];
      real IV32 = y[17];
      
      real N = x_i[1];
      real beta = theta[1];
      real gamma = theta[2];
      real sigma = theta[3];
      real eta = theta[4];
      real week_index = theta[5];
      real epsilon4 = theta[6];
      real epsilon5 = theta[7];
      real alpha1 = theta[8];
      real alpha2 = theta[9];
      real k = theta[10];
      real week_shift = theta[11];
      real d = theta[12];
      real d_old = theta[13];
      
      real c = switch_eta(eta, week_index, k, week_shift); // switch function with k
      
      real dV2_dt = - 0.0129 * V2 - 
                    (1 - 65.5 / 100) * beta * c * (I + IV2 + IV31 + IV32) * V2 / N;
      real dV31_dt = - 0.0078 * V31 - 
                     (1 - 67.2 / 100) * beta * c * (I + IV2 + IV31 + IV32) * V31 / N;
      real dV41_dt = - alpha1 * V41  - 
                     (1 - epsilon4) * beta * c * (I + IV2 + IV31 + IV32) * V41 / N;
      real dV32_dt = - 0.0078 * V32 - 
                     (1 - 67.2 / 100) * beta * c * (I + IV2 + IV31 + IV32) * V32 / N;
      real dV42_dt = - alpha2 * V42 - 
                     (1 - epsilon5) * beta * c * (I + IV2 + IV31 + IV32) * V42 / N;
      
      real dS_dt = - d * beta * c * (I + IV2 + IV31 + IV32) * S / N;
      real dS2_dt = 0.0129 * V2 - beta * c * (I + IV2 + IV31 + IV32) * S2 / N;
      real dS31_dt = 0.0078 * V31 + alpha1 * V41 - 
                     beta * c * (I + IV2 + IV31 + IV32) * S31 / N;
      real dS32_dt = 0.0078 * V32 + alpha2 * V42 - 
                     d_old * beta * c * (I + IV2 + IV31 + IV32) * S32 / N;
      
      real dE_dt = d * beta * c * (I + IV2 + IV31 + IV32) * S / N - sigma * E;
      real dEV2_dt = (1 - 65.5 / 100) * beta * c * (I + IV2 + IV31 + IV32) * V2 / N +
                     beta * c * (I + IV2 + IV31 + IV32) * S2 / N - 
                     sigma * EV2;
      real dEV31_dt = (1 - 67.2 / 100) * beta * c * (I + IV2 + IV31 + IV32) * V31 / N + 
                      (1 - epsilon4) * beta * c * (I + IV2 + IV31 + IV32) * V41 / N +
                      beta * c * (I + IV2 + IV31 + IV32) * S31 / N - 
                      sigma * EV31;
      real dEV32_dt = (1 - 67.2 / 100) * beta * c * (I + IV2 + IV31 + IV32) * V32 / N + 
                      (1 - epsilon5) * beta * c * (I + IV2 + IV31 + IV32) * V42 / N +
                      d_old * beta * c * (I + IV2 + IV31 + IV32) * S32 / N - 
                      sigma * EV32;
                     
      real dI_dt =  sigma * E - gamma * I;
      real dIV2_dt =  sigma * EV2 - gamma * IV2;
      real dIV31_dt =  sigma * EV31 - gamma * IV31;
      real dIV32_dt =  sigma * EV32 - gamma * IV32;
      
      return {dV2_dt, dV31_dt, dV41_dt, dV32_dt, dV42_dt, 
              dS_dt, dS2_dt, dS31_dt, dS32_dt,
              dE_dt, dEV2_dt, dEV31_dt, dEV32_dt, 
              dI_dt, dIV2_dt, dIV31_dt, dIV32_dt};
  }
  
  real[] col_sums(matrix X) {
     real s[4] ;
     int col_index[4] = {14, 15, 16, 17};
     for (j in 1:4){
       s[j] = sum(col(X, col_index[j]));
     } 
	   return s ;
  }
}

data {
  int<lower=1> n_weeks;
  real y0[17]; 
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
  //real<lower=0> beta;
  //real<lower=0> sigma;
  real<lower=0> phi_inv;
  //real week_shift;
  //real<lower=0, upper=1> eta;
  //real<lower=0, upper=1> epsilon1;
  //real<lower=0, upper=1> epsilon2;
  //real<lower=0, upper=1> epsilon3;
  real<lower=0, upper=1> epsilon4;
  real<lower=0, upper=1> epsilon5;
  real<lower=0> alpha1;
  real<lower=0> alpha2;
  real<lower=0> d_old;
  //real<lower=0> k;
}

transformed parameters{
  real y_out[n_weeks, 4];
  real temp[7,17];
  //real gamma = 1./5;
  real phi = 1./phi_inv;
  real theta[13] = {0.95, 1./5, 1./3, 0.14, 1.0,
                   epsilon4, epsilon5, alpha1 / 1000, alpha2 / 1000, 
                   1.30, 1.11, 0.59, d_old};

  temp = integrate_ode_rk45(sir, y0, 0.0, ts, theta, x_r, x_i);
  y_out[1] = col_sums(to_matrix(temp));
  temp[7, 1] = temp[7, 1] + vac_num[1, 1];                                 
  temp[7, 2] = temp[7, 2] + vac_num[1, 2];
  temp[7, 3] = temp[7, 3] + vac_num[1, 4];  
  temp[7, 4] = temp[7, 4] + vac_num[1, 3];
  temp[7, 5] = temp[7, 5] + vac_num[1, 5];
  temp[7, 6] = temp[7, 6] - vac_num[1, 6];
  //print(temp);
  
  for(n in 2:n_weeks){
    theta[5] = n;
    temp = integrate_ode_rk45(sir, temp[7], 0.0, ts, theta, x_r, x_i);
    y_out[n] = col_sums(to_matrix(temp));
    temp[7, 1] = temp[7, 1] + vac_num[n, 1];                                 
    temp[7, 2] = temp[7, 2] + vac_num[n, 2];
    temp[7, 3] = temp[7, 3] + vac_num[n, 4];  
    temp[7, 4] = temp[7, 4] + vac_num[n, 3];
    temp[7, 5] = temp[7, 5] + vac_num[n, 5];
    temp[7, 6] = temp[7, 6] - vac_num[n, 6];
    //print(temp);
  }
}

model {
  //priors
  //beta ~ normal(2, 1) T[0,];
  //k ~ normal(1, 2) T[0,];
  //week_shift ~ normal(6, 3);
  //gamma ~ normal(0.5, 0.3);
  //sigma ~ normal(0.3, 0.3) T[0,1];
  //d ~ normal(0.7, 0.2) T[0,1];
  phi_inv ~ exponential(10);
  d_old ~ normal(0.7, 0.5) T[0,];
  //eta ~ beta(2, 4);
  epsilon4 ~ normal(0.4, 0.1) T[0,1];
  epsilon5 ~ normal(0.7, 0.1) T[0,1];
  alpha1 ~ normal(1, 3) T[0,];
  alpha2 ~ normal(1, 3) T[0,];
  //sampling distribution
  for(i in 4:n_weeks){
    for(j in 1:4){
      target += neg_binomial_2_lpmf(cases[i,j] | y_out[i,j], phi);
    }
  } 
}

generated quantities {
  real protect_rate_10week_young = epsilon4 * (1-alpha1/1000)^70;
  real protect_rate_10week_old = epsilon5 * (1-alpha2/1000)^70;
  real pred_cases[n_weeks, 4];
  for(i in 1:n_weeks){
    for(j in 1:4){
      pred_cases[i,j] = neg_binomial_2_rng(y_out[i,j], phi);
    }
  } 
}

