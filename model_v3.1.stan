
functions {

  real switch_eta(real eta, real week_index, real k, real week_shift) {
    return(eta + (1 - eta) / (1 + exp(k * (week_index - 6 - week_shift))));
  }
  
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      
      real V2 = y[1];
      real V3 = y[2];
      real V4 = y[3];
      real S = y[4];
      real S2 = y[5];
      real S3 = y[6];
      real E = y[7];
      real EV2 = y[8];
      real EV3 = y[9];
      real I = y[10];
      real IV2 = y[11];
      real IV3 = y[12];
      real R = y[13];
      real RV2 = y[14];
      real RV3 = y[15];
      
      real N = x_i[1];
      real beta = theta[1];
      real gamma = theta[2];
      real sigma = theta[3];
      real eta = theta[4];
      real week_index = theta[5];
      real epsilon1 = theta[6];
      real epsilon2 = theta[7];
      real epsilon3 = theta[8];
      real alpha1 = theta[9];
      real k = theta[10];
      real week_shift = theta[11];
      real d = theta[12];
      
      real c = switch_eta(eta, week_index, k, week_shift); // switch function with k
      
      real dV2_dt = - 0.0129 * V2 - 
                    (1 - epsilon1) * beta * c * (I + IV2 + IV3) * V2 / N;
      real dV3_dt = - 0.0078 * V3 - 
                    (1 - epsilon2) * beta * c * (I + IV2 + IV3) * V3 / N;
      real dV4_dt = - alpha1 * V4 - 
                    (1 - epsilon3) * beta * c * (I + IV2 + IV3) * V4 / N;
      
      real dS_dt = 0.0001 * R - d * beta * c * (I + IV2 + IV3) * S / N;
      real dS2_dt = 0.0001 * RV2 + 0.0129 * V2 - 
                    beta * c * (I + IV2 + IV3) * S2 / N;
      real dS3_dt = 0.0001 * RV3 + 0.0078 * V3 + alpha1 * V4 - 
                    beta * c * (I + IV2 + IV3) * S3 / N;
      
      real dE_dt = d * beta * c * (I + IV2 + IV3) * S / N - sigma * E;
      real dEV2_dt = (1 - epsilon1) * beta * c * (I + IV2 + IV3) * V2 / N + 
                     beta * c * (I + IV2 + IV3) * S2 / N -
                     sigma * EV2;
      real dEV3_dt = (1 - epsilon2) * beta * c * (I + IV2 + IV3) * V3 / N + 
                     (1 - epsilon3) * beta * c * (I + IV2 + IV3) * V4 / N +
                     beta * c * (I + IV2 + IV3) * S3 / N - 
                     sigma * EV3;
                     
      real dI_dt =  sigma * E - gamma * I;
      real dIV2_dt =  sigma * EV2 - gamma * IV2;
      real dIV3_dt =  sigma * EV3 - gamma * IV3;
      
      real dR_dt = gamma * I - 0.0001 * R;
      real dRV2_dt = gamma * IV2 - 0.0001 * RV2;
      real dRV3_dt = gamma * IV3 - 0.0001 * RV3;
      
      return {dV2_dt, dV3_dt, dV4_dt, 
              dS_dt, dS2_dt, dS3_dt,
              dE_dt, dEV2_dt, dEV3_dt,  
              dI_dt, dIV2_dt, dIV3_dt,
              dR_dt, dRV2_dt, dRV3_dt};
  }
  
  real[] col_sums(matrix X) {
     real s[3] ;
     int col_index[3] = {10, 11, 12};
     for (j in 1:3) s[j] = sum(col(X, col_index[j])) ;
	   return s ;
  }
}

data {
  int<lower=1> n_weeks;
  real y0[15];
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
  real<lower=0> phi_inv;
  real week_shift;
  real<lower=0, upper=1> eta;
  //real<lower=0, upper=1> epsilon1;
  //real<lower=0, upper=1> epsilon2;
  real<lower=0, upper=1> epsilon3;
  real<lower=0> alpha;
  real<lower=0, upper=1> d;
  //real<lower=0, upper=1> epsilon4;
  //real<lower=0, upper=1> epsilon5;
  real<lower=0> k;
}

transformed parameters{
  real y_out[n_weeks, 3];
  real temp[7,15];
  real gamma = 1./5;
  real phi = 1./phi_inv;
  real theta[12] = {beta, gamma, 1./3, eta, 1.0,
                   65.5 / 100, 67.2 / 100, epsilon3, alpha / 1000,
                   k, week_shift , d};

  temp = integrate_ode_rk45(sir, y0, 0.0, ts, theta, x_r, x_i);
  y_out[1] = col_sums(to_matrix(temp));
  temp[7, 1] = temp[7, 1] + vac_num[1, 1];                                 
  temp[7, 2] = temp[7, 2] + vac_num[1, 2] + vac_num[1, 3];
  temp[7, 3] = temp[7, 3] + vac_num[1, 4] + vac_num[1, 5];                                 
  temp[7, 4] = temp[7, 4] - vac_num[1, 6];
  //print(temp);
  
  for(n in 2:n_weeks){
    theta[5] = n;
    temp = integrate_ode_rk45(sir, temp[7], 0.0, ts, theta, x_r, x_i);
    y_out[n] = col_sums(to_matrix(temp));
    temp[7, 1] = temp[7, 1] + vac_num[n, 1];                                 
    temp[7, 2] = temp[7, 2] + vac_num[n, 2] + vac_num[n, 3];
    temp[7, 3] = temp[7, 3] + vac_num[n, 4] + vac_num[n, 5];                                 
    temp[7, 4] = temp[7, 4] - vac_num[n, 6];
    //print(temp);
  }
}

model {
  //priors
  beta ~ normal(2, 1) T[0,];
  k ~ normal(1, 2) T[0,];
  week_shift ~ normal(0, 2);
  //gamma ~ normal(0.5, 0.3);
  //sigma ~ normal(0.3, 0.3) T[0,1];
  alpha ~ normal(1, 3) T[0,];
  d ~ normal(0.7, 0.2) T[0,1];
  phi_inv ~ exponential(10);
  eta ~ beta(2, 4);
  //epsilon1 ~ normal(0.7, 0.3) T[0,1];
  //epsilon2 ~ normal(0.7, 0.3) T[0,1];
  epsilon3 ~ normal(0.7, 0.1) T[0,1];
  //sampling distribution
  for(i in 4:n_weeks){
    for(j in 1:3){
      target += neg_binomial_2_lpmf(cases[i,j] | y_out[i,j], phi);
    }
  } 
}

generated quantities {
  real R0 = beta / gamma;
  real pred_cases[n_weeks, 3];
  real protect_rate_10week = epsilon3 * (1-alpha/1000)^70;
  for(i in 1:n_weeks){
    for(j in 1:3){
      pred_cases[i,j] = neg_binomial_2_rng(y_out[i,j], phi);
    }
  } 
}

