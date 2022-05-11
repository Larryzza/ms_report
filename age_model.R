rm(list=ls())
## Load packages
library(rstan) 
library(shinystan)
library(tidyverse)
library(truncnorm)
options(mc.cores=parallel::detectCores(logical = F))
rstan_options(auto_write = TRUE)
rstan_options(javascript=FALSE)
## set variable we need and define functions
need_var <- c("first_week_day", "Gender", "age_group", 
              "weekly_cases", "weekly_first_dose",
              "weekly_second_dose", "weekly_third_dose", "weekly_fourth_dose")
start_point <- "2021-11-28"
end_point <- "2022-03-05"
.sort_df <- function(x){
  x <- str_replace_all(x, "<", "")
  x <- ifelse(is.na(x)==T,0,x)
  return(as.numeric(x))
}

## Load data
df_case_vac <- read_csv("cases-among-vaccinated-211 (2).csv",
                        col_types = cols(.default = "c"))
df_case_vac[,c(3:17)] <- apply(as.matrix(df_case_vac[,c(3:17)]), 2, .sort_df)
df_case_vac$pos_2_dose <- rowSums(df_case_vac[,c(3:11)])
df_case_vac$pos_3and4_dose <- rowSums(df_case_vac[,c(12:16)])
df_case_vac$age_group <- "0-59"
temp <- grep("60|70|80|90",df_case_vac$Age_group)
df_case_vac$age_group[temp] <- "60+"
df_case_vac$date <- substr(df_case_vac$Week, 1, 10)

df_dose_case <- df_case_vac %>%
  filter(age_group!="NULL", date<=end_point, date>=start_point) %>% 
  group_by(date, age_group) %>% 
  summarise(weekly_2_dose_pos=sum(pos_2_dose),
            weekly_3and4_dose_pos=sum(pos_3and4_dose),
            weekly_0_dose_pos=sum(Sum_positive_without_vaccination))

case_group1_0 <- df_dose_case$weekly_0_dose_pos[which(df_dose_case$age_group=="0-59")]
case_group1_2 <- df_dose_case$weekly_2_dose_pos[which(df_dose_case$age_group=="0-59")]
case_group1_3and4 <- df_dose_case$weekly_3and4_dose_pos[which(df_dose_case$age_group=="0-59")]

case_group2_0 <- df_dose_case$weekly_0_dose_pos[which(df_dose_case$age_group!="0-59")]
case_group2_2 <- df_dose_case$weekly_2_dose_pos[which(df_dose_case$age_group!="0-59")]
case_group2_3and4 <- df_dose_case$weekly_3and4_dose_pos[which(df_dose_case$age_group!="0-59")]

# model v3
case_grouped <- cbind(case_group1_0+case_group2_0, case_group1_2+case_group2_2, 
                      case_group1_3and4+case_group2_3and4)

# model v4
#case_grouped <- cbind(case_group1_0+case_group2_0, case_group1_2+case_group2_2, 
#                      case_group1_3and4, case_group2_3and4)

df_daily <- read_csv("daily cases.csv") %>% 
  mutate(date=as.Date(date, "%d-%m-%Y")) %>% 
  filter(date<=end_point)
df_age_gender_vac <- read_csv("corona_age_and_gender_ver_00278.csv",
                              col_types = cols(.default = "c")) %>% 
                              select(all_of(need_var))
df_age_gender_vac[,c(4:8)] <- apply(as.matrix(df_age_gender_vac[,c(4:8)]),
                                    2, .sort_df)
df_age_gender_vac$Age_group <- "0-59"
temp <- grep("60|65|70|75|80",df_age_gender_vac$age_group)
df_age_gender_vac$Age_group[temp] <- "60+"
df_age_gender_vac$date <- df_age_gender_vac$first_week_day %>% as.Date

df <- df_age_gender_vac %>%
  filter(Gender!="unknown", age_group!="NULL", date<=end_point) %>% 
  group_by(date, Age_group) %>% 
  summarise(weekly_cases=sum(weekly_cases),
            weekly_first_dose=sum(weekly_first_dose),
            weekly_second_dose=sum(weekly_second_dose),
            weekly_third_dose=sum(weekly_third_dose),
            weekly_fourth_dose=sum(weekly_fourth_dose))

vac11 <- df$weekly_first_dose[df$date<start_point&df$Age_group=="0-59"] %>% sum
vac12 <- df$weekly_first_dose[df$date<start_point&df$Age_group!="0-59"] %>% sum

vac21 <- df$weekly_second_dose[df$date<start_point&df$Age_group=="0-59"] %>% sum
vac22 <- df$weekly_second_dose[df$date<start_point&df$Age_group!="0-59"] %>% sum

vac31 <- df$weekly_third_dose[df$date<start_point&df$Age_group=="0-59"] %>% sum
vac32 <- df$weekly_third_dose[df$date<start_point&df$Age_group!="0-59"] %>% sum

df <- df %>% filter(date>="2021-11-28")
cases_group1 <- df$weekly_cases[df$Age_group=="0-59"]
cases_group2 <- df$weekly_cases[df$Age_group!="0-59"]
cases <- df_daily$case
vac_num <- cbind(df$weekly_second_dose[df$Age_group=="0-59"]+
                   df$weekly_second_dose[df$Age_group!="0-59"],
                 df$weekly_third_dose[df$Age_group=="0-59"],
                 df$weekly_third_dose[df$Age_group!="0-59"],
                 df$weekly_fourth_dose[df$Age_group=="0-59"],
                 df$weekly_fourth_dose[df$Age_group!="0-59"],
                 df$weekly_second_dose[df$Age_group=="0-59"]+
                   df$weekly_second_dose[df$Age_group!="0-59"])
vac_num[,1] <- vac_num[,1]-vac_num[,2]-vac_num[,3]
vac_num[,2] <- vac_num[,2]-vac_num[,4]
vac_num[,3] <- vac_num[,3]-vac_num[,5]
# total population
N <- 9217000; I1 <- 1; E1 <- 5

# times
n_days <- length(cases) 
n_weeks <- length(cases_group1)
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions

## v2
'y0 <- c(V=vac21+vac22-vac31-vac32, 
        V31 = vac31, V31 = vac32, V41 = 0, V42 = 0, 
        S = N-I1-E1-vac21-vac22, E = E1, EV=0, EV1=0, EV2=0, 
        I=I1, IV=0, IV1=0, IV2=0, R = 0)'

## v3
y0 <- c(V2=vac21+vac22-vac31-vac32, 
        V3 = vac31 + vac32, V4 = 0, 
        S = N-I1*3-E1*3-vac21-vac22, E = E1, EV2=E1, EV3=E1, 
        I=I1, IV2=I1, IV3=I1)

## v4
'y0 <- c(V2=vac21+vac22-vac31-vac32, 
        V31 = vac31, V41 = 0, V32 = vac32, V42 = 0, 
        S = N-I1*4-E1*4-vac21-vac22, E = E1, EV2=E1, EV31=E1, EV32=E1,
        I=I1, IV2=I1, IV31=I1, IV32=I1)'

# data for Stan
data_sir <- list(n_weeks = n_weeks, y0 = y0, 
                 N = N, ts = 1:7, 
                 vac_num = vac_num, cases = case_grouped)


# number of MCMC steps
model_version <- "model_v3"
niter <- 10000
model <- stan_model(paste0(model_version, ".stan"))

## v3
initf <- function(chain_id = 1) {
  #set.seed(chain_id)
  list(beta = 1.5, 
       phi_inv = 1/10,
       day_shift = 1,
       eta = 0.1, 
       epsilon1 = 0.1,
       epsilon2 = 0.3,
       epsilon3 = 0.5,
       d = 0.7,
       k = 1)
} 


## v4
'initf <- function(chain_id = 1) {
  #set.seed(chain_id)
  list(beta = 1.5, 
       phi_inv = 1/10,
       day_shift = 6,
       eta = 0.1, 
       epsilon1 = 0.1,
       epsilon2 = 0.3,
       epsilon3 = 0.5,
       epsilon4 = 0.3,
       epsilon5 = 0.5,
       d = 0.7,
       k = 1)
} '

# generate a list of lists to specify initial values
n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

fit_seir <- sampling(model, init = init_ll,
                     data = data_sir,
                     iter = niter,
                     chains = n_chains, 
                     seed = 1996, 
                     control = list(adapt_delta=0.95, max_treedepth=15))

saveRDS(fit_seir, paste0(model_version, "_fit_seir.rds"))

## v3
pars=c("eta", "epsilon1","epsilon2","epsilon3",
       "R0","beta","k", "phi_inv","day_shift", "d")
## v4
'pars=c("eta", "epsilon1","epsilon2","epsilon3","epsilon4","epsilon5",
       "R0","beta","k", "phi_inv","day_shift", "d")'
#pars=c("eta","R0","beta","k")
print(fit_seir, pars = pars)
stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
traceplot(fit_seir, pars = pars)

num_of_group <- dim(case_grouped)[2]

case <- data.frame(case_grouped,1:14)
names(case) <- c(1:num_of_group,"t")
case <- pivot_longer(case,cols = 1:num_of_group)

seir_pred <- cbind(as.data.frame(summary(
  fit_seir, pars = "pred_cases", 
  probs = c(0.05, 0.5, 0.95))$summary), 
  t=rep(1:14,each=num_of_group),
  name=rep(1:num_of_group,14) %>% as.character) %>% 
  left_join(case, by = c("name","t"))
colnames(seir_pred) <- make.names(colnames(seir_pred)) # to remove % in the col names

ggplot(seir_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "blue", alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.) ) + 
  geom_point(mapping = aes(y = value)) +
  facet_wrap(~name) +
  labs(x = "Day", y = "Number of new cases") +
  theme_bw()

#launch_shinystan(fit_seir)
divergent <- get_sampler_params(fit_seir, inc_warmup=FALSE)[[1]]
table(divergent[,5])
#pairs(fit_seir)
