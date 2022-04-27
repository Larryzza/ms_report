rm(list=ls())
## Load packages
library(rstan)
library(shinystan) 
library(tidyverse)
options(mc.cores=parallel::detectCores(logical = F))
rstan_options(auto_write = TRUE)
rstan_options(javascript=FALSE)
## Load data
need_var <- c("first_week_day", "Gender", "age_group", 
              "weekly_cases", "weekly_first_dose",
              "weekly_second_dose", "weekly_third_dose", "weekly_fourth_dose")
df_daily <- read_csv("daily cases.csv")
df_age_gender_vac <- read_csv("corona_age_and_gender_ver_00278.csv",
                              col_types = cols(.default = "c")) %>% 
  select(all_of(need_var))
df_age_gender_vac[,c(4:8)] <- apply(as.matrix(df_age_gender_vac[,c(4:8)]),
                                    2, function(x){
                                      x <- str_replace_all(x, "<", "")
                                      x <- ifelse(is.na(x)==T,0,x)
                                      return(as.numeric(x))})
df_age_gender_vac$Age_group <- "0-59"
temp <- grep("60|65|70|75|80",df_age_gender_vac$age_group)
df_age_gender_vac$Age_group[temp] <- "60+"
df_age_gender_vac$date <- df_age_gender_vac$first_week_day %>% as.Date

df <- df_age_gender_vac %>%
  filter(Gender!="unknown", age_group!="NULL") %>% 
  group_by(date, Age_group) %>% 
  summarise(weekly_cases=sum(weekly_cases),
            weekly_first_dose=sum(weekly_first_dose),
            weekly_second_dose=sum(weekly_second_dose),
            weekly_third_dose=sum(weekly_third_dose),
            weekly_fourth_dose=sum(weekly_fourth_dose))

vac11 <- df$weekly_first_dose[df$date<"2021-12-05"&df$Age_group=="0-59"] %>% sum
vac12 <- df$weekly_first_dose[df$date<"2021-12-05"&df$Age_group!="0-59"] %>% sum

vac21 <- df$weekly_second_dose[df$date<"2021-12-05"&df$Age_group=="0-59"] %>% sum
vac22 <- df$weekly_second_dose[df$date<"2021-12-05"&df$Age_group!="0-59"] %>% sum

vac31 <- df$weekly_third_dose[df$date<"2021-12-05"&df$Age_group=="0-59"] %>% sum
vac32 <- df$weekly_third_dose[df$date<"2021-12-05"&df$Age_group!="0-59"] %>% sum

df <- df %>% filter(date>="2021-12-05")
cases_group1 <- df$weekly_cases[df$Age_group=="0-59"]
cases_group2 <- df$weekly_cases[df$Age_group!="0-59"]
cases <- df_daily$case
vac_num <- cbind(df$weekly_second_dose[df$Age_group=="0-59"]+
                   df$weekly_second_dose[df$Age_group!="0-59"],
                 df$weekly_third_dose[df$Age_group=="0-59"]+
                   df$weekly_third_dose[df$Age_group!="0-59"],
                 df$weekly_fourth_dose[df$Age_group=="0-59"],
                 df$weekly_fourth_dose[df$Age_group!="0-59"],
                 df$weekly_second_dose[df$Age_group=="0-59"]+
                   df$weekly_second_dose[df$Age_group!="0-59"])
vac_num[,1] <- vac_num[,1]-vac_num[,2]
vac_num[,2] <- vac_num[,2]-vac_num[,3]-vac_num[,4]
# total population
N <- 9217000; E1 <- 100

# times
n_days <- length(cases) 
n_weeks <- length(cases_group1)
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions

y0 <- c(V=vac21+vac22-vac31-vac32, 
        V3 = vac31+vac32, V41 = 0, V42 = 0, 
        S = N-E1-vac21-vac22, E = E1, I=0, R = 0)

# data for Stan
data_sir <- list(n_days = n_days, n_weeks = n_weeks, y0 = y0, 
                 N = N, ts = 1:7, vac_num = vac_num,
                 #t0 = t0,   
                 cases = cases)


# number of MCMC steps
niter <- 2000
model <- stan_model("model.stan")
fit_seir <- sampling(model,
                     data = data_sir,
                     iter = niter,
                     chains = 3, 
                     seed = 1)
saveRDS(fit_seir, "fit_seir.rds")
pars=c('beta', 'gamma', "R0","a")
print(fit_seir, pars = pars)
stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
traceplot(fit_seir, pars = c("gamma", "beta", "R0"))


smr_pred <- cbind(as.data.frame(summary(
  fit_seir, pars = "pred_cases", 
  probs = c(0.05, 0.5, 0.95))$summary), t,
  cases = cases)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "blue", alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.) ) + 
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Number of students in bed")
