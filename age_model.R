rm(list=ls())
## Load packages
library(rstan) 
library(tidyverse)
options(mc.cores=parallel::detectCores(logical = F))
rstan_options(auto_write = TRUE)
rstan_options(javascript=FALSE)
## Load data
need_var <- c("first_week_day", "Gender", "age_group", 
              "weekly_cases", "weekly_first_dose",
              "weekly_second_dose", "weekly_third_dose", "weekly_fourth_dose")
start_point <- "2021-11-28"
end_point <- "2022-03-05"

df_daily <- read_csv("daily cases.csv") %>% 
  mutate(date=as.Date(date, "%d-%m-%Y")) %>% 
  filter(date<=end_point)
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
                 df$weekly_third_dose[df$Age_group=="0-59"]+
                   df$weekly_third_dose[df$Age_group!="0-59"],
                 df$weekly_fourth_dose[df$Age_group=="0-59"],
                 df$weekly_fourth_dose[df$Age_group!="0-59"],
                 df$weekly_second_dose[df$Age_group=="0-59"]+
                   df$weekly_second_dose[df$Age_group!="0-59"])
vac_num[,1] <- vac_num[,1]-vac_num[,2]
vac_num[,2] <- vac_num[,2]-vac_num[,3]-vac_num[,4]
# total population
N <- 9217000; I1 <- 1; E1 <- 10

# times
n_days <- length(cases) 
n_weeks <- length(cases_group1)
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions

y0 <- c(V=vac21+vac22-vac31-vac32, 
        V3 = vac31+vac32, V41 = 0, V42 = 0, 
        S = N-I1-E1-vac21-vac22, E = E1, I=I1, R = 0)

# data for Stan
data_sir <- list(n_days = n_days, n_weeks = n_weeks, y0 = y0, 
                 N = N, ts = 1:7, vac_num = vac_num, cases = cases)


# number of MCMC steps
niter <- 100000
model <- stan_model("model.stan")

initf <- function(chain_id = 1) {
  set.seed(chain_id^2)
  list(beta = rnorm(n = 1, mean = 2, sd = 1), 
       eta = rbeta(n = 1, shape1 = 2, shape2 = 4), 
       epsilon1 = 0.3,
       epsilon2 = 0.6,
       epsilon3 = 0.7,
       epsilon4 = 0.7)
} 

# generate a list of lists to specify initial values
n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

fit_seir <- sampling(model, init = init_ll,
                     data = data_sir,
                     iter = niter,
                     chains = 4, 
                     seed = 1)

saveRDS(fit_seir, "fit_seir.rds")

pars=c("eta", "epsilon1","epsilon2","epsilon3","epsilon4","R0","beta")
print(fit_seir, pars = pars)
stan_dens(fit_seir, pars = pars, separate_chains = TRUE)
traceplot(fit_seir, pars = pars)

smr_pred <- cbind(as.data.frame(summary(
  fit_seir, pars = "pred_cases", 
  probs = c(0.05, 0.5, 0.95))$summary), t,
  cases = cases)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "blue", alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.) ) + 
  geom_point(mapping = aes(y = cases)) +
  #scale_y_sqrt()+
  labs(x = "Day", y = "Number of new cases")

