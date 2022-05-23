
rm(list=ls())
library(tidyverse)

df_age_gender_vac <- read_csv("corona_age_and_gender_ver_00278.csv",
                              col_types = cols(.default = "c"))
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
  filter(Gender!="unknown", age_group!="NULL", 
         date>= start_point, date<=end_point) %>%
  group_by(date, Age_group) %>% 
  summarise(weekly_cases=sum(weekly_cases),
            weekly_first_dose=sum(weekly_first_dose),
            weekly_second_dose=sum(weekly_second_dose),
            weekly_third_dose=sum(weekly_third_dose),
            weekly_fourth_dose=sum(weekly_fourth_dose))
df$pop <- ifelse(df$Age_group=="0-59", 7793000, 1424000) 
df$crude_rate <- df$weekly_cases/df$pop
df$date <- as.character(df$date)
ggplot(df, aes(x=date, y=crude_rate, group=Age_group)) +
  geom_line(aes(color=Age_group))+
  #geom_point(aes(color=Age_group))+
  #scale_y_sqrt()+
  scale_x_discrete(breaks=c("2021-11-28", "2022-01-16", "2022-02-27"))+
  theme_bw()+
  ylab(label="new cases / population within age group")+
  xlab(label=NULL)+
  scale_color_discrete(name = "Age group")+
  #labs(title="crude covid rate plot")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
#ggsave("crude_rate.png")
ggsave("crude_rate.png",width = 10,height = 5)
