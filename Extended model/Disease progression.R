# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)

# Define useful paths
path_paper <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/" # Where results from other studies are saved
path_model <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/Vaccination/" #Where model outputs are stored
  
# Read in RSV-illness rates
mild_illness_rate <- read_excel(paste0(path_paper,"RSV illness rates SA.xlsx"), sheet = "Mild illness")
severe_illness_rate <- read_excel(paste0(path_paper,"RSV illness rates SA.xlsx"), sheet = "Severe illness")

# Read in model estimate
conversion_rate_summer <- read.csv(paste0(path_model, "Incidence/strat_summer_conv_inc.csv")) %>% 
  subset(select = -X)
conversion_rate_autumn <- read.csv(paste0(path_model, "Incidence/strat_autumn_conv_inc.csv")) %>% 
  subset(select = -X)
conversion_rate_winter <- read.csv(paste0(path_model, "Incidence/strat_winter_conv_inc.csv")) %>% 
  subset(select = -X)
conversion_rate_spring <- read.csv(paste0(path_model, "Incidence/strat_spring_conv_inc.csv")) %>% 
  subset(select = -X)

# Combine all estimates into one DF
conversion <- conversion_rate_summer %>% mutate(birth_season = "summer") %>%
  rbind(
    conversion_rate_autumn %>% mutate(birth_season = "autumn")
  ) %>%
  rbind(
    conversion_rate_winter %>% mutate(birth_season = "winter")
  ) %>%
  rbind(
    conversion_rate_spring %>% mutate(birth_season = "spring")
  )

# Convert ages to the same units
conversion_formated <- conversion %>% mutate(age_months = trunc(agemid/30.4375)) %>% # turn age into months
  arrange(agemid) # arrange by age instead of by season of birth

# Put ages in the same categories as the disease progression data
conversion_formated <- conversion_formated %>% 
  mutate(age_bracket = case_when(age_months < 1 ~ "<1",
                                 age_months > 11 & age_months < 15 ~ "12-14",
                                 age_months > 14 & age_months < 18 ~ "15-17",
                                 age_months > 17 & age_months < 21 ~ "18-20",
                                 age_months > 20 & age_months < 24 ~ "21-23",
                                 age_months > 23 & age_months < 36 ~ "24-35",
                                 age_months > 35 & age_months < 48 ~ "36-47",
                                 age_months > 48 & age_months < 60 ~ "48-59",
                                 age_months > 59 ~ "≥60",
                                 TRUE ~ as.character(age_months)))

# Define which season of the year the birth cohorts are in
conversion_formated <- conversion_formated %>% 
  group_by(agemid, birth_season) %>%
  mutate(current_season = case_when(birth_season == "summer" & ((agemid <= 30.41*3)          # season in which the cohort was born                    
                                    || ((agemid >= 365) & (agemid <= (365+30.41*3))) 
                                    || ((agemid >= 2*365) & (agemid <= (2*365+30.41*3))) 
                                    || ((agemid >= 3*365) & (agemid <= (3*365+30.41*3))) 
                                    || ((agemid >= 4*365) & (agemid <= (4*365+30.41*3))) 
                                    || ((agemid >= 5*365) & (agemid <= (5*365+30.41*3)))) ~ "summer",
                                    birth_season == "autumn" & ((agemid <= 30.41*3)                              
                                    || ((agemid >= 365) & (agemid <= (365+30.41*3))) 
                                    || ((agemid >= 2*365) & (agemid <= (2*365+30.41*3))) 
                                    || ((agemid >= 3*365) & (agemid <= (3*365+30.41*3))) 
                                    || ((agemid >= 4*365) & (agemid <= (4*365+30.41*3))) 
                                    || ((agemid >= 5*365) & (agemid <= (5*365+30.41*3)))) ~ "autumn",
                                    birth_season == "winter" & ((agemid <= 30.41*3)                              
                                    || ((agemid >= 365) & (agemid <= (365+30.41*3))) 
                                    || ((agemid >= 2*365) & (agemid <= (2*365+30.41*3))) 
                                    || ((agemid >= 3*365) & (agemid <= (3*365+30.41*3))) 
                                    || ((agemid >= 4*365) & (agemid <= (4*365+30.41*3))) 
                                    || ((agemid >= 5*365) & (agemid <= (5*365+30.41*3)))) ~ "winter",
                                    birth_season == "spring" & ((agemid <= 30.41*3)                              
                                    || ((agemid >= 365) & (agemid <= (365+30.41*3))) 
                                    || ((agemid >= 2*365) & (agemid <= (2*365+30.41*3))) 
                                    || ((agemid >= 3*365) & (agemid <= (3*365+30.41*3))) 
                                    || ((agemid >= 4*365) & (agemid <= (4*365+30.41*3))) 
                                    || ((agemid >= 5*365) & (agemid <= (5*365+30.41*3)))) ~ "spring",
                                    
                                    birth_season == "summer" & ((agemid > 30.41*3 & agemid <= 30.41*6 ) # season following the one the cohort was born in
                                    || ((agemid > 365 + 30.41*3) & (agemid <= (365+30.41*6))) 
                                    || ((agemid > 2*365 + 30.41*3) & (agemid <= (2*365+30.41*6))) 
                                    || ((agemid > 3*365 + 30.41*3) & (agemid <= (3*365+30.41*6))) 
                                    || ((agemid > 4*365 + 30.41*3) & (agemid <= (4*365+30.41*6))) 
                                    || ((agemid > 5*365 + 30.41*3) & (agemid <= (5*365+30.41*6)))) ~ "autumn",
                                    birth_season == "autumn" & ((agemid > 30.41*3 & agemid <= 30.41*6 ) 
                                    || ((agemid > 365 + 30.41*3) & (agemid <= (365+30.41*6))) 
                                    || ((agemid > 2*365 + 30.41*3) & (agemid <= (2*365+30.41*6))) 
                                    || ((agemid > 3*365 + 30.41*3) & (agemid <= (3*365+30.41*6))) 
                                    || ((agemid > 4*365 + 30.41*3) & (agemid <= (4*365+30.41*6))) 
                                    || ((agemid > 5*365 + 30.41*3) & (agemid <= (5*365+30.41*6)))) ~ "winter",
                                    birth_season == "winter" & ((agemid > 30.41*3 & agemid <= 30.41*6 ) 
                                    || ((agemid > 365 + 30.41*3) & (agemid <= (365+30.41*6))) 
                                    || ((agemid > 2*365 + 30.41*3) & (agemid <= (2*365+30.41*6))) 
                                    || ((agemid > 3*365 + 30.41*3) & (agemid <= (3*365+30.41*6))) 
                                    || ((agemid > 4*365 + 30.41*3) & (agemid <= (4*365+30.41*6))) 
                                    || ((agemid > 5*365 + 30.41*3) & (agemid <= (5*365+30.41*6)))) ~ "spring",
                                    birth_season == "spring" & ((agemid > 30.41*3 & agemid <= 30.41*6 ) 
                                    || ((agemid > 365 + 30.41*3) & (agemid <= (365+30.41*6))) 
                                    || ((agemid > 2*365 + 30.41*3) & (agemid <= (2*365+30.41*6))) 
                                    || ((agemid > 3*365 + 30.41*3) & (agemid <= (3*365+30.41*6))) 
                                    || ((agemid > 4*365 + 30.41*3) & (agemid <= (4*365+30.41*6))) 
                                    || ((agemid > 5*365 + 30.41*3) & (agemid <= (5*365+30.41*6)))) ~ "summer",
                                    
                                    birth_season == "summer" & ((agemid > 30.41*6 & agemid <= 30.41*9 )  # two seasons after birth
                                    || ((agemid > 365+30.41*6) & (agemid <= (365+30.41*9))) 
                                    || ((agemid > 2*365 + 30.41*6) & (agemid <= (2*365+30.41*9))) 
                                    || ((agemid > 3*365 + 30.41*6) & (agemid <= (3*365+30.41*9))) 
                                    || ((agemid > 4*365 + 30.41*6) & (agemid <= (4*365+30.41*9))) 
                                    || ((agemid > 5*365 + 30.41*6) & (agemid <= (5*365+30.41*9)))) ~ "winter",
                                    birth_season == "autumn" & ((agemid > 30.41*6 & agemid <= 30.41*9 )  
                                    || ((agemid > 365+30.41*6) & (agemid <= (365+30.41*9))) 
                                    || ((agemid > 2*365 + 30.41*6) & (agemid <= (2*365+30.41*9))) 
                                    || ((agemid > 3*365 + 30.41*6) & (agemid <= (3*365+30.41*9))) 
                                    || ((agemid > 4*365 + 30.41*6) & (agemid <= (4*365+30.41*9))) 
                                    || ((agemid > 5*365 + 30.41*6) & (agemid <= (5*365+30.41*9)))) ~ "spring",
                                    birth_season == "winter" & ((agemid > 30.41*6 & agemid <= 30.41*9 )  
                                    || ((agemid > 365+30.41*6) & (agemid <= (365+30.41*9))) 
                                    || ((agemid > 2*365 + 30.41*6) & (agemid <= (2*365+30.41*9))) 
                                    || ((agemid > 3*365 + 30.41*6) & (agemid <= (3*365+30.41*9))) 
                                    || ((agemid > 4*365 + 30.41*6) & (agemid <= (4*365+30.41*9))) 
                                    || ((agemid > 5*365 + 30.41*6) & (agemid <= (5*365+30.41*9)))) ~ "summer",
                                    birth_season == "spring" & ((agemid > 30.41*6 & agemid <= 30.41*9 )  
                                    || ((agemid > 365+30.41*6) & (agemid <= (365+30.41*9))) 
                                    || ((agemid > 2*365 + 30.41*6) & (agemid <= (2*365+30.41*9))) 
                                    || ((agemid > 3*365 + 30.41*6) & (agemid <= (3*365+30.41*9))) 
                                    || ((agemid > 4*365 + 30.41*6) & (agemid <= (4*365+30.41*9))) 
                                    || ((agemid > 5*365 + 30.41*6) & (agemid <= (5*365+30.41*9)))) ~ "autumn",
                                    
                                    birth_season == "summer" & ((agemid > 30.41*9 & agemid <= 30.41*12 ) # 3 seasons after birth                    
                                    || ((agemid > 365+30.41*9) & (agemid <= (365+30.41*12))) 
                                    || ((agemid > 2*365 + 30.41*9) & (agemid <= (2*365+30.41*12))) 
                                    || ((agemid > 3*365 + 30.41*9) & (agemid <= (3*365+30.41*12))) 
                                    || ((agemid > 4*365 + 30.41*9) & (agemid <= (4*365+30.41*12))) 
                                    || ((agemid > 5*365 + 30.41*9) & (agemid <= (5*365+30.41*12)))) ~ "spring",
                                    birth_season == "autumn" & ((agemid > 30.41*9 & agemid <= 30.41*12 )                     
                                    || ((agemid > 365+30.41*9) & (agemid <= (365+30.41*12))) 
                                    || ((agemid > 2*365 + 30.41*9) & (agemid <= (2*365+30.41*12))) 
                                    || ((agemid > 3*365 + 30.41*9) & (agemid <= (3*365+30.41*12))) 
                                    || ((agemid > 4*365 + 30.41*9) & (agemid <= (4*365+30.41*12))) 
                                    || ((agemid > 5*365 + 30.41*9) & (agemid <= (5*365+30.41*12)))) ~ "summer",
                                    birth_season == "winter" & ((agemid > 30.41*9 & agemid <= 30.41*12 )                     
                                    || ((agemid > 365+30.41*9) & (agemid <= (365+30.41*12))) 
                                    || ((agemid > 2*365 + 30.41*9) & (agemid <= (2*365+30.41*12))) 
                                    || ((agemid > 3*365 + 30.41*9) & (agemid <= (3*365+30.41*12))) 
                                    || ((agemid > 4*365 + 30.41*9) & (agemid <= (4*365+30.41*12))) 
                                    || ((agemid > 5*365 + 30.41*9) & (agemid <= (5*365+30.41*12)))) ~ "autumn",
                                    birth_season == "spring" & ((agemid > 30.41*9 & agemid <= 30.41*12 )                   
                                    || ((agemid > 365+30.41*9) & (agemid <= (365+30.41*12))) 
                                    || ((agemid > 2*365 + 30.41*9) & (agemid <= (2*365+30.41*12))) 
                                    || ((agemid > 3*365 + 30.41*9) & (agemid <= (3*365+30.41*12))) 
                                    || ((agemid > 4*365 + 30.41*9) & (agemid <= (4*365+30.41*12))) 
                                    || ((agemid > 5*365 + 30.41*9) & (agemid <= (5*365+30.41*12)))) ~ "winter"
                                    
  ))

# Convert illness rates to proportions (prop = 1-exp(1/rate*T))
mild_illness_prop <- mild_illness_rate %>%
  mutate(total_MI_prop = 1-exp(-`Total mild illness rate`/100000),
         total_MI_prop_low = 1-exp(-`Total mild illness lower CI`/100000),
         total_MI_prop_up = 1-exp(-`Total mild illness upper CI`/100000),
         total_MA_MI_prop = 1-exp(-`MA mild illness rate`/100000),
         total_MA_MI_prop_low = 1-exp(-`MA mild illness lower CI`/100000),
         total_MA_MI_prop_up = 1-exp(-`MA mild illness upper CI`/100000),
         total_nMA_MI_prop = 1-exp(-`Non-MA mild illness rate`/100000),
         total_nMA_MI_prop_low = 1-exp(-`Non-MA mild illness lower CI`/100000),
         total_nMA_MI_prop_up = 1-exp(-`Non-MA mild illness upper CI`/100000)) %>%
  select(`Age (months)`,
         total_MI_prop,
         total_MI_prop_low,
         total_MI_prop_up,
         total_MA_MI_prop,
         total_MA_MI_prop_low,
         total_MA_MI_prop_up,
         total_nMA_MI_prop,
         total_nMA_MI_prop_low,
         total_nMA_MI_prop_up)

severe_illness_prop <- severe_illness_rate %>%
  mutate(total_SI_prop = 1-exp(-`Total severe illness rate`/100000),
         total_SI_prop_low = 1-exp(-`Total severe illness lower CI`/100000),
         total_SI_prop_up = 1-exp(-`Total severe illness upper CI`/100000),
         total_MA_SI_prop = 1-exp(-`MA severe illness rate`/100000),
         total_MA_SI_prop_low = 1-exp(-`MA severe illness lower CI`/100000),
         total_MA_SI_prop_up = 1-exp(-`MA severe illness upper CI`/100000),
         total_nMA_SI_prop = 1-exp(-`Non-MA severe illness rate`/100000),
         total_nMA_SI_prop_low = 1-exp(-`Non-MA severe illness lower CI`/100000),
         total_nMA_SI_prop_up = 1-exp(-`Non-MA severe illness upper CI`/100000)) %>%
  select(`Age (months)`,
         total_SI_prop,
         total_SI_prop_low,
         total_SI_prop_up,
         total_MA_SI_prop,
         total_MA_SI_prop_low,
         total_MA_SI_prop_up,
         total_nMA_SI_prop,
         total_nMA_SI_prop_low,
         total_nMA_SI_prop_up)

# Combine seroconversion data with illness data
# By age
mild_illness_progression <- conversion_formated %>%
  merge(mild_illness_prop, by.x = "age_bracket", by.y = "Age (months)") %>%
  mutate(total_mild_cases = median*total_MI_prop,
         total_MA_mild_cases = median*total_MA_MI_prop,
         total_nonMA_mild_cases = median*total_nMA_MI_prop) %>%
  select(age_bracket, age_months, agemid, birth_season, current_season, median, 
         total_mild_cases, total_MA_mild_cases, total_nonMA_mild_cases)

severe_illness_progression <- conversion_formated %>%
  merge(severe_illness_prop, by.x = "age_bracket", by.y = "Age (months)") %>%
  mutate(total_severe_cases = median*total_SI_prop,
         total_MA_severe_cases = median*total_MA_SI_prop,
         total_nonMA_severe_cases = median*total_nMA_SI_prop) %>%
  select(age_bracket, age_months, agemid, birth_season, current_season, median, 
         total_severe_cases, total_MA_severe_cases, total_nonMA_severe_cases)

# By season in the first year of life
mild_illness_progression_season <- mild_illness_progression %>%
  # get the contribution of each birth cohort to the proportion of sick children
  group_by(current_season, agemid, birth_season) %>%
  mutate(seasonal_contribution = case_when(birth_season == "spring" ~ 0.26*total_mild_cases, 
                                           birth_season == "summer" ~ 0.29*total_mild_cases,
                                           birth_season == "autumn" ~ 0.24*total_mild_cases, 
                                           birth_season == "winter" ~ 0.20*total_mild_cases
                                            )) %>%
  ungroup() %>%
  filter(agemid <= 365) %>% # look at first year of life only
  # get total proportion of sick children by season
  group_by(current_season) %>%
  mutate(mild_illness_season = sum(seasonal_contribution)) %>% 
  select(age_months, agemid, birth_season, current_season, seasonal_contribution, mild_illness_season)

severe_illness_progression_season <- severe_illness_progression %>%
  group_by(current_season, agemid, birth_season) %>%
  mutate(seasonal_contribution = case_when(birth_season == "spring" ~ 0.26*total_severe_cases,
                                           birth_season == "summer" ~ 0.29*total_severe_cases,
                                           birth_season == "autumn" ~ 0.24*total_severe_cases, 
                                           birth_season == "winter" ~ 0.20*total_severe_cases
  )) %>%
  ungroup() %>%
  filter(agemid <= 365) %>%
  group_by(current_season) %>%
  mutate(severe_illness_season = sum(seasonal_contribution)) %>%
  select(age_months, agemid, birth_season, current_season, seasonal_contribution, severe_illness_season)

# Plots
mild_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_mild_cases, colour = birth_season)) +
  labs(title = "Proportion of all mild illness", x = "Age (days)",
        y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() 

mild_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_MA_mild_cases, colour = birth_season)) +
  labs(title = "Proportion of medically assisted mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

mild_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_nonMA_mild_cases, colour = birth_season)) +
  labs(title = "Proportion of non-medically-assisted mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_severe_cases, colour = birth_season)) +
  labs(title = "Proportion of all severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() 

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_MA_severe_cases, colour = birth_season)) +
  labs(title = "Proportion of medically assisted severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_nonMA_severe_cases, colour = birth_season)) +
  labs(title = "Proportion of non-medically-assisted severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

# Check seasonal trends by season of birth
mild_illness_progression %>% filter (agemid <= 365) %>%
  ggplot(aes(x = agemid, y = total_mild_cases)) +
  geom_point(aes(colour = current_season)) + 
  geom_smooth() +
  labs(title = "Proportion of all mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Current season") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~birth_season)

severe_illness_progression %>% filter (agemid <= 365) %>%
  ggplot(aes(x = agemid, y = total_severe_cases)) +
  geom_point(aes(colour = current_season)) + 
  geom_smooth() +
  labs(title = "Proportion of all severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Current season") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~birth_season)

# Seasonal trends in the first year of life
mild_illness_progression_season %>% ggplot() +
  geom_point(aes(x = current_season, y = mild_illness_season)) +
  labs(title = "Proportion of all mild illness in the first year of life", x = "Season",
       y = "Proportion of illness") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression_season %>% ggplot() +
  geom_point(aes(x = current_season, y = severe_illness_season)) +
  labs(title = "Proportion of all severe illness in the first year of life", x = "Season",
       y = "Proportion of illness") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

# Save outputs
mild_illness_progression %>% write.csv(paste0(path_model, "Proportion mild illness by age.csv"), row.names = F)
severe_illness_progression %>% write.csv(paste0(path_model, "Proportion severe illness by age.csv"), row.names = F)
