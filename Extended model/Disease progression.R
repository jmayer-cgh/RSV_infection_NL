# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)

# Define useful paths
path_paper <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/" # Where results from other studies are saved
path_model <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/mcstate/"  #Where model outputs are stored
  
# Read in RSV-illness rates
mild_illness_rate <- read_excel(paste0(path_paper,"RSV illness rates SA.xlsx"), sheet = "Mild illness")
severe_illness_rate <- read_excel(paste0(path_paper,"RSV illness rates SA.xlsx"), sheet = "Severe illness")

# Read in model estimate
conversion_rate <- read.csv(paste0(path_model, "incidence by age.csv")) 

# Convert into a long format
conversion <- conversion_rate %>% 
  select(age_midpoint, low95_sp, incidence_sp, up95_sp) %>%
  rename(low95 = "low95_sp", incidence = "incidence_sp", up95 = "up95_sp") %>%
  mutate(season_birth = "spring") %>%
  rbind(
    conversion_rate %>% 
      select(age_midpoint, low95_sm, incidence_sm, up95_sm) %>%
      rename(low95 = "low95_sm", incidence = "incidence_sm", up95 = "up95_sm") %>%
      mutate(season_birth = "summer")
  ) %>%
  rbind(
    conversion_rate %>% 
      select(age_midpoint, low95_au, incidence_au, up95_au) %>%
      rename(low95 = "low95_au", incidence = "incidence_au", up95 = "up95_au") %>%
      mutate(season_birth = "autumn")
  ) %>%
  rbind(
    conversion_rate %>% 
    select(age_midpoint, low95_wt, incidence_wt, up95_wt) %>%
    rename(low95 = "low95_wt", incidence = "incidence_wt", up95 = "up95_wt") %>%
    mutate(season_birth = "winter")
  )

# Convert ages to the same units
conversion_formated <- conversion %>% mutate(age_months = trunc(age_midpoint/30.4375)) %>% # turn age into months
  arrange(age_midpoint) # arrange by age instead of by season of birth

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
  group_by(age_midpoint, season_birth) %>%
  mutate(current_season = case_when(season_birth == "summer" & ((age_midpoint <= 30.41*3)          # season in which the cohort was born                    
                                    || ((age_midpoint >= 365) & (age_midpoint <= (365+30.41*3))) 
                                    || ((age_midpoint >= 2*365) & (age_midpoint <= (2*365+30.41*3))) 
                                    || ((age_midpoint >= 3*365) & (age_midpoint <= (3*365+30.41*3))) 
                                    || ((age_midpoint >= 4*365) & (age_midpoint <= (4*365+30.41*3))) 
                                    || ((age_midpoint >= 5*365) & (age_midpoint <= (5*365+30.41*3)))) ~ "summer",
                                    season_birth == "autumn" & ((age_midpoint <= 30.41*3)                              
                                    || ((age_midpoint >= 365) & (age_midpoint <= (365+30.41*3))) 
                                    || ((age_midpoint >= 2*365) & (age_midpoint <= (2*365+30.41*3))) 
                                    || ((age_midpoint >= 3*365) & (age_midpoint <= (3*365+30.41*3))) 
                                    || ((age_midpoint >= 4*365) & (age_midpoint <= (4*365+30.41*3))) 
                                    || ((age_midpoint >= 5*365) & (age_midpoint <= (5*365+30.41*3)))) ~ "autumn",
                                    season_birth == "winter" & ((age_midpoint <= 30.41*3)                              
                                    || ((age_midpoint >= 365) & (age_midpoint <= (365+30.41*3))) 
                                    || ((age_midpoint >= 2*365) & (age_midpoint <= (2*365+30.41*3))) 
                                    || ((age_midpoint >= 3*365) & (age_midpoint <= (3*365+30.41*3))) 
                                    || ((age_midpoint >= 4*365) & (age_midpoint <= (4*365+30.41*3))) 
                                    || ((age_midpoint >= 5*365) & (age_midpoint <= (5*365+30.41*3)))) ~ "winter",
                                    season_birth == "spring" & ((age_midpoint <= 30.41*3)                              
                                    || ((age_midpoint >= 365) & (age_midpoint <= (365+30.41*3))) 
                                    || ((age_midpoint >= 2*365) & (age_midpoint <= (2*365+30.41*3))) 
                                    || ((age_midpoint >= 3*365) & (age_midpoint <= (3*365+30.41*3))) 
                                    || ((age_midpoint >= 4*365) & (age_midpoint <= (4*365+30.41*3))) 
                                    || ((age_midpoint >= 5*365) & (age_midpoint <= (5*365+30.41*3)))) ~ "spring",
                                    
                                    season_birth == "summer" & ((age_midpoint > 30.41*3 & age_midpoint <= 30.41*6 ) # season following the one the cohort was born in
                                    || ((age_midpoint > 365 + 30.41*3) & (age_midpoint <= (365+30.41*6))) 
                                    || ((age_midpoint > 2*365 + 30.41*3) & (age_midpoint <= (2*365+30.41*6))) 
                                    || ((age_midpoint > 3*365 + 30.41*3) & (age_midpoint <= (3*365+30.41*6))) 
                                    || ((age_midpoint > 4*365 + 30.41*3) & (age_midpoint <= (4*365+30.41*6))) 
                                    || ((age_midpoint > 5*365 + 30.41*3) & (age_midpoint <= (5*365+30.41*6)))) ~ "autumn",
                                    season_birth == "autumn" & ((age_midpoint > 30.41*3 & age_midpoint <= 30.41*6 ) 
                                    || ((age_midpoint > 365 + 30.41*3) & (age_midpoint <= (365+30.41*6))) 
                                    || ((age_midpoint > 2*365 + 30.41*3) & (age_midpoint <= (2*365+30.41*6))) 
                                    || ((age_midpoint > 3*365 + 30.41*3) & (age_midpoint <= (3*365+30.41*6))) 
                                    || ((age_midpoint > 4*365 + 30.41*3) & (age_midpoint <= (4*365+30.41*6))) 
                                    || ((age_midpoint > 5*365 + 30.41*3) & (age_midpoint <= (5*365+30.41*6)))) ~ "winter",
                                    season_birth == "winter" & ((age_midpoint > 30.41*3 & age_midpoint <= 30.41*6 ) 
                                    || ((age_midpoint > 365 + 30.41*3) & (age_midpoint <= (365+30.41*6))) 
                                    || ((age_midpoint > 2*365 + 30.41*3) & (age_midpoint <= (2*365+30.41*6))) 
                                    || ((age_midpoint > 3*365 + 30.41*3) & (age_midpoint <= (3*365+30.41*6))) 
                                    || ((age_midpoint > 4*365 + 30.41*3) & (age_midpoint <= (4*365+30.41*6))) 
                                    || ((age_midpoint > 5*365 + 30.41*3) & (age_midpoint <= (5*365+30.41*6)))) ~ "spring",
                                    season_birth == "spring" & ((age_midpoint > 30.41*3 & age_midpoint <= 30.41*6 ) 
                                    || ((age_midpoint > 365 + 30.41*3) & (age_midpoint <= (365+30.41*6))) 
                                    || ((age_midpoint > 2*365 + 30.41*3) & (age_midpoint <= (2*365+30.41*6))) 
                                    || ((age_midpoint > 3*365 + 30.41*3) & (age_midpoint <= (3*365+30.41*6))) 
                                    || ((age_midpoint > 4*365 + 30.41*3) & (age_midpoint <= (4*365+30.41*6))) 
                                    || ((age_midpoint > 5*365 + 30.41*3) & (age_midpoint <= (5*365+30.41*6)))) ~ "summer",
                                    
                                    season_birth == "summer" & ((age_midpoint > 30.41*6 & age_midpoint <= 30.41*9 )  # two seasons after birth
                                    || ((age_midpoint > 365+30.41*6) & (age_midpoint <= (365+30.41*9))) 
                                    || ((age_midpoint > 2*365 + 30.41*6) & (age_midpoint <= (2*365+30.41*9))) 
                                    || ((age_midpoint > 3*365 + 30.41*6) & (age_midpoint <= (3*365+30.41*9))) 
                                    || ((age_midpoint > 4*365 + 30.41*6) & (age_midpoint <= (4*365+30.41*9))) 
                                    || ((age_midpoint > 5*365 + 30.41*6) & (age_midpoint <= (5*365+30.41*9)))) ~ "winter",
                                    season_birth == "autumn" & ((age_midpoint > 30.41*6 & age_midpoint <= 30.41*9 )  
                                    || ((age_midpoint > 365+30.41*6) & (age_midpoint <= (365+30.41*9))) 
                                    || ((age_midpoint > 2*365 + 30.41*6) & (age_midpoint <= (2*365+30.41*9))) 
                                    || ((age_midpoint > 3*365 + 30.41*6) & (age_midpoint <= (3*365+30.41*9))) 
                                    || ((age_midpoint > 4*365 + 30.41*6) & (age_midpoint <= (4*365+30.41*9))) 
                                    || ((age_midpoint > 5*365 + 30.41*6) & (age_midpoint <= (5*365+30.41*9)))) ~ "spring",
                                    season_birth == "winter" & ((age_midpoint > 30.41*6 & age_midpoint <= 30.41*9 )  
                                    || ((age_midpoint > 365+30.41*6) & (age_midpoint <= (365+30.41*9))) 
                                    || ((age_midpoint > 2*365 + 30.41*6) & (age_midpoint <= (2*365+30.41*9))) 
                                    || ((age_midpoint > 3*365 + 30.41*6) & (age_midpoint <= (3*365+30.41*9))) 
                                    || ((age_midpoint > 4*365 + 30.41*6) & (age_midpoint <= (4*365+30.41*9))) 
                                    || ((age_midpoint > 5*365 + 30.41*6) & (age_midpoint <= (5*365+30.41*9)))) ~ "summer",
                                    season_birth == "spring" & ((age_midpoint > 30.41*6 & age_midpoint <= 30.41*9 )  
                                    || ((age_midpoint > 365+30.41*6) & (age_midpoint <= (365+30.41*9))) 
                                    || ((age_midpoint > 2*365 + 30.41*6) & (age_midpoint <= (2*365+30.41*9))) 
                                    || ((age_midpoint > 3*365 + 30.41*6) & (age_midpoint <= (3*365+30.41*9))) 
                                    || ((age_midpoint > 4*365 + 30.41*6) & (age_midpoint <= (4*365+30.41*9))) 
                                    || ((age_midpoint > 5*365 + 30.41*6) & (age_midpoint <= (5*365+30.41*9)))) ~ "autumn",
                                    
                                    season_birth == "summer" & ((age_midpoint > 30.41*9 & age_midpoint <= 30.41*12 ) # 3 seasons after birth                    
                                    || ((age_midpoint > 365+30.41*9) & (age_midpoint <= (365+30.41*12))) 
                                    || ((age_midpoint > 2*365 + 30.41*9) & (age_midpoint <= (2*365+30.41*12))) 
                                    || ((age_midpoint > 3*365 + 30.41*9) & (age_midpoint <= (3*365+30.41*12))) 
                                    || ((age_midpoint > 4*365 + 30.41*9) & (age_midpoint <= (4*365+30.41*12))) 
                                    || ((age_midpoint > 5*365 + 30.41*9) & (age_midpoint <= (5*365+30.41*12)))) ~ "spring",
                                    season_birth == "autumn" & ((age_midpoint > 30.41*9 & age_midpoint <= 30.41*12 )                     
                                    || ((age_midpoint > 365+30.41*9) & (age_midpoint <= (365+30.41*12))) 
                                    || ((age_midpoint > 2*365 + 30.41*9) & (age_midpoint <= (2*365+30.41*12))) 
                                    || ((age_midpoint > 3*365 + 30.41*9) & (age_midpoint <= (3*365+30.41*12))) 
                                    || ((age_midpoint > 4*365 + 30.41*9) & (age_midpoint <= (4*365+30.41*12))) 
                                    || ((age_midpoint > 5*365 + 30.41*9) & (age_midpoint <= (5*365+30.41*12)))) ~ "summer",
                                    season_birth == "winter" & ((age_midpoint > 30.41*9 & age_midpoint <= 30.41*12 )                     
                                    || ((age_midpoint > 365+30.41*9) & (age_midpoint <= (365+30.41*12))) 
                                    || ((age_midpoint > 2*365 + 30.41*9) & (age_midpoint <= (2*365+30.41*12))) 
                                    || ((age_midpoint > 3*365 + 30.41*9) & (age_midpoint <= (3*365+30.41*12))) 
                                    || ((age_midpoint > 4*365 + 30.41*9) & (age_midpoint <= (4*365+30.41*12))) 
                                    || ((age_midpoint > 5*365 + 30.41*9) & (age_midpoint <= (5*365+30.41*12)))) ~ "autumn",
                                    season_birth == "spring" & ((age_midpoint > 30.41*9 & age_midpoint <= 30.41*12 )                   
                                    || ((age_midpoint > 365+30.41*9) & (age_midpoint <= (365+30.41*12))) 
                                    || ((age_midpoint > 2*365 + 30.41*9) & (age_midpoint <= (2*365+30.41*12))) 
                                    || ((age_midpoint > 3*365 + 30.41*9) & (age_midpoint <= (3*365+30.41*12))) 
                                    || ((age_midpoint > 4*365 + 30.41*9) & (age_midpoint <= (4*365+30.41*12))) 
                                    || ((age_midpoint > 5*365 + 30.41*9) & (age_midpoint <= (5*365+30.41*12)))) ~ "winter"
                                    
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
  mutate(total_mild_cases_prop = incidence*total_MI_prop,
         total_MA_mild_cases_prop = incidence*total_MA_MI_prop,
         total_nonMA_mild_cases_prop = incidence*total_nMA_MI_prop) %>%
  select(age_bracket, age_months, age_midpoint, season_birth, current_season, incidence, 
         total_mild_cases_prop, total_MA_mild_cases_prop, total_nonMA_mild_cases_prop)

severe_illness_progression <- conversion_formated %>%
  merge(severe_illness_prop, by.x = "age_bracket", by.y = "Age (months)") %>%
  mutate(total_severe_cases_prop = incidence*total_SI_prop,
         total_MA_severe_cases_prop = incidence*total_MA_SI_prop,
         total_nonMA_severe_cases_prop = incidence*total_nMA_SI_prop) %>%
  select(age_bracket, age_months, age_midpoint, season_birth, current_season, incidence, 
         total_severe_cases_prop, total_MA_severe_cases_prop, total_nonMA_severe_cases_prop)

# By season in the first year of life
mild_illness_progression_season <- mild_illness_progression %>%
  # get the contribution of each birth cohort to the proportion of sick children
  group_by(current_season, age_midpoint, season_birth) %>%
  mutate(seasonal_contribution = case_when(season_birth == "spring" ~ 0.26*total_mild_cases_prop, 
                                           season_birth == "summer" ~ 0.29*total_mild_cases_prop,
                                           season_birth == "autumn" ~ 0.24*total_mild_cases_prop, 
                                           season_birth == "winter" ~ 0.20*total_mild_cases_prop
                                            )) %>%
  ungroup() %>%
  filter(age_midpoint <= 365) %>% # look at first year of life only
  # get total proportion of sick children by season
  group_by(current_season) %>%
  mutate(mild_illness_season = sum(seasonal_contribution)) %>% 
  select(age_months, age_midpoint, season_birth, current_season, seasonal_contribution, mild_illness_season)

severe_illness_progression_season <- severe_illness_progression %>%
  group_by(current_season, age_midpoint, season_birth) %>%
  mutate(seasonal_contribution = case_when(season_birth == "spring" ~ 0.26*total_severe_cases_prop,
                                           season_birth == "summer" ~ 0.29*total_severe_cases_prop,
                                           season_birth == "autumn" ~ 0.24*total_severe_cases_prop, 
                                           season_birth == "winter" ~ 0.20*total_severe_cases_prop
  )) %>%
  ungroup() %>%
  filter(age_midpoint <= 365) %>%
  group_by(current_season) %>%
  mutate(severe_illness_season = sum(seasonal_contribution)) %>%
  select(age_months, age_midpoint, season_birth, current_season, seasonal_contribution, severe_illness_season)

# Plots
mild_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = total_mild_cases_prop, colour = season_birth)) +
  labs(title = "Proportion of all mild illness", x = "Age (days)",
        y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() 

mild_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = total_MA_mild_cases_prop, colour = season_birth)) +
  labs(title = "Proportion of medically assisted mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

mild_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = total_nonMA_mild_cases_prop, colour = season_birth)) +
  labs(title = "Proportion of non-medically-assisted mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = total_severe_cases_prop, colour = season_birth)) +
  labs(title = "Proportion of all severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() 

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = total_MA_severe_cases_prop, colour = season_birth)) +
  labs(title = "Proportion of medically assisted severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = total_nonMA_severe_cases_prop, colour = season_birth)) +
  labs(title = "Proportion of non-medically-assisted severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

# Check seasonal trends by season of birth
mild_illness_progression %>% filter (age_midpoint <= 365) %>%
  ggplot(aes(x = age_midpoint, y = total_mild_cases_prop)) +
  geom_point(aes(colour = current_season)) + 
  geom_smooth() +
  labs(title = "Proportion of all mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Current season") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth)

severe_illness_progression %>% filter (age_midpoint <= 365) %>%
  ggplot(aes(x = age_midpoint, y = total_severe_cases_prop)) +
  geom_point(aes(colour = current_season)) + 
  geom_smooth() +
  labs(title = "Proportion of all severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Current season") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth)

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
