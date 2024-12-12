# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)

# Define useful paths
path_paper <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/" # Where results from other studies are saved
path_model <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/ExtensioL/CSV files/Vaccination/" #Where model outputs are stored
  
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
conversion <- conversion_rate_summer %>% mutate(season = "summer") %>%
  rbind(
    conversion_rate_autumn %>% mutate(season = "autumn")
  ) %>%
  rbind(
    conversion_rate_winter %>% mutate(season = "winter")
  ) %>%
  rbind(
    conversion_rate_spring %>% mutate(season = "spring")
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
mild_illness_progression <- conversion_formated %>%
  merge(mild_illness_prop, by.x = "age_bracket", by.y = "Age (months)") %>%
  mutate(total_mild_cases = median*total_MI_prop,
         total_MA_mild_cases = median*total_MA_MI_prop,
         total_nonMA_mild_cases = median*total_nMA_MI_prop) %>%
  select(age_bracket, age_months, agemid, season, median, total_mild_cases,
         total_MA_mild_cases, total_nonMA_mild_cases)

severe_illness_progression <- conversion_formated %>%
  merge(severe_illness_prop, by.x = "age_bracket", by.y = "Age (months)") %>%
  mutate(total_severe_cases = median*total_SI_prop,
         total_MA_severe_cases = median*total_MA_SI_prop,
         total_nonMA_severe_cases = median*total_nMA_SI_prop) %>%
  select(age_bracket, age_months, agemid, season, median, total_severe_cases,
         total_MA_severe_cases, total_nonMA_severe_cases)

# Plots
mild_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_mild_cases, colour = season)) +
  labs(title = "Proportion of all mild illness", x = "Age (days)",
        y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() 

mild_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_MA_mild_cases, colour = season)) +
  labs(title = "Proportion of medically assisted mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

mild_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_nonMA_mild_cases, colour = season)) +
  labs(title = "Proportion of non-medically-assisted mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_severe_cases, colour = season)) +
  labs(title = "Proportion of all severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() 

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_MA_severe_cases, colour = season)) +
  labs(title = "Proportion of medically assisted severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

severe_illness_progression %>% ggplot() +
  geom_point(aes(x = agemid, y = total_nonMA_severe_cases, colour = season)) +
  labs(title = "Proportion of non-medically-assisted severe illness", x = "Age (days)",
       y = "Proportion of severe illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

# Save outputs
mild_illness_progression %>% write.csv(paste0(path_model, "Proportion mild illness by age.csv"), row.names = F)
severe_illness_progression %>% write.csv(paste0(path_model, "Proportion severe illness by age.csv"), row.names = F)

# Check by season of birth
mild_illness_progression %>% filter (agemid <= 365) %>%
  ggplot() +
  geom_point(aes(x = agemid, y = total_mild_cases, colour = season)) +
  labs(title = "Proportion of all mild illness", x = "Age (days)",
       y = "Proportion of mild illness", colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season)
