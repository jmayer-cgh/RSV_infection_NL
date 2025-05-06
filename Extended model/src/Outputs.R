# The whole script takes just over 3 hours (3:07)

# Housekeeping
rm(list=ls())

# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
options(dplyr.summarise.inform = FALSE)

# Path to files
path_output <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"  #Where model outputs are stored
path_pop <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/German population estimates/"
path_code <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Extended model/src/"
path_paper <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/" # Where results from other studies are saved
path_model <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"  #Where model outputs are stored

# --------- Read in files ------------------------------------------------------
# Read in birth numbers for 2023
suppressMessages(births_de <- read_excel(paste0(path_pop, "Births.xlsx"))) 
colnames(births_de)[2] <- "Month"
colnames(births_de)[5] <- "Value"
births_de <- births_de %>% select(Month, Value) %>% 
  filter(grepl(".0", Value)) %>%
  mutate(Value = as.numeric(Value))

births_de <- births_de %>%
  mutate(season = case_when(
    Month == "January" | Month == "February" | Month == "December" ~ "winter",
    Month == "March" | Month == "April" | Month == "May" ~ "spring",
    Month == "June" | Month == "July" | Month == "August" ~ "summer",
    Month == "September" | Month == "October" | Month == "November" ~ "autumn"))

# Get % births by season
births_de <- births_de %>% 
  mutate(N = sum(Value)) %>%
  group_by(season, N) %>%
  reframe(total = sum(Value),
          prop = total/N) %>%
  unique()

# Read in RSV-illness estimates
mild_illness <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Mild illness numbers") %>% 
  janitor::clean_names()
severe_illness <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Severe illness numbers") %>% 
  janitor::clean_names()
# ------------------------------------------------------------------------------

# --------- Define our functions ----------------------------------------------
# Newly seroconverted at a given age
new_seroconv_age <- function (index){
  converted_all_rand <- converted_all[,index]
  converted_sp_rand <- converted_sp[,index]
  converted_sm_rand <- converted_sm[,index]
  converted_au_rand <- converted_au[,index]
  converted_wt_rand <- converted_wt[,index]
  
  converted <- data.frame(age_midpoint = incidence_data_season_wide$time,
                          sero_all = converted_all_rand,
                          sero_sp = converted_sp_rand,
                          sero_sm = converted_sm_rand,
                          sero_au = converted_au_rand,
                          sero_wt = converted_wt_rand)
  
  incidence_df <- data.frame(age_midpoint = converted$age_midpoint,
                             incidence = c(0, diff(converted$sero_sp)),
                             season_birth = "spring") %>%
    rbind(
      data.frame(age_midpoint = converted$age_midpoint,
                 incidence = c(0, diff(converted$sero_sm)),
                 season_birth = "summer")
    ) %>%
    rbind(
      data.frame(age_midpoint = converted$age_midpoint,
                 incidence = c(0, diff(converted$sero_au)),
                 season_birth = "autumn")
    ) %>%
    rbind(
      data.frame(age_midpoint = converted$age_midpoint,
                 incidence = c(0, diff(converted$sero_wt)),
                 season_birth = "winter")
    ) %>%
    rbind(
      data.frame(age_midpoint = converted$age_midpoint,
                 incidence = c(0, diff(converted$sero_all)),
                 season_birth = "all")
    )
}

# Get proportion of cases that are mild, severe, MA, non-MA
illness_type <- mild_illness %>% 
  merge(severe_illness, by = "age_months") %>%
  group_by(age_months) %>% # total number of cases in each category first
  summarise(total_cases = total_mild_illness_number + total_severe_illness_number,
            total_ma_cases = ma_mild_illness_number + ma_severe_illness_number,
            total_non_ma_cases = non_ma_mild_illness_number + non_ma_severe_illness_number,
            
            # prop of each type of mild cases
            prop_mild_cases = total_mild_illness_number / total_cases,
            prop_ma_mild_cases = ma_mild_illness_number / total_cases,
            prop_non_ma_mild_cases = non_ma_mild_illness_number / total_cases,
            
            # prop of each type of severe cases
            prop_severe_cases = total_severe_illness_number / total_cases,
            prop_ma_severe_cases = ma_severe_illness_number / total_cases,
            prop_non_ma_severe_cases = non_ma_severe_illness_number / total_cases,
            
            # prop of cases by MA status
            prop_ma_cases = total_ma_cases / total_cases,
            prop_non_ma_cases = total_non_ma_cases / total_cases) %>%
  ungroup()

# Read in model estimates
# conversion <- read.csv(paste0(path_model, "incidence by age.csv")) 

# ----- Data processing
progression_format <- function (conversion, illness_type) {
  # Convert ages to the same units
  conversion_formated <- conversion %>% mutate(age_months = trunc(age_midpoint/30)) %>% # turn age into months
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
  
  # Combine seroconversion data with illness data
  mild_illness_progression <- conversion_formated %>%
    merge(illness_type, by.x = "age_bracket", by.y = "age_months") %>%
    mutate(total_mild_cases_prop = incidence * prop_mild_cases,
           total_ma_mild_cases_prop = incidence * prop_ma_mild_cases,
           total_non_ma_mild_cases_prop = incidence * prop_non_ma_mild_cases) %>%
    select(age_bracket, age_months, age_midpoint, season_birth, incidence, 
           total_mild_cases_prop, total_ma_mild_cases_prop, total_non_ma_mild_cases_prop)
  
  severe_illness_progression <- conversion_formated %>%
    merge(illness_type, by.x = "age_bracket", by.y = "age_months") %>%
    mutate(total_severe_cases_prop = incidence * prop_severe_cases,
           total_ma_severe_cases_prop = incidence * prop_ma_severe_cases,
           total_non_ma_severe_cases_prop = incidence * prop_non_ma_severe_cases) %>%
    select(age_bracket, age_months, age_midpoint, season_birth, incidence, 
           total_severe_cases_prop, total_ma_severe_cases_prop, total_non_ma_severe_cases_prop)
  
  return (severe_illness_progression)
}

cases_format <- function (severe_illness, mild_illness) {
  # Restrict to < 1 year
  severe_illness_u1 <- severe_illness %>% filter (age_months <= 12)
  mild_illness_u1 <- mild_illness %>% filter (age_months <= 12)
  
  # Add pop size in <1 y.o by season
  severe_cases_u1_de <- severe_illness_u1 %>%
    merge(births_de %>% select(total, season),
          by.x = "season_birth" , by.y = "season")
  
  mild_cases_u1_de <- mild_illness_u1 %>%
    merge(births_de %>% select(total, season),
          by.x = "season_birth" , by.y = "season")
  
  # Get number of cases
  severe_cases_u1_de <- severe_cases_u1_de %>%
    mutate(n_severe = total_severe_cases_prop * total,
           n_ma_severe = total_ma_severe_cases_prop * total,
           n_non_ma_severe = total_non_ma_severe_cases_prop * total)
  
  mild_cases_u1_de <- mild_cases_u1_de %>%
    mutate(n_mild = total_mild_cases_prop * total,
           n_ma_mild = total_ma_mild_cases_prop * total,
           n_non_ma_mild = total_non_ma_mild_cases_prop * total)
  
  # Get totals
  severe_cases_u1_de <- severe_cases_u1_de %>% rbind(
    severe_cases_u1_de %>% summarise (season_birth = "all",
                                      age_bracket = "all",
                                      age_months = NA,
                                      age_midpoint = NA,
                                      incidence = NA,
                                      total_severe_cases_prop = NA,
                                      total_ma_severe_cases_prop = NA,
                                      total_non_ma_severe_cases_prop = NA,
                                      total = births_de$N[1],
                                      n_severe = sum(n_severe),
                                      n_ma_severe = sum(n_ma_severe),
                                      n_non_ma_severe = sum(n_non_ma_severe))
  )
  
  mild_cases_u1_de <- mild_cases_u1_de %>% rbind(
    mild_cases_u1_de %>% summarise (season_birth = "all",
                                    age_bracket = "all",
                                    age_months = NA,
                                    age_midpoint = NA,
                                    incidence = NA,
                                    total_mild_cases_prop = NA,
                                    total_ma_mild_cases_prop = NA,
                                    total_non_ma_mild_cases_prop = NA,
                                    total = births_de$N[1],
                                    n_mild = sum(n_mild),
                                    n_ma_mild = sum(n_ma_mild),
                                    n_non_ma_mild = sum(n_non_ma_mild))
  )
  
  # Scale the numbers to have them match Fabienne's paper
  # About 22,000 hospitalisations in 2019, of which about 12% in year 1 --> 2,640 
  # Model predicts about 65,000 hospitalisations --> scale by 0.04
  severe_cases_u1_de <- severe_cases_u1_de %>% 
    mutate (n_severe_scaled = n_severe * 0.04,
            n_ma_severe_scaled = n_ma_severe * 0.04,
            n_non_ma_severe_scaled = n_non_ma_severe * 0.04)
  
  mild_cases_u1_de <- mild_cases_u1_de %>% 
    mutate (n_mild_scaled = n_mild * 0.04,
            n_ma_mild_scaled = n_ma_mild * 0.04,
            n_non_ma_mild_scaled = n_non_ma_mild * 0.04)
  
  return (severe_cases_u1_de)
  
}

# ----- Interventions
hosp_vacc <- function (cases, VE_distribution){
  # Add the estimated VE to the number of cases
  hosp_prevented_vacc <- cases %>% filter(!is.na(age_months)) %>%
    select(season_birth, age_bracket, age_months, age_midpoint,
           n_ma_severe_scaled) %>%
    merge(
      VE_distribution %>% filter(group == "severe") %>%
        select(t, VE_t),
      by.x = "age_midpoint", by.y = "t"
    )
  
  # Get number of prevented cases
  hosp_prevented_vacc <- hosp_prevented_vacc %>% 
    group_by(season_birth, age_bracket, age_months, age_midpoint, VE_t) %>%
    summarise(n_ma_severe_vacc = n_ma_severe_scaled * (1-VE_t), # number of cases despite vaccination
              n_ma_severe_averted = n_ma_severe_scaled * VE_t)
  
  hosp_prevented_vacc <- hosp_prevented_vacc %>% 
    rbind(
      hosp_prevented_vacc %>% 
        group_by(season_birth) %>%
        summarise(age_bracket = 'all', 
                  age_months = NA,
                  age_midpoint = NA, 
                  n_ma_severe_vacc = sum(n_ma_severe_vacc),
                  n_ma_severe_averted = sum(n_ma_severe_averted))
    )
  return (hosp_prevented_vacc)
}


hosp_mAB <- function (hosp_prevented_vacc) {
  nirs_eff_hosp <- runif(1, 0.623, 0.852)
  
  hosp_prevented <- hosp_prevented_vacc %>% filter(!is.na(age_months)) %>%
    select(season_birth, age_bracket, age_months, age_midpoint,
           n_ma_severe_vacc, n_ma_severe_averted) %>%
    mutate(n_ma_severe_nirs = n_ma_severe_vacc * (1-nirs_eff_hosp),
           n_ma_severe_averted_nirs = n_ma_severe_vacc * nirs_eff_hosp)
  
  # Get number of cases averted or not by giving nirvesimab at a given age
  hosp_prevented_age <- hosp_prevented %>% group_by(season_birth) %>%
    # averted hospitalisations
    summarise (prev_u1 = sum(n_ma_severe_averted_nirs),
               prev_1 = sum(n_ma_severe_averted_nirs[age_months >= 1], na.rm = T),
               prev_3 = sum(n_ma_severe_averted_nirs[age_months >= 3], na.rm = T),
               prev_4 = sum(n_ma_severe_averted_nirs[age_months >= 4], na.rm = T),
               prev_7 = sum(n_ma_severe_averted_nirs[age_months >= 7], na.rm = T),
               prev_8 = sum(n_ma_severe_averted_nirs[age_months >= 8], na.rm = T),
               prev_10 = sum(n_ma_severe_averted_nirs[age_months >= 10], na.rm = T),
               prev_11 = sum(n_ma_severe_averted_nirs[age_months >= 11], na.rm = T),
               # prevented hospitalisations when immunising spring cohort only
               prev_spring_10 = sum(n_ma_severe_averted_nirs[age_months >= 10 & season_birth == "spring"], na.rm = T),
               prev_summer_4 = sum(n_ma_severe_averted_nirs[age_months >= 4 & season_birth == "summer"], na.rm = T),
               # remaining hospitalisations
               hosp_u1 = sum(n_ma_severe_nirs),
               hosp_1 = sum(n_ma_severe_nirs[age_months >= 1], na.rm = T) + sum(n_ma_severe_vacc[age_months < 1], na.rm = T),
               hosp_3 = sum(n_ma_severe_nirs[age_months >= 3], na.rm = T) + sum(n_ma_severe_vacc[age_months < 3], na.rm = T),
               hosp_4 = sum(n_ma_severe_nirs[age_months >= 4], na.rm = T) + sum(n_ma_severe_vacc[age_months < 4], na.rm = T),
               hosp_7 = sum(n_ma_severe_nirs[age_months >= 7], na.rm = T) + sum(n_ma_severe_vacc[age_months < 7], na.rm = T),
               hosp_8 = sum(n_ma_severe_nirs[age_months >= 8], na.rm = T) + sum(n_ma_severe_vacc[age_months < 8], na.rm = T),
               hosp_10 = sum(n_ma_severe_nirs[age_months >= 10], na.rm = T) + sum(n_ma_severe_vacc[age_months < 10], na.rm = T),
               hosp_11 = sum(n_ma_severe_nirs[age_months >= 11], na.rm = T) + sum(n_ma_severe_vacc[age_months < 11], na.rm = T),
    ) %>%
    ungroup() %>%
    # hospitalisations when immunising spring cohort only or summer cohort only (hosp. are averted in the spring or summer cohort only)
    mutate(hosp_spring_10 = case_when (season_birth == "spring" ~ hosp_10,
                                       TRUE ~ prev_10 + hosp_10),
           hosp_summer_4 = case_when (season_birth == "summer" ~ hosp_4,
                                      TRUE ~ prev_4 + hosp_4)) %>%
    rbind(
      data.frame(season_birth = "all",
                 # averted hospitalisations
                 prev_u1 = sum(hosp_prevented$n_ma_severe_averted_nirs),
                 prev_1 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 1], na.rm = T),
                 prev_3 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 3], na.rm = T),
                 prev_4 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 4], na.rm = T),
                 prev_7 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 7], na.rm = T),
                 prev_8 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 8], na.rm = T),
                 prev_10 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 10], na.rm = T),
                 prev_11 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 11], na.rm = T),
                 # prevented hospitalisations when immunising spring or summer cohort only
                 prev_spring_10 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 10 & hosp_prevented$season_birth == "spring"], na.rm = T),
                 prev_summer_4 = sum(hosp_prevented$n_ma_severe_averted_nirs[hosp_prevented$age_months >= 4 & hosp_prevented$season_birth == "summer"], na.rm = T),
                 # remaining hospitalisations
                 hosp_u1 = sum(hosp_prevented$n_ma_severe_nirs),
                 hosp_1 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 1], na.rm = T),
                 hosp_3 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 3], na.rm = T),
                 hosp_4 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 4], na.rm = T),
                 hosp_7 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 7], na.rm = T),
                 hosp_8 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 8], na.rm = T),
                 hosp_10 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 10], na.rm = T),
                 hosp_11 = sum(hosp_prevented$n_ma_severe_nirs[hosp_prevented$age_months >= 11], na.rm = T),
                 hosp_spring_10 = 0,
                 hosp_summer_4 = 0) %>%# placeholder
        ungroup()
    ) %>%
    mutate(hosp_spring_10 = case_when (season_birth == "all" ~ sum(hosp_spring_10),
                                       TRUE ~ hosp_spring_10),
           hosp_summer_4 = case_when (season_birth == "all" ~ sum(hosp_summer_4),
                                      TRUE ~ hosp_summer_4))
  
  return (hosp_prevented_age)
}

hosp_intervention <- function (cases, hosp_prevented_vacc, hosp_prevented_age){
  # Get the numbers
  total_hosp_intervention <- cases %>%
    filter(season_birth == "all") %>%
    select(season_birth, age_bracket, n_ma_severe_scaled) %>%
    rename(n_hospitalisations = n_ma_severe_scaled) %>%
    mutate(prevented_hospitalisations = 0,
           intervention = "No immunisation") %>%
    rbind(
      hosp_prevented_vacc %>% filter(age_bracket == "all") %>%
        ungroup() %>%
        select(season_birth, age_bracket, n_ma_severe_vacc, n_ma_severe_averted) %>%
        summarise(season_birth = "all",
                  age_bracket = "all",
                  n_hospitalisations = sum(n_ma_severe_vacc),
                  prevented_hospitalisations = sum(n_ma_severe_averted),
                  intervention = "Maternal vaccination")
    ) %>%
    rbind(
      hosp_prevented_age %>% filter(season_birth == "all") %>%
        ungroup() %>%
        select(hosp_spring_10, prev_spring_10) %>%
        summarise(season_birth = "all",
                  age_bracket = "all",
                  n_hospitalisations = hosp_spring_10,
                  prevented_hospitalisations = prev_spring_10,
                  intervention = "Nirsevimab for spring births")
    ) %>%
    rbind(
      hosp_prevented_age %>% filter(season_birth == "all") %>%
        ungroup() %>%
        select(hosp_summer_4, prev_summer_4) %>%
        summarise(season_birth = "all",
                  age_bracket = "all",
                  n_hospitalisations = hosp_summer_4,
                  prevented_hospitalisations = prev_summer_4,
                  intervention = "Nirsevimab for summer births")
    ) %>%
    rbind(
      hosp_prevented_age %>% filter(season_birth == "all") %>%
        ungroup() %>%
        select(hosp_summer_4, prev_summer_4, hosp_spring_10, prev_spring_10) %>%
        summarise(season_birth = "all",
                  age_bracket = "all",
                  n_hospitalisations = 0, # we need to calculate this later
                  prevented_hospitalisations = prev_summer_4 + prev_spring_10,
                  intervention = "Nirsevimab for spring and summer births")
    ) %>% mutate(intervention = factor(intervention, 
                                       levels = c("No immunisation", "Maternal vaccination", 
                                                  "Nirsevimab for spring births",
                                                  "Nirsevimab for summer births",
                                                  "Nirsevimab for spring and summer births")))
  
  total_hosp_intervention <- total_hosp_intervention %>%
    mutate(n_hospitalisations = case_when (intervention == "Nirsevimab for spring and summer births" ~
                                             n_hospitalisations[intervention == "Maternal vaccination"] - prevented_hospitalisations[intervention == "Nirsevimab for spring and summer births"],
                                           TRUE ~ n_hospitalisations))
  
  return (total_hosp_intervention)
}
# -----------------------------------------------------------------------------

# --------- Run the models once -----------------------------------------------
source(paste0(path_code, "VE estimates.R")) # takes about 10 min 
source(paste0(path_code, "Seroconversion fit.R")) # takes about 3 hours 
# -----------------------------------------------------------------------------

# Combine the estimates
# Define where we're going to store the outputs
total_hosp_intervention <- list()

for (i in 1:100){
  # ------ Calculate number of RSV cases by age in Germany ------------
  # Proportion of newly seroconverted children at a given age
  # Take one random iteration of the model
  rand <- sample(1:10000, 1)
  # And extract the values for it
  conversion <- new_seroconv_age (rand)

  # Proportion of children with RSV at a given age (seroconversion to disease)
  severe_illness <- progression_format (conversion, illness_type)
  mild_illness <- read.csv(paste0(path_output, "Proportion mild illness by age.csv")) %>%
    select (!current_season) # Change this later
  # Number of cases in Germany
  severe_illness_de <- cases_format(severe_illness, mild_illness)
  
  # ------ Maternal vaccination -----------------------------------------------
  # Take one random iteration of the model
  draw <- sample(1:10000, 1)
  # And extract the values for it
  VE_distribution <- model_ve_runs_t %>% filter(iter == draw)
  
  # Calculate number of hospitalisations despite the vaccine
  hosp_prevented_vacc <- hosp_vacc(severe_illness_de, VE_distribution)
  
  # Calculate number of hospitalisations despite the vaccine and the mAB
  hosp_prevented_age <- hosp_mAB(hosp_prevented_vacc)
  
  # Get total hospitalisations by intervention
  total_hosp_intervention[[i]] <- hosp_intervention(severe_illness_de, hosp_prevented_vacc, hosp_prevented_age)
  total_hosp_intervention[[i]]$iter <- rep(i, each=5) # save index in case we want to check one iteration

  
}
  
total_hosp_intervention_df <- do.call("rbind", total_hosp_intervention)

# Compute 95% CI
total_hosp_intervention_int <- total_hosp_intervention_df %>% 
  group_by(intervention) %>%
  summarise(hosp_low_95 = quantile(n_hospitalisations, 0.05),
            hosp_median = quantile(n_hospitalisations, 0.5),
            hosp_up_95 = quantile(n_hospitalisations, 0.95),
            prev_low_95 = quantile(prevented_hospitalisations, 0.05),
            prev_median = quantile(prevented_hospitalisations, 0.5),
            prev_up_95 = quantile(prevented_hospitalisations, 0.95)) %>%
  ungroup()

# Get NNV
nnv <- total_hosp_intervention_int %>% select(intervention, prev_low_95, 
                                              prev_median, prev_up_95) %>%
  mutate(NNV_low_95 = case_when (intervention == "No immusination" ~ NA,
                                 intervention == "Maternal vaccination" ~ births_de$N[1]/prev_low_95,
                                 intervention == "Nirsevimab for spring births" ~ births_de$total[births_de$season == "spring"]/prev_low_95,
                                 intervention == "Nirsevimab for summer births" ~ births_de$total[births_de$season == "summer"]/prev_low_95,
                                 intervention == "Nirsevimab for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_low_95),
         NNV_median = case_when (intervention == "No immusination" ~ NA,
                                 intervention == "Maternal vaccination" ~ births_de$N[1]/prev_median,
                                 intervention == "Nirsevimab for spring births" ~ births_de$total[births_de$season == "spring"]/prev_median,
                                 intervention == "Nirsevimab for summer births" ~ births_de$total[births_de$season == "summer"]/prev_median,
                                 intervention == "Nirsevimab for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_median),
         NNV_up_95 = case_when (intervention == "No immusination" ~ NA,
                                 intervention == "Maternal vaccination" ~ births_de$N[1]/prev_up_95,
                                 intervention == "Nirsevimab for spring births" ~ births_de$total[births_de$season == "spring"]/prev_up_95,
                                 intervention == "Nirsevimab for summer births" ~ births_de$total[births_de$season == "summer"]/prev_up_95,
                                 intervention == "Nirsevimab for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_up_95))


# ------ Plot and save outputs -----------------------------------------------
# VE_distribution %>% ggplot(aes(x = t, y = VE_t, col = group)) +
#   geom_line()

palette <- c("No immunisation" = "#9C964A", "Maternal vaccination" = "#6A9D96",
             "Nirsevimab for spring births" = "#88BBA0",
             "Nirsevimab for summer births" = "#85D4E3",
             "Nirsevimab for spring and summer births" = "#B39BC8")

# One iteration
total_hosp_intervention_df %>% filter (iter == 3) %>% ggplot() +
  geom_col(aes(x = intervention, y = n_hospitalisations, fill = intervention), position="dodge") +
  labs(x = "\nIntervention",
       y = "Number of RSV-associated hospitalisations in children under the age of 1 year\n") +
  theme_light() +
  theme (#axis.text.y=element_blank(), 
    axis.ticks.y=element_blank(),
    legend.position = "None",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    title = element_text(size = 20)) +
  scale_fill_manual(values = palette)

# All with CI
plt <- total_hosp_intervention_int %>% ggplot() +
  geom_col(aes(x = intervention, y = hosp_median, fill = intervention), position = "dodge") +
  geom_errorbar(aes(x = intervention, ymin = hosp_low_95, ymax = hosp_up_95), 
                width = 0.2, position = "dodge") +
  labs(x = "\nIntervention",
       y = "Number of RSV-associated hospitalisations in children under the age of 1 year\n") +
  theme_light() +
  theme (axis.ticks.y=element_blank(),
    legend.position = "None",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    title = element_text(size = 20)) +
  scale_fill_manual(values = palette)

plt

plt_nnv <- nnv %>% filter (intervention != "No immunisation") %>% 
  ggplot() +
  geom_col(aes(x = intervention, y = NNV_median, fill = intervention), position = "dodge") +
  geom_errorbar(aes(x = intervention, ymin = NNV_low_95, ymax = NNV_up_95), 
                width = 0.2, position = "dodge") +
  labs(x = "\nIntervention",
       y = "NNV to prevent one hospitalisation\n") +
  theme_light() +
  theme (axis.ticks.y=element_blank(),
         legend.position = "None",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         title = element_text(size = 20)) +
  scale_fill_manual(values = palette)

plt_nnv

# Save files
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"
write.csv(total_hosp_intervention_int, paste0(path, "Final outputs.csv"), row.names = F)
write.csv(nnv, paste0(path, "Final outputs NNV hosp.csv"), row.names = F)

plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Hosp by intervention CI.png",
           width = 14, height = 16, units = "in", 
           device='png')
plt_nnv %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/NNV hosp by intervention CI.png",
               width = 14, height = 16, units = "in", 
               device='png')
