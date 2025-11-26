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

births_de <- births_de %>%
  rbind(
    births_de %>%
      summarise(season = "all",
                total = sum(total),
                prop = sum(prop),
                N = N[1])
  )

# Read in RSV-illness estimates
mild_illness <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Mild illness numbers") %>% 
  janitor::clean_names()
severe_illness <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Severe illness numbers") %>%
  janitor::clean_names()

severe_illness_new <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Severe illness rates") %>% 
  janitor::clean_names()

mild_illness_new <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Mild illness rates") %>% 
  janitor::clean_names()

# ----------- Data -------------------------------------------------------------
# Read in and the data and put it in the right format
data <- #read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/refs/heads/master/data/infection_status.csv")
read.csv("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Original code/RSV_infection_NL/Data/infection_status.csv")
# Group age into intervals 
# bi-monthly for 0-2 years and 6-monthly for 2-5 years
data$age_grp <- cut(data$age_days,
                    breaks = c(seq(0,730, by = 30.25*2),
                               seq(909,2000, by = 30.25*6)), 
                    include.lowest = T, right = F)

# Divide by season of birth
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)

data <- data %>%
  mutate(
    Birth_mo = format(as.Date(Birth_doy, origin = "2020-12-31"), "%m") %>% as.numeric(),
    season_birth = case_when (Birth_mo %in% spring ~ "spring",
                              Birth_mo %in% summer ~ "summer",
                              Birth_mo %in% autumn ~ "autumn",
                              Birth_mo %in% winter ~ "winter"))

get_midpoint <- function(cut_label) {
  round(mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ",")))))
}

data$xMidpoint <- sapply(data$age_grp, get_midpoint)

# Different groupings
# Get number of cases by age
incidence_data <- data %>% select (age_grp, age_days, infection) %>%
  mutate(N_tot = n()) %>%
  group_by(age_grp) %>%
  reframe(time = round(median(age_days)), 
          N = n(),
          n_infection = sum(infection),
          prop_seroconv = n_infection/N) %>%
  ungroup() %>% 
  distinct() %>%
  mutate(cum_pop = cumsum(N),
         incidence = n_infection/cum_pop)

# Calculate seroprevalence and binomial confidence intervals
incidence_data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom::binom.confint(incidence_data$n_infection, 
                                                                                             incidence_data$N,
                                                                                             method="exact")[,c("mean","lower","upper")]

# Calculate incidence and binomial confidence intervals
incidence_data[,c("incidence_mean","incidence_low95","incidence_up95")] <- binom::binom.confint(incidence_data$n_infection, 
                                                                                                incidence_data$cum_pop,
                                                                                                method="exact")[,c("mean","lower","upper")]


incidence_data_season <- data %>% select (age_grp, age_days, infection, season_birth, xMidpoint) %>%
  mutate(N_tot = n()) %>%
  group_by(age_grp, xMidpoint, season_birth) %>%
  summarise(age_mid = round(median(age_days)), 
            N = n(),
            n_infection = sum(infection),
            prop_seroconv = n_infection/N) %>%
  ungroup() %>% 
  distinct() %>%
  group_by(season_birth) %>%
  mutate(cum_pop = cumsum(N),
         incidence = n_infection/cum_pop)

incidence_data_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom::binom.confint(incidence_data_season$n_infection, 
                                                                                                    incidence_data_season$N, 
                                                                                                    method="exact")[,c("mean","lower","upper")]

# Calculate incidence and binomial confidence intervals
incidence_data_season[,c("incidence_mean","incidence_low95","incidence_up95")] <- binom::binom.confint(incidence_data_season$n_infection, 
                                                                                                       incidence_data_season$cum_pop,
                                                                                                       method="exact")[,c("mean","lower","upper")]


incidence_data_season_wide <- incidence_data_season %>% 
  select (!c(age_mid, age_grp, seroprev_mean, incidence_mean, cum_pop)) %>%
  pivot_wider(
    names_from = season_birth,
    values_from = c(N , n_infection, prop_seroconv, seroprev_low95, seroprev_up95, 
                    incidence, incidence_low95,incidence_up95),
    values_fill = 0
  ) %>%
  rename(time = "xMidpoint")

# --------- Define our functions ----------------------------------------------
# Newly seroconverted at a given age
new_seroconv_age <- function (index){
  converted_all_rand <- as.numeric(t(converted_all[index,0:366]))
  converted_sp_rand <- as.numeric(t(converted_sp[index,0:366]))
  converted_sm_rand <- as.numeric(t(converted_sm[index,0:366]))
  converted_au_rand <- as.numeric(t(converted_au[index,0:366]))
  converted_wt_rand <- as.numeric(t(converted_wt[index,0:366]))
  
  converted <- data.frame(age_midpoint = 0:365,
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
illness_type <- mild_illness_new %>%
  merge(severe_illness_new, by = "age_months") %>%
  group_by(age_months) %>% # total rate of cases in each category first
  mutate(total_rate = total_mild_illness_rate + total_severe_illness_rate,
            total_ma_rate = ma_mild_illness_rate + ma_severe_illness_rate,
            total_non_ma_rate = non_ma_mild_illness_rate + non_ma_severe_illness_rate,

            # prop of each type of mild cases
            prop_mild_rate = total_mild_illness_rate / total_rate,
            prop_ma_mild_rate = ma_mild_illness_rate / total_rate,
            prop_non_ma_mild_rate = non_ma_mild_illness_rate / total_rate,

            # prop of each type of severe cases
            prop_severe_rate = total_severe_illness_rate / total_rate,
            prop_ma_severe_rate = ma_severe_illness_rate / total_rate,
            prop_non_ma_severe_rate = non_ma_severe_illness_rate / total_rate,

            # prop of cases by MA status
            prop_ma_rate = total_ma_rate / total_rate,
            prop_non_ma_rate = total_non_ma_rate / total_rate) %>%
  ungroup()

# # Get proportion of cases that are hospitalised
# prop_hosp <- severe_illness %>%
#   merge(mild_illness, by = "age_months") %>%
#   group_by(age_months) %>%
#   mutate(total_cases = total_severe_illness_number + total_mild_illness_number) %>%
#   group_by(age_months) %>%
#   summarise(prop_hosp = total_severe_illness_number / total_cases)

# ----- Data processing
# Number of infections
infections <- function(conversion, births){
  # Convert ages to the same units
  conversion_formated <- conversion %>% mutate(age_months = trunc(age_midpoint/30)) %>% # turn age into months
    arrange(age_midpoint) # arrange by age instead of by season of birth
  
  # Combine seroconversion data with birth data
  infections <- conversion_formated %>%
    merge(births %>% select(total, season),
          by.x = "season_birth" , by.y = "season")
  
  # Get number of infections
  infections <- infections %>%
    mutate(n_infections = incidence * total)
  
  return (infections)
}

# Get number of hospitalisations
hosp <- function (infections){
  hospitalisations <- infections %>%
    merge(prop_hosp %>% select(age_months, prop_hosp),
          by = "age_months") %>%
    mutate(n_severe = n_infections * prop_hosp,
           n_severe = case_when(n_severe < 0 ~ 0,
                                T ~ n_severe))
  return (hospitalisations)
}


s_progression_format <- function (infections, severe_illness) {
  # Put ages in the same categories as the disease progression data
  infections_formated <- infections %>% 
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
  
  # severe_illness_progression <- conversion_formated %>%
  #   merge(illness_type, by.x = "age_bracket", by.y = "age_months") %>%
  #   mutate(severe_cases = case_when((total_severe_illness_rate/100000) * incidence < 0 ~ 0,
  #                            is.infinite((total_severe_illness_rate/100000) * incidence ) ~ NA,
  #                            T ~ (total_severe_illness_rate/100000) * incidence ),
  #          mild_cases = case_when((total_mild_illness_rate/100000) * incidence < 0 ~ 0,
  #                                   is.infinite((total_mild_illness_rate/100000) * incidence ) ~ NA,
  #                                   T ~ (total_mild_illness_rate/100000) * incidence ),
  #          ma_severe_cases = case_when((ma_severe_illness_rate/100000) * incidence < 0 ~ 0,
  #                                   is.infinite((ma_severe_illness_rate/100000) * incidence ) ~ NA,
  #                                   T ~ (ma_severe_illness_rate/100000) * incidence ),
  #          non_ma_severe_cases = case_when((non_ma_severe_illness_rate/100000) * incidence < 0 ~ 0,
  #                                   is.infinite((non_ma_severe_illness_rate/100000) * incidence ) ~ NA,
  #                                   T ~ (non_ma_severe_illness_rate/100000) * incidence )) %>%
  #   select(age_midpoint, age_months, age_bracket, season_birth, 
  #          severe_cases, mild_cases, ma_severe_cases, non_ma_severe_cases)
  
  # Get proportion of severe cases in SA
  severe_illness_progression <- infections_formated %>%
    filter(season_birth != "all") %>%
    group_by(age_bracket) %>%
    summarise(n_infections = sum(n_infections)) %>%
    ungroup() %>%
    merge(severe_illness, by.x = "age_bracket", by.y = "age_months") %>%
    group_by(age_bracket) %>%
    summarise(incidence_hosp = case_when(is.infinite(total_severe_illness_number/n_infections) ~ 0,
                                      T ~ total_severe_illness_number/n_infections)) %>%
    filter(age_bracket != "12-14") # remove this age group as we don't have it in our data
  
 
  # severe_illness_progression %>%
  #   mutate(age_bracket = factor(age_bracket,
  #                               levels = c("<1", "1", "2", "3", "4", "5",
  #                                          "6", "7", "8", "9", "10", "11"))) %>%
  #   group_by(age_bracket) %>%
  #   summarise(incidence_hosp = sum(incidence_hosp)) %>%
  #   ggplot(aes(x = age_bracket, y = incidence_hosp)) +
  #   geom_bar(stat = "identity")
  
  # Apply that to the birth seasons
  severe_illness_progression <- severe_illness_progression %>%
    merge(infections_formated, by = c("age_bracket")) %>%
    mutate(n_hosp = n_infections * incidence_hosp) %>%
    select(age_midpoint, age_months, age_bracket, season_birth, n_hosp)
    
  # severe_illness_progression %>%
  #   mutate(age_bracket = factor(age_bracket,
  #                               levels = c("<1", "1", "2", "3", "4", "5",
  #                                          "6", "7", "8", "9", "10", "11"))) %>%
  #   group_by(season_birth, age_bracket) %>%
  #   summarise(n_hosp = sum(n_hosp)) %>%
  #   ggplot(aes(x = age_bracket, y = n_hosp, fill = season_birth)) +
  #   geom_bar(position = position_dodge(), stat = "identity")
  
  return (severe_illness_progression)
}

m_progression_format <- function (conversion, illness_type) {
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
    mutate(mild_cases = case_when((total_mild_illness_rate/100000) * incidence < 0 ~ 0,
                             is.infinite((total_mild_illness_rate/100000) * incidence ) ~ NA,
                             T ~ (total_mild_illness_rate/100000) * incidence ),
           ma_mild_cases = case_when((ma_mild_illness_rate/100000) * incidence < 0 ~ 0,
                             is.infinite((ma_mild_illness_rate/100000) * incidence ) ~ NA,
                             T ~ (ma_mild_illness_rate/100000) * incidence ),
           non_ma_mild_cases = case_when((non_ma_mild_illness_rate/100000) * incidence < 0 ~ 0,
                             is.infinite((non_ma_mild_illness_rate/100000) * incidence ) ~ NA,
                             T ~ (non_ma_mild_illness_rate/100000) * incidence )) %>%
    select(age_midpoint, age_months, age_bracket, season_birth, 
           mild_cases, ma_mild_cases, non_ma_mild_cases)
  
  return (mild_illness_progression)
}

scale_cases <- function (severe_illness) {
  # Get totals
  severe_cases_u1_de <- severe_illness 
  
  # Scale the numbers to have them match Fabienne's paper
  # About 22,000 hospitalisations in 2019, of which about 12% in year 1 --> 2,640 
  # But we had 10,564 in 2024 according to InEKDatenBrowser for 28 days - 1 year
  # Model predicts about 79,000 hospitalisations --> scale by 0.033 for Fabienne and by 0.133 for InEKDatenBrowser
  severe_cases_u1_de <- severe_cases_u1_de %>% 
    mutate (n_hosp_scaled = n_hosp * 0.133) # 0.04)
  
  return (severe_cases_u1_de)
}

m_cases_format <- function (mild_illness) {
  # Restrict to < 1 year
  mild_illness_u1 <- mild_illness %>% filter (age_months <= 12)
  
  # Add pop size in <1 y.o by season
  mild_cases_u1_de <- mild_illness_u1 %>%
    merge(births_de %>% select(total, season),
          by.x = "season_birth" , by.y = "season")
  
  # Get number of cases
  mild_cases_u1_de <- mild_cases_u1_de %>%
    # mutate(n_mild = total_mild_cases_prop * total,
    #        n_ma_mild = total_ma_mild_cases_prop * total,
    #        n_non_ma_mild = total_non_ma_mild_cases_prop * total)
    mutate(n_mild = mild_cases * total,
           total_mild_cases_prop = NA,
           total_ma_mild_cases_prop = NA,
           total_non_ma_mild_cases_prop = NA,
           n_ma_mild = ma_mild_cases * total,
           n_non_ma_mild = non_ma_mild_cases * total)

  # Get totals
  mild_cases_u1_de <- mild_cases_u1_de %>% rbind(
    mild_cases_u1_de %>% summarise (season_birth = "all",
                                    age_bracket = "all",
                                    age_months = NA,
                                    age_midpoint = NA,
                                    mild_cases = NA,
                                    ma_mild_cases = NA,
                                    non_ma_mild_cases = NA,
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
  # But we had 10,564 in 2024 according to InEKDatenBrowser for 28 days - 1 year
  # Model predicts about 65,000 hospitalisations --> scale by 0.04 for Fabienne and by 0.163 for InEKDatenBrowser
  mild_cases_u1_de <- mild_cases_u1_de %>% 
    mutate (n_mild_scaled = n_mild * 0.088, # 0.04,
            n_ma_mild_scaled = n_ma_mild * 0.088, #0.04,
            n_non_ma_mild_scaled = n_non_ma_mild * 0.088) #0.04)
  
  return (mild_cases_u1_de)
}

# We need to combine mild and severe MA cases
count_ma <- function (mild_cases, severe_cases){
  ma_cases_c <- mild_cases %>% filter(!is.na(age_months)) %>%
    select (season_birth, age_months, age_midpoint, age_bracket, n_ma_mild_scaled) %>%
    mutate(severity = "mild") %>%
    rename(n_ma_cases = "n_ma_mild_scaled") %>%
    rbind(
      severe_cases %>% filter(!is.na(age_months)) %>%
        select (season_birth, age_months, age_midpoint, age_bracket, n_ma_severe_scaled) %>%
        mutate(severity = "severe") %>%
        rename(n_ma_cases = "n_ma_severe_scaled")
    )
  
  ma_cases <- ma_cases_c %>%
    rbind(
      ma_cases_c %>% group_by(season_birth, age_months, age_midpoint, age_bracket) %>%
        summarise (n_ma_cases = sum(n_ma_cases),
                   severity = "all")
    ) %>%
    rbind(
      ma_cases_c %>% group_by(age_months, age_midpoint, age_bracket) %>%
        summarise (n_ma_cases = sum(n_ma_cases),
                   severity = "all",
                   season_birth = "all")
    ) %>%
    rbind(
      ma_cases_c %>% group_by(severity, age_months, age_midpoint, age_bracket) %>%
        summarise (n_ma_cases = sum(n_ma_cases),
                   season_birth = "all")
    )
  
  return (ma_cases)
}

# ----- Interventions
# Effect of maternal vaccination on hospitalisations
hosp_vacc <- function (cases, VE_distribution){
  # Add the estimated VE to the number of cases
  hosp_prevented_vacc <- cases %>% filter(!is.na(age_months)) %>%
    select(season_birth, age_bracket, age_months, age_midpoint, n_hosp_scaled) %>%
    merge(
      VE_distribution %>% filter(group == "severe") %>% select(t, VE_t),
      by.x = "age_midpoint", by.y = "t")
  
  # Get number of prevented cases
  hosp_prevented_vacc <- hosp_prevented_vacc %>% 
    group_by(season_birth, age_bracket, age_months, age_midpoint, VE_t) %>%
    # summarise(n_ma_severe_vacc = n_ma_severe_scaled * (1-VE_t), # number of cases despite vaccination
    #           n_ma_severe_averted = n_ma_severe_scaled * VE_t)
    summarise(n_vacc = n_hosp_scaled * (1-VE_t), # number of cases despite vaccination
              n_averted = n_hosp_scaled * VE_t)
  
  hosp_prevented_vacc <- hosp_prevented_vacc %>% 
    filter(season_birth != "all") %>%
    rbind(
      hosp_prevented_vacc %>% 
        filter(season_birth != "all") %>%
        group_by(season_birth) %>%
        summarise(age_bracket = 'all', 
                  age_months = NA,
                  age_midpoint = NA, 
                  n_vacc = sum(n_vacc),
                  n_averted = sum(n_averted))
    )
  return (hosp_prevented_vacc)
}

# Effect of maternal vaccination on MA cases
ma_vacc <- function (ma_cases, VE_distribution){
  # Add the estimated VE
  ma_prevented_vacc <- ma_cases %>% filter(!is.na(age_months)) %>%
    select(season_birth, age_months, age_midpoint, n_ma_cases, severity) %>%
    merge(
      VE_distribution %>%
        select(t, VE_t, group) %>%
        mutate(group = case_when (group == "LRTI" ~ "mild",
                                  TRUE ~ group)),
      by.x = c("age_midpoint","severity"), by.y = c("t", "group")
    )
  
  ma_prevented_vacc <- ma_prevented_vacc %>% 
    group_by(season_birth, age_midpoint, age_months, severity, VE_t) %>%
    summarise(n_ma_cases_vacc = n_ma_cases * (1-VE_t), # number of cases despite vaccination
              n_ma_cases_averted = n_ma_cases * VE_t)
  
  ma_prevented_vacc_ <- ma_prevented_vacc %>% 
    mutate (age_bracket = as.character(age_months)) %>%
    rbind(
      ma_prevented_vacc %>% 
        group_by(season_birth) %>%
        summarise(age_months = NA,
                  age_midpoint = NA,
                  age_bracket = "all",
                  severity = "all",
                  n_ma_cases_vacc = sum(n_ma_cases_vacc, na.rm = T),
                  n_ma_cases_averted = sum(n_ma_cases_averted, na.rm = T))
    ) %>%
    rbind(
      ma_prevented_vacc %>% 
        group_by(season_birth, severity) %>%
        summarise(age_months = NA,
                  age_midpoint = NA,
                  age_bracket = "all",
                  n_ma_cases_vacc = sum(n_ma_cases_vacc, na.rm = T),
                  n_ma_cases_averted = sum(n_ma_cases_averted, na.rm = T))
    ) %>%
    rbind(
      ma_prevented_vacc %>% 
        group_by(season_birth, age_midpoint, age_months) %>%
        summarise(age_bracket = as.character(age_months),
                  severity = "all",
                  n_ma_cases_vacc = sum(n_ma_cases_vacc, na.rm = T),
                  n_ma_cases_averted = sum(n_ma_cases_averted, na.rm = T)) %>%
        unique()
    ) %>% ungroup()
  
  return (ma_prevented_vacc_)
}

# Effect of mAB on hospitalisations
hosp_mAB <- function (hosp_prevented_vacc, IE_distribution){
  hosp_prevented <- hosp_prevented_vacc %>% filter(!is.na(age_months)) %>%
    ungroup() %>%
    select(season_birth, age_bracket, age_months, age_midpoint,
           n_vacc, n_averted) %>%
    merge(
      IE_distribution %>% 
        filter(group == "hospitalisation") %>%
        select(t, IE_t),
      by.x = "age_midpoint", by.y = "t")
  
  # Get number of prevented cases
  hosp_prevented_mAB <- hosp_prevented %>%
    group_by(season_birth, age_bracket, age_months, age_midpoint, IE_t) %>%
    mutate(n_nirs = n_vacc * (1-IE_t),
            n_averted_nirs = n_vacc * IE_t)
  
  # Get number of cases averted or not by giving nirvesimab at a given age
  hosp_prevented_age <- hosp_prevented_mAB %>% group_by(season_birth) %>%
    # averted hospitalisations
    summarise (prev_u1 = sum(n_averted_nirs),
               prev_1 = sum(n_averted_nirs[age_months >= 1], na.rm = T),
               prev_3 = sum(n_averted_nirs[age_months >= 3], na.rm = T),
               prev_4 = sum(n_averted_nirs[age_months >= 4], na.rm = T),
               prev_7 = sum(n_averted_nirs[age_months >= 7], na.rm = T),
               prev_8 = sum(n_averted_nirs[age_months >= 8], na.rm = T),
               prev_10 = sum(n_averted_nirs[age_months >= 10], na.rm = T),
               prev_11 = sum(n_averted_nirs[age_months >= 11], na.rm = T),
               # prevented hospitalisations when immunising each birth cohort just before their first winter
               prev_spring_8 = sum(n_averted_nirs[age_months >= 8 & season_birth == "spring"], na.rm = T),
               prev_summer_4 = sum(n_averted_nirs[age_months >= 4 & season_birth == "summer"], na.rm = T),
               prev_autumn_1 = sum(n_averted_nirs[age_months >= 1 & season_birth == "autumn"], na.rm = T),
               prev_winter_11 = sum(n_averted_nirs[age_months >= 11 & season_birth == "winter"], na.rm = T),
               # remaining hospitalisations
               hosp_u1 = sum(n_nirs),
               hosp_1 = sum(n_nirs[age_months >= 1], na.rm = T) + sum(n_vacc[age_months < 1], na.rm = T),
               hosp_3 = sum(n_nirs[age_months >= 3], na.rm = T) + sum(n_vacc[age_months < 3], na.rm = T),
               hosp_4 = sum(n_nirs[age_months >= 4], na.rm = T) + sum(n_vacc[age_months < 4], na.rm = T),
               hosp_7 = sum(n_nirs[age_months >= 7], na.rm = T) + sum(n_vacc[age_months < 7], na.rm = T),
               hosp_8 = sum(n_nirs[age_months >= 8], na.rm = T) + sum(n_vacc[age_months < 8], na.rm = T),
               hosp_10 = sum(n_nirs[age_months >= 10], na.rm = T) + sum(n_vacc[age_months < 10], na.rm = T),
               hosp_11 = sum(n_nirs[age_months >= 11], na.rm = T) + sum(n_vacc[age_months < 11], na.rm = T),
    ) %>%
    ungroup() %>%
    # hospitalisations when immunising one cohort only (hosp. are averted in the cohort only)
    mutate(hosp_spring_8 = case_when (season_birth == "spring" ~ hosp_8,
                                      TRUE ~ prev_8 + hosp_8),
           hosp_summer_4 = case_when (season_birth == "summer" ~ hosp_4,
                                      TRUE ~ prev_4 + hosp_4),
           hosp_autumn_1 = case_when (season_birth == "autumn" ~ hosp_1,
                                      TRUE ~ prev_1 + hosp_1),
           hosp_winter_11 = case_when (season_birth == "winter" ~ hosp_11,
                                       TRUE ~ prev_11 + hosp_11)) 
  
  hosp_prevented_age <- hosp_prevented_age %>% 
    filter(season_birth != "all") %>%
    rbind(
      data.frame(season_birth = "all",
                 # averted hospitalisations
                 prev_u1 = sum(hosp_prevented_mAB$n_averted_nirs),
                 prev_1 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 1], na.rm = T),
                 prev_3 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 3], na.rm = T),
                 prev_4 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 4], na.rm = T),
                 prev_7 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 7], na.rm = T),
                 prev_8 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 8], na.rm = T),
                 prev_10 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 10], na.rm = T),
                 prev_11 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 11], na.rm = T),
                 # prevented hospitalisations when immunising one cohort only
                 prev_spring_8 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 8 & hosp_prevented_mAB$season_birth == "spring"], na.rm = T),
                 prev_summer_4 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 4 & hosp_prevented_mAB$season_birth == "summer"], na.rm = T),
                 prev_autumn_1 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 1 & hosp_prevented_mAB$season_birth == "autumn"], na.rm = T),
                 prev_winter_11 = sum(hosp_prevented_mAB$n_averted_nirs[hosp_prevented_mAB$age_months >= 11 & hosp_prevented_mAB$season_birth == "winter"], na.rm = T),
                 # remaining hospitalisations
                 hosp_u1 = sum(hosp_prevented_mAB$n_nirs),
                 hosp_1 = sum(hosp_prevented_age$hosp_1, na.rm = T),
                 hosp_3 = sum(hosp_prevented_age$hosp_3, na.rm = T),
                 hosp_4 = sum(hosp_prevented_age$hosp_4, na.rm = T),
                 hosp_7 = sum(hosp_prevented_age$hosp_7, na.rm = T),
                 hosp_8 = sum(hosp_prevented_age$hosp_8, na.rm = T),
                 hosp_10 = sum(hosp_prevented_age$hosp_10, na.rm = T),
                 hosp_11 = sum(hosp_prevented_age$hosp_11, na.rm = T),
                 hosp_spring_8 = 0,
                 hosp_summer_4 = 0,
                 hosp_autumn_1 = 0,
                 hosp_winter_11 = 0) %>%# placeholder
        ungroup()
    ) %>%
    mutate(hosp_spring_8 = case_when (season_birth == "all" ~ sum(hosp_spring_8),
                                      TRUE ~ hosp_spring_8),
           hosp_summer_4 = case_when (season_birth == "all" ~ sum(hosp_summer_4),
                                      TRUE ~ hosp_summer_4),
           hosp_autumn_1 = case_when (season_birth == "all" ~ sum(hosp_autumn_1),
                                      TRUE ~ hosp_autumn_1),
           hosp_winter_11 = case_when (season_birth == "all" ~ sum(hosp_winter_11),
                                       TRUE ~ hosp_winter_11)) %>%
    ungroup()
  
  return (hosp_prevented_age)
}

# Effect of mAB on MA case
ma_mAB <- function (ma_prevented_vacc, IE_distribution) {
  # Get number of prevented cases
  ma_prevented_nirs <- ma_prevented_vacc %>% filter(!is.na(age_months)) %>%
    select(season_birth, age_midpoint, age_months, n_ma_cases_vacc, severity) %>%
    merge(
      IE_distribution %>% filter(group == "medically-attended") %>%
        select(t, IE_t),
      by.x = "age_midpoint", by.y = "t")
  
  ma_prevented_nirs <- ma_prevented_nirs %>% 
    group_by(season_birth, age_midpoint, severity, IE_t) %>%
    mutate(n_ma_cases_nirs = n_ma_cases_vacc * (1-IE_t), # number of cases despite vaccination
              n_ma_cases_averted_nirs = n_ma_cases_vacc * IE_t) %>%
    ungroup()
  
  # Get number of MA cases averted by given nirvesimab at a given age
  ma_prevented_age_nirs <- ma_prevented_nirs %>% group_by(season_birth, severity) %>%
    summarise (prev_u1 = sum(n_ma_cases_averted_nirs),
               prev_1 = sum(n_ma_cases_averted_nirs[age_months >= 1], na.rm = T),
               prev_3 = sum(n_ma_cases_averted_nirs[age_months >= 3], na.rm = T),
               prev_4 = sum(n_ma_cases_averted_nirs[age_months >= 4], na.rm = T),
               prev_7 = sum(n_ma_cases_averted_nirs[age_months >= 7], na.rm = T),
               prev_8 = sum(n_ma_cases_averted_nirs[age_months >= 8], na.rm = T),
               prev_10 = sum(n_ma_cases_averted_nirs[age_months >= 10], na.rm = T),
               prev_11 = sum(n_ma_cases_averted_nirs[age_months >= 11], na.rm = T),
               # prevented cases when immunising one cohort only
               prev_spring_8 = sum(n_ma_cases_averted_nirs[age_months >= 8 & season_birth == "spring"], na.rm = T),
               prev_summer_4 = sum(n_ma_cases_averted_nirs[age_months >= 4 & season_birth == "summer"], na.rm = T),
               prev_autumn_1 = sum(n_ma_cases_averted_nirs[age_months >= 1 & season_birth == "autumn"], na.rm = T),
               prev_winter_11 = sum(n_ma_cases_averted_nirs[age_months >= 11 & season_birth == "winter"], na.rm = T),
               # remaining cases
               ma_u1 = sum(n_ma_cases_nirs),
               ma_1 = sum(n_ma_cases_nirs[age_months >= 1], na.rm = T) + sum(n_ma_cases_vacc[age_months < 1], na.rm = T),
               ma_3 = sum(n_ma_cases_nirs[age_months >= 3], na.rm = T) + sum(n_ma_cases_vacc[age_months < 3], na.rm = T),
               ma_4 = sum(n_ma_cases_nirs[age_months >= 4], na.rm = T) + sum(n_ma_cases_vacc[age_months < 4], na.rm = T),
               ma_7 = sum(n_ma_cases_nirs[age_months >= 7], na.rm = T) + sum(n_ma_cases_vacc[age_months < 7], na.rm = T),
               ma_8 = sum(n_ma_cases_nirs[age_months >= 8], na.rm = T) + sum(n_ma_cases_vacc[age_months < 8], na.rm = T),
               ma_10 = sum(n_ma_cases_nirs[age_months >= 10], na.rm = T) + sum(n_ma_cases_vacc[age_months < 10], na.rm = T),
               ma_11 = sum(n_ma_cases_nirs[age_months >= 11], na.rm = T) + sum(n_ma_cases_vacc[age_months < 11], na.rm = T)) %>%
    ungroup() %>%
    # MA cases when immunising one cohort only (cases are averted in the cohort only)
    mutate(ma_spring_8 = case_when (season_birth == "spring"  ~ ma_8,
                                     season_birth == "all"  ~ 0,
                                     TRUE ~ prev_8 + ma_8),
           ma_summer_4 = case_when (season_birth == "summer" ~ ma_4,
                                    season_birth == "all"  ~ 0,
                                    TRUE ~ prev_4 + ma_4),
           ma_autumn_1 = case_when (season_birth == "autumn" ~ ma_1,
                                    season_birth == "all"  ~ 0,
                                    TRUE ~ prev_1 + ma_1),
           ma_winter_11 = case_when (season_birth == "winter" ~ ma_11,
                                    season_birth == "all"  ~ 0,
                                    TRUE ~ prev_11 + ma_11)) %>%
        ungroup() %>%
    group_by(severity) %>%
    mutate(ma_spring_8 = case_when (season_birth == "all" ~ sum(ma_spring_8),
                                       TRUE ~ ma_spring_8),
           ma_summer_4 = case_when (season_birth == "all" ~ sum(ma_summer_4),
                                      TRUE ~ ma_summer_4),
           ma_autumn_1 = case_when (season_birth == "all" ~ sum(ma_autumn_1),
                                    TRUE ~ ma_autumn_1),
           ma_winter_11 = case_when (season_birth == "all" ~ sum(ma_winter_11),
                                    TRUE ~ ma_winter_11),
           prev_spring_8 = case_when (season_birth == "all" ~ prev_spring_8[season_birth == "spring"],
                                     TRUE ~ prev_spring_8),
           prev_summer_4 = case_when (season_birth == "all" ~ prev_summer_4[season_birth == "summer"],
                                    TRUE ~ prev_summer_4),
           prev_autumn_1 = case_when (season_birth == "all" ~ prev_autumn_1[season_birth == "autumn"],
                                      TRUE ~ prev_autumn_1),
           prev_winter_11 = case_when (season_birth == "all" ~ prev_winter_11[season_birth == "winter"],
                                      TRUE ~ prev_winter_11)
    ) %>%
    ungroup()
  
  return(ma_prevented_age_nirs)
}

# Summarising the outputs
hosp_intervention <- function (cases, hosp_prevented_vacc, hosp_prevented_age){
  # Get the numbers
  total_hosp_intervention <- cases %>%
    # select(season_birth, age_bracket, n_ma_severe_scaled) %>%
    # rename(n_hospitalisations = n_ma_severe_scaled) %>%
    select(season_birth, age_bracket, n_hosp_scaled) %>%
    rename(n_hospitalisations = n_hosp_scaled) %>%
    mutate(prevented_hospitalisations = 0,
           intervention = "No immunisation") %>%
    rbind(
      hosp_prevented_vacc %>% filter(age_bracket == "all") %>%
        ungroup() %>%
        # select(season_birth, age_bracket, n_ma_severe_vacc, n_ma_severe_averted) %>%
        select(season_birth, age_bracket, n_vacc, n_averted) %>%
        group_by(season_birth) %>%
        summarise(age_bracket = "all",
                   n_hospitalisations = sum(n_vacc),
                   prevented_hospitalisations = sum(n_averted),
                   intervention = "MV")
    ) %>%
    rbind(
      hosp_prevented_vacc %>% 
        filter(age_bracket != "all") %>%
        ungroup() %>%
        select(season_birth, age_bracket, n_vacc, n_averted) %>%
        group_by(season_birth, age_bracket) %>%
        summarise(n_hospitalisations = sum(n_vacc),
                  prevented_hospitalisations = sum(n_averted),
                  intervention = "MV")
    ) %>%
    rbind(
      hosp_prevented_vacc %>% filter(age_bracket == "all") %>%
        ungroup() %>%
        # select(age_bracket, n_ma_severe_vacc, n_ma_severe_averted) %>%
        select(age_bracket, n_vacc, n_averted) %>%
        summarise(season_birth = "all",
                  age_bracket = "all",
                  n_hospitalisations = sum(n_vacc),
                  prevented_hospitalisations = sum(n_averted),
                  intervention = "MV")
    ) %>%
    rbind(
      hosp_prevented_age %>%
        select(hosp_spring_8, prev_spring_8, season_birth) %>%
        group_by(season_birth) %>%
        summarise(age_bracket = "all",
                  n_hospitalisations = hosp_spring_8,
                  prevented_hospitalisations = prev_spring_8,
                  intervention = "+ mAB for spring births")
    ) %>%
    rbind(
      hosp_prevented_age %>%
        select(season_birth, hosp_summer_4, prev_summer_4) %>%
        group_by(season_birth) %>%
        summarise(age_bracket = "all",
                  n_hospitalisations = hosp_summer_4,
                  prevented_hospitalisations = prev_summer_4,
                  intervention = "+ mAB for summer births")
    ) %>%
    rbind(
      hosp_prevented_age %>% 
        select(season_birth, hosp_autumn_1, prev_autumn_1) %>%
        group_by(season_birth) %>%
        summarise(age_bracket = "all",
                  n_hospitalisations = hosp_autumn_1,
                  prevented_hospitalisations = prev_autumn_1,
                  intervention = "+ mAB for autumn births")
    ) %>%
    rbind(
      hosp_prevented_age %>% 
        select(season_birth, hosp_winter_11, prev_winter_11) %>%
        group_by(season_birth) %>%
        summarise(age_bracket = "all",
                  n_hospitalisations = hosp_winter_11,
                  prevented_hospitalisations = prev_winter_11,
                  intervention = "+ mAB for winter births")
    ) %>%
    rbind(
      hosp_prevented_age %>% 
        select(season_birth, hosp_summer_4, prev_summer_4, 
               hosp_spring_8, prev_spring_8) %>%
        group_by(season_birth) %>%
        summarise(age_bracket = "all",
                  n_hospitalisations = 0, # we need to calculate this later
                  prevented_hospitalisations = prev_summer_4 + prev_spring_8,
                  intervention = "+ mAB for spring and summer births")
    ) %>% mutate(intervention = factor(intervention, 
                                       levels = c("No immunisation", "MV", 
                                                  "+ mAB for spring births",
                                                  "+ mAB for summer births",
                                                  "+ mAB for autumn births",
                                                  "+ mAB for winter births",
                                                  "+ mAB for spring and summer births")))
  
  total_hosp_intervention <- total_hosp_intervention %>%
    group_by(season_birth) %>%
    mutate(n_hospitalisations = case_when (intervention == "+ mAB for spring and summer births" ~
                                             n_hospitalisations[intervention == "MV" & age_bracket == "all"] - prevented_hospitalisations[intervention == "+ mAB for spring and summer births"],
                                           TRUE ~ n_hospitalisations))
  
  return (total_hosp_intervention)
}

ma_intervention <- function (cases, ma_prevented_vacc, ma_prevented_age){
  # Get the numbers
  total_ma_intervention <- cases %>%
    filter(severity == "all") %>%
    select(season_birth, n_ma_cases, severity) %>%
    rename(n_cases = n_ma_cases) %>%
    group_by(season_birth) %>% 
    summarise(#season_birth = "all",
              age_bracket = "all",
              n_cases = sum(n_cases),
              prevented_cases = 0,
              severity = "all",
              intervention = "No immunisation") %>%
    rbind(
      ma_prevented_vacc %>% filter(age_bracket == "all" & severity == "all") %>%
        ungroup() %>%
        select(season_birth, severity, n_ma_cases_vacc, n_ma_cases_averted) %>%
        group_by(season_birth) %>%
        summarise(#season_birth = "all",
                  age_bracket = "all",
                  n_cases = sum(n_ma_cases_vacc),
                  prevented_cases = sum(n_ma_cases_averted),
                  severity = "all",
                  intervention = "MV")
    ) %>%
    rbind(
      ma_prevented_age %>% filter(severity == "all") %>%
        ungroup() %>%
        select(season_birth, ma_spring_8, prev_spring_8) %>%
        group_by(season_birth) %>%
        summarise(#season_birth = "all",
                  age_bracket = "all",
                  n_cases = ma_spring_8,
                  prevented_cases = prev_spring_8,
                  severity = "all",
                  intervention = "+ mAB for spring births")
    ) %>%
    rbind(
      ma_prevented_age %>% filter(severity == "all") %>%
        ungroup() %>%
        select(ma_summer_4, prev_summer_4, season_birth) %>%
        group_by(season_birth) %>%
        summarise(#season_birth = "all",
                  age_bracket = "all",
                  n_cases = ma_summer_4,
                  prevented_cases = prev_summer_4,
                  severity = "all",
                  intervention = "+ mAB for summer births")
    ) %>%
    rbind(
      ma_prevented_age %>% filter(severity == "all") %>%
        ungroup() %>%
        select(, season_birth,ma_autumn_1, prev_autumn_1) %>%
        group_by(season_birth) %>%
        summarise(#season_birth = "all",
                  age_bracket = "all",
                  n_cases = ma_autumn_1,
                  prevented_cases = prev_autumn_1,
                  severity = "all",
                  intervention = "+ mAB for autumn births")
    ) %>%
    rbind(
      ma_prevented_age %>% filter(severity == "all") %>%
        ungroup() %>%
        select(season_birth, ma_winter_11, prev_winter_11) %>%
        group_by(season_birth) %>%
        summarise(#season_birth = "all",
                  age_bracket = "all",
                  n_cases = ma_winter_11,
                  prevented_cases = prev_winter_11,
                  severity = "all",
                  intervention = "+ mAB for winter births")
    ) %>%
    rbind(
      ma_prevented_age %>% filter(severity == "all") %>%
        ungroup() %>%
        select(season_birth, ma_summer_4, prev_summer_4, ma_spring_8, prev_spring_8) %>%
        group_by(season_birth) %>%
        summarise(#season_birth = "all",
                  age_bracket = "all",
                  n_cases = 0, # we need to calculate this later
                  prevented_cases = prev_summer_4 + prev_spring_8,
                  severity = "all",
                  intervention = "+ mAB for spring and summer births")
    ) %>% mutate(intervention = factor(intervention, 
                                       levels = c("No immunisation", "MV", 
                                                  "+ mAB for spring births",
                                                  "+ mAB for summer births",
                                                  "+ mAB for autumn births",
                                                  "+ mAB for winter births",
                                                  "+ mAB for spring and summer births")))
  
  total_ma_intervention <- total_ma_intervention %>%
    group_by(season_birth) %>%
    mutate(n_cases = case_when (intervention == "+ mAB for spring and summer births" ~
                                             n_cases[intervention == "MV"] - prevented_cases[intervention == "+ mAB for spring and summer births"],
                                           TRUE ~ n_cases))
  
  return (total_ma_intervention)
}
# -----------------------------------------------------------------------------

# --------- Run the models once -----------------------------------------------
# Get VE
source(paste0(path_code, "VE estimates.R")) # takes about 10 min 
# Get IE of mAB
source(paste0(path_code, "mAB IE estimates.R")) # takes about X min
# Get seroconversion by age
# source(paste0(path_code, "Main.R")) # takes about 3 hours 
converted_summary <- readRDS("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/RDS files/seroconversion simulation outputs.rds")
converted_all <- do.call(rbind.data.frame, converted_summary[1])
converted_sp <- do.call(rbind.data.frame, converted_summary[2])
converted_sm <- do.call(rbind.data.frame, converted_summary[3])
converted_au <- do.call(rbind.data.frame, converted_summary[4])
converted_wt <- do.call(rbind.data.frame, converted_summary[5])

# Replace negative values with 0
converted_all[converted_all<0] <- 0
converted_sp[converted_sp<0] <- 0
converted_sm[converted_sm<0] <- 0
converted_au[converted_au<0] <- 0
converted_wt[converted_wt<0] <- 0


# -----------------------------------------------------------------------------

# --------- Combine the estimates ---------------------------------------------
# Define where we're going to store the outputs
total_hosp_intervention <- list()
# total_ma_intervention <- list()

for (i in 1:100){
  # ------ Calculate number of RSV cases by age in Germany ------------
  # Proportion of newly seroconverted children at a given age
  # Take one random iteration of the model
  set.seed(i)
  rand <- base::sample(1:10000, 1)
  # And extract the values for it
  conversion <- new_seroconv_age (rand)
  
  # Number of infections
  n_infections <- infections(conversion, births_de)

  # Proportion of children with RSV at a given age (seroconversion to disease)
  severe_illness_gen <- s_progression_format (n_infections, severe_illness)
  # mild_illness <- m_progression_format (conversion, illness_type)
  # Number of cases in Germany
  severe_illness_de <- scale_cases (severe_illness_gen)
  # mild_illness_de <- m_cases_format(mild_illness)
  
  # ------ Maternal vaccination -----------------------------------------------
  # Take one random iteration of the model
  set.seed(42*i)
  draw <- base::sample(1:10000, 1)
  # And extract the values for it
  VE_distribution <- model_ve_runs_t %>% filter(iter == draw)
  
  # Calculate number of hospitalisations despite the vaccine
  hosp_prevented_vacc <- hosp_vacc(severe_illness_de, VE_distribution)
  
  # Calculate number of MA cases despite the vaccine
  # ma_cases <- count_ma(mild_illness_de, severe_illness_de)
  # ma_prevented_vacc <- ma_vacc(ma_cases, VE_distribution)
  
  # ------ Monoclonal antibody -----------------------------------------------
  # Take one random iteration of the model
  draw <- base::sample(1:10000, 1)
  # And extract the values for it
  IE_distribution <- model_mab_runs_t %>% filter(iter == draw)
  # Calculate number of hospitalisations despite the vaccine and the mAB
  hosp_prevented_age <- hosp_mAB(hosp_prevented_vacc, IE_distribution)
  
  # Calculate number of MA cases despite the vaccine and the mAB
  # ma_prevented_age <- ma_mAB(ma_prevented_vacc, IE_distribution)
  
  # Get total hospitalisations by intervention
  total_hosp_intervention[[i]] <- hosp_intervention(severe_illness_de, hosp_prevented_vacc, hosp_prevented_age)
  total_hosp_intervention[[i]]$iter <- rep(i, each = 1878) # save index in case we want to check one iteration
  
  # Get total MA cases by intervention
  # total_ma_intervention[[i]] <- ma_intervention(ma_cases, ma_prevented_vacc, ma_prevented_age)
  # total_ma_intervention[[i]]$iter <- rep(i, each = 35) # save index in case we want to check one iteration
}
  
total_hosp_intervention_df <- do.call("rbind", total_hosp_intervention)
# total_ma_intervention_df <- do.call("rbind", total_ma_intervention)

total_hosp_intervention_df <- total_hosp_intervention_df %>%
  rbind(
    total_hosp_intervention_df %>% 
      filter (intervention == "No immunisation") %>%
      group_by(intervention, season_birth, iter) %>%
      summarise(age_bracket = "all",
                n_hospitalisations = sum(n_hospitalisations),
                prevented_hospitalisations = 0) %>%
      ungroup()
  )

incidence_hosp <- total_hosp_intervention_df %>%
  merge(births_de %>% select(season, total), by.x = "season_birth", by.y = "season") %>%
  mutate(incidence_per_100000 = n_hospitalisations/total*100000)

# Get reduction in incidence compared to no immunisation
incidence_reduction_no_imm <- incidence_hosp %>%
  filter (intervention != "No immunisation" & age_bracket == "all") %>%
  merge(
    incidence_hosp %>%
      filter (intervention == "No immunisation" & age_bracket == "all") %>%
      select(season_birth, age_bracket, n_hospitalisations, iter) %>%
      rename(n_hospitalisations_no_int = n_hospitalisations),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(incidence_reduction_perc = (n_hospitalisations_no_int - n_hospitalisations)/n_hospitalisations_no_int*100)

# Get reduction in incidence for mAB compared to MV
incidence_reduction_mab_mv <- incidence_hosp %>%
  filter (intervention != "No immunisation" & intervention != "MV" & age_bracket == "all") %>%
  merge(
    incidence_hosp %>%
      filter (intervention == "MV" & age_bracket == "all") %>%
      select(season_birth, age_bracket, n_hospitalisations, iter) %>%
      rename(n_hospitalisations_mv = n_hospitalisations),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(incidence_reduction_perc = (n_hospitalisations_mv - n_hospitalisations)/n_hospitalisations_mv*100)

# Compute 95% CI
total_hosp_intervention_int <- total_hosp_intervention_df %>% 
  filter (age_bracket == "all") %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
            hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
            hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T),
            prev_low_95 = quantile(prevented_hospitalisations, 0.025, na.rm = T),
            prev_median = quantile(prevented_hospitalisations, 0.5, na.rm = T),
            prev_up_95 = quantile(prevented_hospitalisations, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Winter", "Spring", "Summer", "Autumn",
                                          "All")))

incidence_hosp_int <- incidence_hosp %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(incidence_up_95 = quantile(incidence_per_100000, 0.025, na.rm = T),
            incidence_median = quantile(incidence_per_100000, 0.5, na.rm = T),
            incidence_low_95 = quantile(incidence_per_100000, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Winter", "Spring", "Summer", "Autumn",
                                          "All")))

incidence_reduction_mab_mv_int <- incidence_reduction_mab_mv %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(incidence_reduction_low_95 = quantile(incidence_reduction_perc, 0.025, na.rm = T),
            incidence_reduction_median = quantile(incidence_reduction_perc, 0.5, na.rm = T),
            incidence_reduction_up_95 = quantile(incidence_reduction_perc, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Winter", "Spring", "Summer", "Autumn",
                                          "All")))


incidence_reduction_no_imm_int <- incidence_reduction_no_imm %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(incidence_reduction_low_95 = quantile(incidence_reduction_perc, 0.025, na.rm = T),
            incidence_reduction_median = quantile(incidence_reduction_perc, 0.5, na.rm = T),
            incidence_reduction_up_95 = quantile(incidence_reduction_perc, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Winter", "Spring", "Summer", "Autumn",
                                          "All")))

# Get absolute reduction in incidence for mAB compared to MV
abs_inc_reduction <- incidence_hosp %>%
  filter (intervention == "MV" & age_bracket == "all" & season_birth == "spring") %>%
  select(season_birth, age_bracket, incidence_per_100000, iter) %>%
  merge(
    incidence_hosp %>%
      filter (intervention == "+ mAB for spring births") %>%
      select(season_birth, age_bracket, incidence_per_100000, iter) %>%
      rename(incidence_mAB = incidence_per_100000),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
         season_birth = "spring") %>%
  rbind(
    incidence_hosp %>%
      filter (intervention == "MV" & age_bracket == "all" & season_birth == "summer") %>%
      select(season_birth, age_bracket, incidence_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter (intervention == "+ mAB for summer births") %>%
          select(season_birth, age_bracket, incidence_per_100000, iter) %>%
          rename(incidence_mAB = incidence_per_100000),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
             season_birth = "summer")
  ) %>%
  rbind(
    incidence_hosp %>%
      filter (intervention == "MV" & age_bracket == "all" & season_birth == "autumn") %>%
      select(season_birth, age_bracket, incidence_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter (intervention == "+ mAB for autumn births") %>%
          select(season_birth, age_bracket, incidence_per_100000, iter) %>%
          rename(incidence_mAB = incidence_per_100000),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
             season_birth = "autumn")
  ) %>%
  rbind(
    incidence_hosp %>%
      filter (intervention == "MV" & age_bracket == "all" & season_birth == "winter") %>%
      select(season_birth, age_bracket, incidence_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter (intervention == "+ mAB for winter births") %>%
          select(season_birth, age_bracket, incidence_per_100000, iter) %>%
          rename(incidence_mAB = incidence_per_100000),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
             season_birth = "winter")
  )

abs_inc_reduction_ci <- abs_inc_reduction %>%
  group_by(season_birth, age_bracket) %>%
  summarise(abs_incidence_reduction_low_95 = quantile(abs_incidence_reduction, 0.025, na.rm = T),
            abs_incidence_reduction_median = quantile(abs_incidence_reduction, 0.5, na.rm = T),
            abs_incidence_reduction_up_95 = quantile(abs_incidence_reduction, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Winter", "Spring", "Summer", "Autumn",
                                          "All")))

# total_ma_intervention_int <- total_ma_intervention_df %>% 
#   group_by(intervention, season_birth) %>%
#   summarise(cases_low_95 = quantile(n_cases, 0.025, na.rm = T),
#             cases_median = quantile(n_cases, 0.5, na.rm = T),
#             cases_up_95 = quantile(n_cases, 0.975, na.rm = T),
#             prev_low_95 = quantile(prevented_cases, 0.025, na.rm = T),
#             prev_median = quantile(prevented_cases, 0.5, na.rm = T),
#             prev_up_95 = quantile(prevented_cases, 0.975, na.rm = T)) %>%
#   ungroup() %>%
#   mutate(season_birth = str_to_sentence(season_birth)) %>%
#   mutate(season_birth = factor(season_birth, 
#                                levels = c("Winter", "Spring", "Summer", "Autumn",
#                                           "All")))

# Get NNV to prevent one hospitalisation
nnv <- total_hosp_intervention_int %>%
  select(intervention, season_birth, prev_low_95, prev_median, prev_up_95) %>%
  mutate(NNV_low_95 = case_when (intervention == "No immunisation" ~ NA,
                                 intervention == "MV" ~ births_de$N[1]/prev_low_95,
                                 intervention == "+ mAB for spring births" ~ births_de$total[births_de$season == "spring"]/prev_low_95,
                                 intervention == "+ mAB for summer births" ~ births_de$total[births_de$season == "summer"]/prev_low_95,
                                 intervention == "+ mAB for autumn births" ~ births_de$total[births_de$season == "autumn"]/prev_low_95,
                                 intervention == "+ mAB for winter births" ~ births_de$total[births_de$season == "winter"]/prev_low_95,
                                 intervention == "+ mAB for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_low_95),
         NNV_median = case_when (intervention == "No immunisation" ~ NA,
                                 intervention == "MV" ~ births_de$N[1]/prev_median,
                                 intervention == "+ mAB for spring births" ~ births_de$total[births_de$season == "spring"]/prev_median,
                                 intervention == "+ mAB for summer births" ~ births_de$total[births_de$season == "summer"]/prev_median,
                                 intervention == "+ mAB for autumn births" ~ births_de$total[births_de$season == "autumn"]/prev_median,
                                 intervention == "+ mAB for winter births" ~ births_de$total[births_de$season == "winter"]/prev_median,
                                 intervention == "+ mAB for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_median),
         NNV_up_95 = case_when (intervention == "No immunisation" ~ NA,
                                 intervention == "MV" ~ births_de$N[1]/prev_up_95,
                                 intervention == "+ mAB for spring births" ~ births_de$total[births_de$season == "spring"]/prev_up_95,
                                 intervention == "+ mAB for summer births" ~ births_de$total[births_de$season == "summer"]/prev_up_95,
                                intervention == "+ mAB for autumn births" ~ births_de$total[births_de$season == "autumn"]/prev_up_95,
                                intervention == "+ mAB for winter births" ~ births_de$total[births_de$season == "winter"]/prev_up_95,
                                 intervention == "+ mAB for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_up_95))


# Get NNV to prevent one MA cases
# nnv_ma <- total_ma_intervention_int %>% select(intervention, season_birth, 
#                                                prev_low_95, prev_median, prev_up_95) %>%
#   mutate(NNV_low_95 = case_when (intervention == "No immunisation" ~ NA,
#                                  intervention == "MV" ~ births_de$N[1]/prev_low_95,
#                                  intervention == "+ mAB for spring births" ~ births_de$total[births_de$season == "spring"]/prev_low_95,
#                                  intervention == "+ mAB for summer births" ~ births_de$total[births_de$season == "summer"]/prev_low_95,
#                                  intervention == "+ mAB for autumn births" ~ births_de$total[births_de$season == "autumn"]/prev_low_95,
#                                  intervention == "+ mAB for winter births" ~ births_de$total[births_de$season == "winter"]/prev_low_95,
#                                  intervention == "+ mAB for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_low_95),
#          NNV_median = case_when (intervention == "No immunisation" ~ NA,
#                                  intervention == "MV" ~ births_de$N[1]/prev_median,
#                                  intervention == "+ mAB for spring births" ~ births_de$total[births_de$season == "spring"]/prev_median,
#                                  intervention == "+ mAB for summer births" ~ births_de$total[births_de$season == "summer"]/prev_median,
#                                  intervention == "+ mAB for autumn births" ~ births_de$total[births_de$season == "autumn"]/prev_median,
#                                  intervention == "+ mAB for winter births" ~ births_de$total[births_de$season == "winter"]/prev_median,
#                                  intervention == "+ mAB for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_median),
#          NNV_up_95 = case_when (intervention == "No immunisation" ~ NA,
#                                 intervention == "MV" ~ births_de$N[1]/prev_up_95,
#                                 intervention == "+ mAB for spring births" ~ births_de$total[births_de$season == "spring"]/prev_up_95,
#                                 intervention == "+ mAB for summer births" ~ births_de$total[births_de$season == "summer"]/prev_up_95,
#                                 intervention == "+ mAB for autumn births" ~ births_de$total[births_de$season == "autumn"]/prev_up_95,
#                                 intervention == "+ mAB for winter births" ~ births_de$total[births_de$season == "winter"]/prev_up_95,
#                                 intervention == "+ mAB for spring and summer births" ~ (births_de$total[births_de$season == "spring"] + births_de$total[births_de$season == "summer"])/prev_up_95))
# 

# ------ Plot and save outputs -----------------------------------------------
# VE_distribution %>% ggplot(aes(x = t, y = VE_t, col = group)) +
#   geom_line()

palette <- c("No immunisation" = "#9C964A", "MV" = "#6A9D96",
             "+ mAB for spring births" = "#88BBA0",
             "+ mAB for summer births" = "#85D4E3",
             "+ mAB for autumn births" = "#13D4C7",
             "+ mAB for winter births" = "#B479E0",
             "+ mAB for spring and summer births" = "#B39BC8")

palette_season <- c("Autumn" = "#88BBA0", "Winter" = "#13D4C7", 
                    "Spring" = "#B39BC8", "Summer" = "#B479E0",
                    "All" = "#9C964A")

# One iteration
total_hosp_intervention_df %>% 
  filter (iter == 1 & season_birth != "all") %>% 
  ggplot() +
  geom_col(aes(x = intervention, y = n_hospitalisations, fill = season_birth), 
           position="stack") +
  labs(x = "\nIntervention",
       y = "Number of RSV-associated hospitalisations in children under the age of 1 year\n") +
  theme_light() +
  theme ( axis.ticks.y = element_blank(),
          legend.position = "None",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25)) +
  # scale_fill_manual(values = palette_season) +
  scale_y_continuous(labels=scales::comma)

# total_ma_intervention_df %>% filter (iter == 3 & season_birth != "all") %>% 
#   ggplot() +
#   geom_col(aes(x = intervention, y = n_cases, fill = season_birth), position="stack") +
#   labs(x = "\nIntervention",
#        y = "Number of RSV-associated MA cases in children under the age of 1 year\n") +
#   theme_light() +
#   theme (#axis.text.y=element_blank(), 
#     axis.ticks.y=element_blank(),
#     legend.position = "right",
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
#     axis.title.x = element_text(size = 20),
#     axis.title.y = element_text(size = 20),
#     title = element_text(size = 20)) +
#   scale_fill_manual(values = palette_season)

# All with CI
# We need to recalculate the CI for No immunisation
total_hosp_intervention_plt <- total_hosp_intervention_df %>%
  filter(!(intervention == "No immunisation" & season_birth == "all")) %>%
  rbind(
    total_hosp_intervention_df %>% 
      filter (intervention == "No immunisation" & season_birth != "all") %>%
      group_by(intervention,season_birth, iter) %>%
      summarise(age_bracket = "all",
                n_hospitalisations = sum(n_hospitalisations),
                prevented_hospitalisations = 0) %>%
      ungroup()
  )

total_hosp_intervention_plt %>%
  filter(intervention == 'No immunisation') %>%
  group_by(intervention, age_bracket, iter) %>%
  summarise(season_birth = "all",
            total_hosp = sum(n_hospitalisations)) %>%
  ungroup()

# Compute 95% CI
total_hosp_intervention_plt_int <- total_hosp_intervention_plt %>% 
  filter (age_bracket == "all") %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
            hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
            hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T),
            prev_low_95 = quantile(prevented_hospitalisations, 0.025, na.rm = T),
            prev_median = quantile(prevented_hospitalisations, 0.5, na.rm = T),
            prev_up_95 = quantile(prevented_hospitalisations, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Winter", "Spring", "Summer", "Autumn",
                                          "All")))


plt <- total_hosp_intervention_int %>% 
  filter (season_birth != "All") %>%
  ggplot() +
  geom_col(aes(x = intervention, y = hosp_median, fill = season_birth), position = "stack") +
  geom_errorbar(data = subset(total_hosp_intervention_int, 
                              season_birth == "All"),
                aes(x = intervention, ymin = hosp_low_95, ymax = hosp_up_95), 
                width = 0.2, position = "dodge") +
  labs(x = "\nIntervention",
       y = "Number of RSV-associated hospitalisations in\nchildren under the age of 1 year\n",
       fill = "Season of birth") +
  theme_light() +
  theme (axis.ticks.y=element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20)) +
  scale_fill_manual(values = palette_season) +
  scale_y_continuous(labels=scales::comma)

plt

# plt_ma <- total_ma_intervention_int %>% filter (season_birth != "all") %>% 
#   ggplot() +
#   geom_col(aes(x = intervention, y = cases_median, fill = season_birth), position = "stack") +
#   geom_errorbar(data = subset(total_ma_intervention_int, season_birth == "all"),
#                 aes(x = intervention, ymin = cases_low_95, ymax = cases_up_95), 
#                 width = 0.2, position = "dodge") +
#   labs(x = "\nIntervention",
#        y = "Number of RSV-associated MA cases in children under the age of 1 year\n",
#        fill = "Season of birth") +
#   theme_light() +
#   theme (axis.ticks.y=element_blank(),
#          legend.position = "right",
#          axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
#          axis.text.y = element_text(size = 18),
#          axis.title.x = element_text(size = 20),
#          axis.title.y = element_text(size = 20),
#          title = element_text(size = 20),
#          legend.text = element_text(size = 18)) +
#   scale_fill_manual(values = palette_season)
# 
# plt_ma # We have more MA cases than hospitalisations because we take severe = hospitalisation
# as opposed to MA severe = hospitalisation

plt_nnv <- nnv %>% filter (intervention != "No immunisation" & season_birth == "All") %>% 
  ggplot() +
  geom_col(aes(x = intervention, y = NNV_median, fill = intervention), position = "dodge") +
  geom_errorbar(aes(x = intervention, ymin = NNV_low_95, ymax = NNV_up_95), 
                width = 0.2, position = "dodge") +
  labs(x = "\nIntervention",
       y = "NNI to prevent one hospitalisation\n") +
  theme_light() +
  theme (axis.ticks.y=element_blank(),
         legend.position = "None",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         axis.text.y = element_text(size = 20)) +
    scale_fill_manual(values = palette) +
  scale_y_continuous(labels=scales::comma)

plt_nnv

# plt_nnv_cases <- nnv_ma %>% filter (intervention != "No immunisation" & season_birth == "all") %>% 
#   ggplot() +
#   geom_col(aes(x = intervention, y = NNV_median, fill = intervention), 
#            position = "stack") +
#   geom_errorbar(aes(x = intervention, ymin = NNV_low_95, ymax = NNV_up_95), 
#                 width = 0.2, position = "dodge") +
#   labs(x = "\nIntervention",
#        y = "NNV to prevent one MA case\n") +
#   theme_light() +
#   theme (axis.ticks.y=element_blank(),
#          legend.position = "none",
#          axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
#          axis.title.x = element_text(size = 20),
#          axis.title.y = element_text(size = 20),
#          title = element_text(size = 20)) +
#   scale_fill_manual(values = palette)
# 
# plt_nnv_cases

# Plot outcomes by age
total_hosp_no_int_age <- total_hosp_intervention_df %>%
  filter(age_bracket != "all" & intervention == 'No immunisation') %>%
  group_by(age_bracket, iter, intervention) %>%
  summarise(season_birth = "all",
         n_hospitalisations = sum(n_hospitalisations)) %>%
  rbind(total_hosp_intervention_df %>%
          filter(age_bracket != "all" & intervention == 'No immunisation') %>%
  group_by(age_bracket, iter, season_birth, intervention) %>%
  summarise(n_hospitalisations = sum(n_hospitalisations, na.rm = T))
  ) %>%
  ungroup() %>%
  group_by(age_bracket, intervention, season_birth) %>%
  summarise(hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
            hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
            hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("winter", "spring", "summer", "autumn",
                                          "all")),
         age_bracket = factor(age_bracket, 
                               levels = c("<1", "1", "2", "3", "4", "5", "6", "7", "8", 
                               "9", "10", "11", "12-14")))

plt_hosp <- total_hosp_no_int_age %>% 
  ggplot() +
            geom_point(aes(x = age_bracket, y = hosp_median, col = season_birth),
                       size = 3) +
            geom_errorbar(aes(x = age_bracket, ymin = hosp_low_95, ymax = hosp_up_95,
                              col = season_birth), width = 0.5, linewidth = 1) +
            labs(x = "\nAge (months)",
                 y = "Hospitalisations\n",
                 title = "Total hospitalisations without intervention",
                 col = "Season of birth") +
            theme_light() +
            facet_wrap(~season_birth) +
  theme (axis.ticks.y=element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
         axis.text.y = element_text(size = 25),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 20, color = "black")) 
plt_hosp

total_hosp_mv_age <- total_hosp_intervention_df %>%
  filter(age_bracket != "all" & intervention == 'MV') %>%
  group_by(age_bracket, iter, intervention) %>%
  summarise(season_birth = "all",
            n_hospitalisations = sum(n_hospitalisations)) %>%
  rbind(total_hosp_intervention_df %>%
          filter(age_bracket != "all" & intervention == 'MV') %>%
  group_by(age_bracket, iter, season_birth, intervention) %>%
  summarise(n_hospitalisations = sum(n_hospitalisations, na.rm = T))
  ) %>%
  ungroup() %>%
  group_by(age_bracket, intervention, season_birth) %>%
  summarise(hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
            hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
            hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("winter", "spring", "summer", "autumn",
                                          "all")),
         age_bracket = factor(age_bracket, 
                              levels = c("<1", "1", "2", "3", "4", "5", "6", "7", "8", 
                                         "9", "10", "11", "12-14")))

plt_hosp_mv <- total_hosp_mv_age %>% 
  ggplot() +
  geom_point(aes(x = age_bracket, y = hosp_median, col = season_birth),
             size = 3) +
  geom_errorbar(aes(x = age_bracket, ymin = hosp_low_95, ymax = hosp_up_95,
                    col = season_birth), width = 0.5, linewidth = 1) +
  labs(x = "\nAge (months)",
       y = "Hospitalisations\n",
       title = "Total hospitalisations with MV",
       col = "Season of birth") +
  theme_light() +
  facet_wrap(~season_birth) +
  theme (axis.ticks.y=element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
         axis.text.y = element_text(size = 25),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 20, color = "black")) 
plt_hosp_mv

# Look at total burden by season of birth
total_hosp_1_y <- total_hosp_intervention_df %>%
  filter(age_bracket == "all") %>%
  group_by(season_birth, intervention) %>%
  summarise(hosp_low_95 = quantile(n_hospitalisations, 0.025),
            hosp_median = quantile(n_hospitalisations, 0.5),
            hosp_up_95 = quantile(n_hospitalisations, 0.975)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth))

total_hosp_1_y_plt <- total_hosp_1_y %>% ggplot() +
  geom_col(aes(x = intervention, y = hosp_median, fill = season_birth), 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = intervention, ymin = hosp_low_95, ymax = hosp_up_95, group = season_birth), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "Number of RSV-associated hospitalisations in children under the age of 1 year",
       x = "\nIntervention",
       y = "Hospitalisations\n",
       fill = "Season of birth") +
  theme_light() +
  theme (axis.ticks.y=element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         legend.text = element_text(size = 25)) +
  scale_fill_manual(values = palette_season)

total_hosp_1_y_plt

# Save files
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"
write.csv(total_hosp_intervention_int, paste0(path, "Final outputs hospitalisations.csv"), row.names = F)
write.csv(nnv, paste0(path, "Final outputs NNV hosp.csv"), row.names = F)
write.csv(incidence_hosp_int, paste0(path, "Final outputs hosp incidence.csv"), row.names = F)
write.csv(incidence_reduction_int, paste0(path, "Final outputs incidence vs no.csv"), row.names = F)
write.csv(incidence_reduction_mab_mv_int, paste0(path, "Final outputs incidence vs MV.csv"), row.names = F)
write.csv(abs_inc_reduction_ci, paste0(path, "Final outputs absolute incidence reduction.csv"), row.names = F)
write.csv(total_ma_intervention_int, paste0(path, "Final outputs MA cases.csv"), row.names = F)
write.csv(nnv_ma, paste0(path, "Final outputs NNV MA cases.csv"), row.names = F)

plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Hosp by intervention CI.png",
           width = 14, height = 16, units = "in", 
           device='png')

plt_ma %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Cases by intervention CI.png",
                  width = 14, height = 16, units = "in", 
                   device='png')

plt_nnv %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/NNV hosp by intervention CI.png",
               width = 14, height = 16, units = "in", 
               device='png')

plt_nnv_cases %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/NNV cases by intervention CI.png",
                   width = 14, height = 16, units = "in", 
                   device='png')

plt_hosp %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Hospitalisations no intervention.png",
                         width = 14, height = 16, units = "in", 
                         device='png')

total_hosp_1_y_plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Total hospitalisations.png",
                              width = 14, height = 16, units = "in", 
                              device='png')

# Plot seroconversions with CI
converted_ci <- list()

for (index in seq(1:1000)){
  converted_all_rand <- as.numeric(t(converted_all[index,0:366]))
  converted_sp_rand <- as.numeric(t(converted_sp[index,0:366]))
  converted_sm_rand <- as.numeric(t(converted_sm[index,0:366]))
  converted_au_rand <- as.numeric(t(converted_au[index,0:366]))
  converted_wt_rand <- as.numeric(t(converted_wt[index,0:366]))
  
  converted <- data.frame(age_midpoint = 0:365,
                          All = converted_all_rand,
                          Spring = converted_sp_rand,
                          Summer = converted_sm_rand,
                          Autumn = converted_au_rand,
                          Winter = converted_wt_rand)
  
  converted_ci[[index]] <- converted
  converted_ci[[index]]$iter <- index # save index in case we want to check one iteration
}

converted_ci_df <- do.call("rbind", converted_ci)
converted_ci_long <- converted_ci_df %>%
  pivot_longer(cols = !c(age_midpoint, iter),
               names_to = "season_birth",
               values_to = "seroconverted") %>%
  mutate(season_birth = factor(season_birth, levels = c("Autumn", "Winter",
                                                         "Spring", "Summer",
                                                         "All")))

converted_ci <- converted_ci_long %>% 
  group_by(season_birth, age_midpoint) %>%
  summarise(sero_low_95 = quantile(seroconverted, 0.025, na.rm = T),
            sero_median = quantile(seroconverted, 0.5, na.rm = T),
            sero_up_95 = quantile(seroconverted, 0.975, na.rm = T)) %>%
  ungroup()

converted_plt <- converted_ci %>%
  ggplot() +
  geom_line(aes(x = age_midpoint, y = sero_median, col = season_birth),
            linewidth = 1.2) +
  geom_ribbon(aes(x = age_midpoint, ymin = sero_low_95, ymax = sero_up_95, 
                    fill = season_birth),
              alpha = 0.4) +
  facet_wrap (~season_birth) +
  labs(title = "Proportion of seroconverted children under the age of 1 year",
       x = "\nAge (days)",
       y = " Proportion seroconverted\n",
       colour = "Season of birth",
       fill = "Season of birth") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 25),
         axis.text.y = element_text(size = 18),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         legend.text = element_text(size = 25),
         strip.text.x = element_text(size = 20, color = "black")) +
  scale_y_continuous(labels = scales::percent_format())

converted_plt

converted_plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Prop seroconverted.png",
                              width = 18, height = 16, units = "in", 
                              device='png')

# Plot VE with CI
VE_ci <- model_ve_runs_t %>%
  filter(group == "severe") %>%
  group_by(t) %>%
  summarise(ve_low_95 = quantile(VE_t, 0.05),
            ve_median = quantile(VE_t, 0.5),
            ve_up_95 = quantile(VE_t, 0.95)) %>%
  ungroup()

VE_ci_plt <- VE_ci %>%
  ggplot() +
  geom_line(aes(x = t, y = ve_median),
            linewidth = 1.2) +
  geom_ribbon(aes(x = t, ymin = ve_low_95, ymax = ve_up_95),
              alpha = 0.4) +
  labs(title = "Efficacy of the maternal vaccine",
       x = "\nAge (days)",
       y = "VE\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 25),
         axis.text.y = element_text(size = 18),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         legend.text = element_text(size = 25)) +
  scale_y_continuous(labels = scales::percent_format())

VE_ci_plt 

VE_ci_plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/VE waning.png",
                         width = 18, height = 16, units = "in", 
                         device='png')

VE_inv_plt <- VE_ci %>%
  filter(t/30.25 <= 11) %>%
  ggplot() +
  geom_line(aes(x = t/30.25, y = 1-ve_median),
            linewidth = 1.5, colour = "blue4") +
  geom_ribbon(aes(x = t/30.25, ymax = 1-ve_low_95, ymin = 1-ve_up_95),
              alpha = 0.4, fill = "darkslategray1") +
  labs(x = "\nAge (months)",
       y = "Relative risk of hospitalisation despite receiving the MV\n(1-VE)\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
         axis.title.x = element_text(size = 30),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 30),
         title = element_blank()) +
  scale_x_continuous(breaks = seq(0,12,1))

VE_inv_plt 

VE_inv_plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/OR VE.png",
                     width = 20, height = 12, units = "in", 
                     device='png')

# Do seroconversion * (1-VE) 100 times and plot it
odds_MV <- list()

for (i in 1:100){
  # Take one random iteration of the model
  rand <- base::sample(1:10000, 1)
  # And extract the values for it
  converted_all_rand <- as.numeric(t(converted_all[rand,0:366]))
  converted_sp_rand <- as.numeric(t(converted_sp[rand,0:366]))
  converted_sm_rand <- as.numeric(t(converted_sm[rand,0:366]))
  converted_au_rand <- as.numeric(t(converted_au[rand,0:366]))
  converted_wt_rand <- as.numeric(t(converted_wt[rand,0:366]))
  
  converted <- data.frame(age_midpoint = 0:365,
                          All = converted_all_rand,
                          Spring = converted_sp_rand,
                          Summer = converted_sm_rand,
                          Autumn = converted_au_rand,
                          Winter = converted_wt_rand)

  
  # Maternal vaccination
  # Take one random iteration of the model
  draw <- base::sample(1:10000, 1)
  # And extract the values for it
  VE_distribution <- model_ve_runs_t %>% filter(iter == draw)
  # Get OR
  OR_dist <- VE_distribution %>%
    mutate(OR_t = 1 - VE_t) %>%
    select(t, OR_t)
  
  # Combined the two
  risk <- converted %>%
    left_join(OR_dist, by = c("age_midpoint" = "t")) %>%
    mutate(adjusted_seroconverted_all = All * OR_t,
           adjusted_seroconverted_sp = Spring * OR_t,
           adjusted_seroconverted_sm = Summer * OR_t,
           adjusted_seroconverted_au = Autumn * OR_t,
           adjusted_seroconverted_wt = Winter * OR_t)
  
  odds_MV[[i]] <- risk
  odds_MV[[i]]$iter <- i # save index in case we want to check one iteration
}

odds_MV_df <- do.call("rbind", odds_MV)
odds_MV_ci <- odds_MV_df %>%
  pivot_longer(cols = c(adjusted_seroconverted_all, adjusted_seroconverted_sp,
                        adjusted_seroconverted_sm, adjusted_seroconverted_au,
                        adjusted_seroconverted_wt),
               names_to = "season_birth",
               values_to = "adjusted_seroconverted") %>%
  mutate(season_birth = case_when(season_birth == "adjusted_seroconverted_all" ~ "All",
                                 season_birth == "adjusted_seroconverted_sp" ~ "Spring",
                                 season_birth == "adjusted_seroconverted_sm" ~ "Summer",
                                 season_birth == "adjusted_seroconverted_au" ~ "Autumn",
                                 season_birth == "adjusted_seroconverted_wt" ~ "Winter")) %>%
  group_by(season_birth, age_midpoint) %>%
  summarise(adj_sero_low_95 = quantile(adjusted_seroconverted, 0.025, na.rm = T),
            adj_sero_median = quantile(adjusted_seroconverted, 0.5, na.rm = T),
            adj_sero_up_95 = quantile(adjusted_seroconverted, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = factor(season_birth, levels = c("Autumn", "Winter",
                                                         "Spring", "Summer",
                                                         "All")))

odds_MV_plt <- odds_MV_ci %>%
  ggplot() +
  geom_line(aes(x = age_midpoint, y = adj_sero_median, colour = season_birth),
            linewidth = 1.2) +
  geom_ribbon(aes(x = age_midpoint, ymin = adj_sero_low_95, ymax = adj_sero_up_95, 
                    fill = season_birth),
              alpha = 0.4) +
  facet_wrap (~season_birth) +
  labs(x = "\nAge (days)",
       y = "Risk of hospitalisation despite receiving the MV\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
         axis.title.x = element_text(size = 30),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 30),
         strip.text.x = element_text(size = 30, color = "black")) 

odds_MV_plt 

odds_MV_plt  %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/OR with risk and VE.png",
                      width = 18, height = 16, units = "in", 
                      device='png')

# Age dependent disease progression
severe_illness_plt <- severe_illness_new %>%
  filter(age_months != "12-14" & age_months != "15-17" & age_months != "18-20" & 
           age_months != "21-23" & age_months != "24-35" & age_months != "36-47" &
           age_months != "48-59" & age_months != "<5 years" & age_months !="12" &
           age_months != "< 1 year") %>%
  mutate(age_months = factor(age_months, levels = c("<1", "1", "2", "3", "4", 
                                                    "5", "6", "7", "8", "9", "10",
                                                    "11"))) %>%
  ggplot() +
  geom_point(aes(x = age_months, y = total_severe_illness_rate), size = 2) +
  geom_errorbar(aes(x = age_months, ymin = total_severe_illness_lower_ci, 
                    ymax = total_severe_illness_upper_ci)) +
 labs(x = "\nAge (months)",
       y = "Severe RSV illness rates under the age of 1 year\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
         axis.title.x = element_text(size = 30),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 30)) +
  scale_y_continuous(labels = scales::comma)

severe_illness_plt

severe_illness_plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Severe illness rates.png",
                              width = 20, height = 12, units = "in", 
                              device='png')

# Get average age at hospitalisation
av_age_hosp <- total_hosp_intervention_df %>%
  ungroup() %>%
  filter(age_bracket != "all" & age_bracket != "12-14") %>%
  mutate (age_bracket = case_when(age_bracket == "<1" ~ 0,
                                  T ~ as.numeric(age_bracket))) %>%
  group_by(iter, intervention, season_birth) %>%
  summarise(av_age = sum(age_bracket * n_hospitalisations)/sum(n_hospitalisations)) %>%
  ungroup() %>%
  group_by(intervention, season_birth) %>%
  summarise(av_age_low_95 = quantile(av_age, 0.025),
            av_age_median = quantile(av_age, 0.5),
            av_age_up_95 = quantile(av_age, 0.975)) %>%
  ungroup() %>%
  rbind(
    total_hosp_intervention_df %>%
      ungroup() %>%
      filter(age_bracket != "all" & age_bracket != "12-14") %>%
      mutate (age_bracket = case_when(age_bracket == "<1" ~ 0,
                                      T ~ as.numeric(age_bracket))) %>%
      group_by(iter, intervention) %>%
      summarise(season_birth = "all",
                av_age = sum(age_bracket * n_hospitalisations)/sum(n_hospitalisations)) %>%
      ungroup() %>%
      group_by(intervention, season_birth) %>%
      summarise(av_age_low_95 = quantile(av_age, 0.025),
                av_age_median = quantile(av_age, 0.5),
                av_age_up_95 = quantile(av_age, 0.975)) %>%
      ungroup()
  ) %>% arrange(season_birth)

# Get median age at hospitalisation
median_age_hosp <- total_hosp_intervention_df %>%
  ungroup() %>%
  filter(age_bracket != "all" & age_bracket != "12-14") %>%
  mutate (age_bracket = case_when(age_bracket == "<1" ~ 0,
                                  T ~ as.numeric(age_bracket))) %>%
  group_by(iter, intervention, season_birth) %>%
  summarise(median_age = spatstat.univar::weighted.median(age_bracket, n_hospitalisations)) %>%
  ungroup() %>%
  group_by(intervention, season_birth) %>%
  summarise(median_age_low_95 = quantile(median_age, 0.025),
            median_age_median = quantile(median_age, 0.5),
            median_age_up_95 = quantile(median_age, 0.975)) %>%
  ungroup() %>%
  rbind(
    total_hosp_intervention_df %>%
      ungroup() %>%
      filter(age_bracket != "all" & age_bracket != "12-14") %>%
      mutate (age_bracket = case_when(age_bracket == "<1" ~ 0,
                                      T ~ as.numeric(age_bracket))) %>%
      group_by(iter, intervention) %>%
      summarise(season_birth = "all",
                median_age = spatstat.univar::weighted.median(age_bracket, n_hospitalisations)) %>%
      ungroup() %>%
      group_by(intervention, season_birth) %>%
      summarise(median_age_low_95 = quantile(median_age, 0.025),
                median_age_median = quantile(median_age, 0.5),
                median_age_up_95 = quantile(median_age, 0.975)) %>%
      ungroup()
  ) %>% arrange(season_birth)


severe_cases_u1_de %>%
      group_by(age_months) %>%
      summarise(n_infections = sum(n_total_scaled),
                n_severe = sum(n_severe_scaled),
                prop_hosp_age = n_severe/n_infections) %>%
      mutate(season_birth = "all") %>%
  ggplot(aes(x = age_months, y = prop_hosp_age)) +
  geom_col(fill = "darkblue") +
  labs(x = "\nAge (months)",
       y = "Proportion of infections leading to hospitalisation\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
         axis.title.x = element_text(size = 30),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 30)) +
  scale_y_continuous(labels = scales::percent_format())

severe_cases_u1_de %>%
  group_by(season_birth, age_months) %>%
  summarise(n_infections = sum(n_total_scaled),
            n_severe = sum(n_severe_scaled)) %>%
  rbind(
    severe_cases_u1_de %>%
      group_by(age_months) %>%
      summarise(n_infections = sum(n_total_scaled),
            n_severe = sum(n_severe_scaled)) %>%
      mutate(season_birth = "all") 
    ) %>%
  ggplot(aes(x = age_months, y = n_severe, fill = season_birth)) +
  geom_col() +
  facet_wrap(~season_birth) +
  labs(x = "\nAge (months)",
       y = "Number of hospitalisations\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.title.x = element_text(size = 25),
         axis.text.y = element_text(size = 20),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 20, color = "black")) +
  scale_y_continuous(labels = scales::comma)

severe_cases_u1_de %>%
  group_by(season_birth, age_months) %>%
  summarise(n_infections = sum(n_total_scaled),
            n_severe = sum(n_severe_scaled)) %>%
  rbind(
    severe_cases_u1_de %>%
      group_by(age_months) %>%
      summarise(n_infections = sum(n_total_scaled),
                n_severe = sum(n_severe_scaled)) %>%
      mutate(season_birth = "all") 
  ) %>%
  ggplot(aes(x = age_months, y = n_infections, fill = season_birth)) +
  geom_col() +
  facet_wrap(~season_birth) +
  labs(x = "\nAge (months)",
       y = "Number of infections\n") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.title.x = element_text(size = 25),
         axis.text.y = element_text(size = 20),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 20, color = "black")) +
  scale_y_continuous(labels = scales::comma)

prop_hosp %>%
  filter (age_months == "<1" |
            age_months == "1" |
            age_months == "2" |
            age_months == "3" |
            age_months == "4" |
            age_months == "5" |
            age_months == "6" |
            age_months == "7" |
            age_months == "8" |
            age_months == "9" |
            age_months == "10" |
            age_months == "11" ) %>%
  ggplot(aes(x = age_months, y = prop_age)) +
  geom_point()

# % case that get hospitalised by age bracket
prop_hosp_age_bracket <- total_hosp_intervention_df %>%
  group_by(season_birth, age_bracket, iter, intervention, severity) %>%
  summarise(n_cases = sum(n_cases)) %>%
  ungroup() %>%
  select(season_birth, age_bracket, iter, n_cases, intervention, severity) %>%
  pivot_wider(names_from = severity,
              values_from = n_cases) %>%
  group_by(season_birth, age_bracket, iter, intervention) %>%
  mutate(total = severe + LRTI,
         prop_hosp_age = severe/total) %>%
  ungroup() %>%
  group_by(season_birth, age_bracket, intervention) %>%
  summarise(prop_hosp_low_95 = quantile(prop_hosp_age, 0.025, na.rm = T),
            prop_hosp_median = quantile(prop_hosp_age, 0.5, na.rm = T),
            prop_hosp_up_95 = quantile(prop_hosp_age, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth),
         season_birth = factor(season_birth, 
                              levels = c("Winter", "Spring", "Summer", "Autumn",
                                         "All")),
         age_bracket = factor(age_bracket, 
                             levels = c("<1", "1", "2", "3", "4", "5", "6", "7", "8", 
                                        "9", "10", "11", "12-14")))

plt_prop_hosp <- prop_hosp_age_bracket %>%
  filter(age_bracket != "all" & age_bracket != "12-14" & 
           intervention == "No immunisation") %>%
  ggplot(aes(x = age_bracket, y = prop_hosp_median, fill = season_birth)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = prop_hosp_low_95, ymax = prop_hosp_up_95, group = season_birth), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "\nAge (months)",
       y = "Proportion of cases leading to hospitalisation\n",
       fill = "Season of birth") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         legend.text = element_text(size = 25)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_season) +
  facet_wrap(~season_birth)

plt_prop_hosp 

plt_prop_hosp %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Prop hosp cases.png",
                    width = 14, height = 16, units = "in", 
                    device='png')

# Get % reduction in hosp rate between 0 and 6 months with CI
hosp_rate_diff <- severe_illness_new %>%
  filter(age_months == "<1" | age_months == "6") %>%
  select(age_months, total_severe_illness_rate) %>%
  pivot_wider(names_from = age_months,
              values_from = total_severe_illness_rate) %>%
  mutate(reduction = (`<1` - `6`)/`<1`) %>%
  rename(hosp_rate_u1 = `<1`,
         hosp_rate_6m = `6`) 
  # get 95% CI
  
  popEpi::rate_ratio(x = c(hosp_rate_diff$hosp_rate_6m, 100000), 
                     y = c(hosp_rate_diff$hosp_rate_u1, 100000), SE.method = FALSE)
  