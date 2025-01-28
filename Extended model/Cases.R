# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)

# Path to files
path_output <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/mcstate/"  #Where model outputs are stored
path_pop <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/"

# Read in model outputs
severe_illness <- read.csv(paste0(path_output, "Proportion severe illness by age.csv"))
mild_illness <- read.csv(paste0(path_output, "Proportion mild illness by age.csv"))

# Restrict to < 1 year
severe_illness_u1 <- severe_illness %>% filter (age_months <= 12)
mild_illness_u1 <- mild_illness %>% filter (age_months <= 12)

# Read in German population estimates and format it
# pop_de <- read_excel(paste0(path_pop, "Population estimates.xlsx")) 
# colnames(pop_de)[1] <- "Age"
# colnames(pop_de)[2] <- "Population size"
# 
# pop_de <- pop_de %>% filter(grepl("year", Age)) %>%
#   mutate(Age = case_when(Age == "under 1 year" ~ 0,
#                          Age == "1 year" ~ 1,
#                          TRUE ~ as.numeric(gsub("years", "", Age))),
#          Age = as.numeric(Age),
#          `Population size` = as.numeric(`Population size`))

# Read in birth numbers for 2023
births_de <- read_excel(paste0(path_pop, "Births.xlsx")) 
colnames(births_de)[2] <- "Month"
colnames(births_de)[5] <- "Value"
births_de <- births_de %>% select(Month, Value) %>% 
  filter(grepl(".0", Value)) %>%
  mutate(Value = as.numeric(Value))

# Define birth cohorts
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)
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
         n_MA_severe = total_MA_severe_cases_prop * total,
         n_nonMA_severe = total_nonMA_severe_cases_prop * total)

mild_cases_u1_de <- mild_cases_u1_de %>%
  mutate(n_severe = total_mild_cases_prop * total,
         n_MA_severe = total_MA_mild_cases_prop * total,
         n_nonMA_severe = total_nonMA_mild_cases_prop * total)
