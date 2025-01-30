# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)

# Path to files
path_output <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/mcstate/"  #Where model outputs are stored
path_pop <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/German population estimates/"

# ----------- Step 1: get number of cases by age and season in Germany ----------
# We use the number of births by month in Germany for 2023 and combine them with 
# our previously obtained estimates for the proportion of children who will
# develop RSV illness at a given age. 
# ------------------------------------------------------------------------------

# Read in model outputs
severe_illness <- read.csv(paste0(path_output, "Proportion severe illness by age.csv"))
mild_illness <- read.csv(paste0(path_output, "Proportion mild illness by age.csv"))

# Restrict to < 1 year
severe_illness_u1 <- severe_illness %>% filter (age_months <= 12)
mild_illness_u1 <- mild_illness %>% filter (age_months <= 12)

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

# Plot the numbers
severe_cases_u1_de %>% ggplot() +
  geom_point(aes(x = age_months, y = n_severe, colour = season_birth)) +
  labs(title = "Number of severe RSV cases", x = "Age (months)",
       y = "Cases", colour = "Season of birth") +
  theme_light()

mild_cases_u1_de %>% ggplot() +
  geom_point(aes(x = age_months, y = n_mild, colour = season_birth)) +
  labs(title = "Number of severe RSV cases", x = "Age (months)",
       y = "Cases", colour = "Season of birth") +
  theme_light()

severe_cases_u1_de %>% ggplot(aes(x = age_midpoint, y = n_severe)) +
  geom_point(aes(colour = current_season)) + 
  geom_smooth() +
  labs(title = "Number of severe RSV cases", x = "Age (days)",
       y = "Cases", colour = "Current season") +
  theme_light() +
  facet_wrap(~season_birth)

# ------------- Nirvesimab ----------------------------------------------------
# From Hammit et al, nirvesimab prevents 76.8% of RSV-associated hospitalisations.
# From the SA data, MA severe illness = hospitalisation
# Combine the two estimates
# ------------------------------------------------------------------------------

hosp_prevented <- severe_cases_u1_de %>% 
  select(season_birth, age_bracket, age_months, age_midpoint, current_season,
         n_MA_severe) %>%
  mutate(n_MA_severe_nirv = n_MA_severe * (1-0.768),
         n_MA_severe_averted = n_MA_severe * 0.768)

# Get number of cases averted by given nirvesimab at a given age
hosp_prevented_age <- hosp_prevented %>% group_by(season_birth) %>%
  summarise (prev_u1 = sum(n_MA_severe_averted),
             prev_1 = sum(n_MA_severe_averted[age_months >= 1]),
             prev_3 = sum(n_MA_severe_averted[age_months >= 3]),
             prev_4 = sum(n_MA_severe_averted[age_months >= 4]),
             prev_7 = sum(n_MA_severe_averted[age_months >= 7]),
             prev_8 = sum(n_MA_severe_averted[age_months >= 8]),
             prev_10 = sum(n_MA_severe_averted[age_months >= 10]),
             prev_12 = sum(n_MA_severe_averted[age_months >= 12])) %>%
  ungroup() %>%
  rbind(
    data.frame(season_birth = "all",
               prev_u1 = sum(hosp_prevented$n_MA_severe_averted),
               prev_1 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 1]),
               prev_3 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 3]),
               prev_4 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 4]),
               prev_7 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 7]),
               prev_8 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 8]),
               prev_10 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 10]),
               prev_12 = sum(hosp_prevented$n_MA_severe_averted[hosp_prevented$age_months >= 12]))
  )

# Pivot longer for easier plotting
hosp_prevented_age_long <- hosp_prevented_age %>% 
  rename(`0` = "prev_u1",
         `1` = "prev_1",
         `3` = "prev_3",
         `4` = "prev_4",
         `7` = "prev_7",
         `8` = "prev_8",
         `10` = "prev_10",
         `12` = "prev_12") %>%
  pivot_longer(!season_birth,
               names_to = "immun_age", values_to = "hosp_prev") %>%
  mutate(immun_age = factor(immun_age, levels = c("0", "1", "3", "4", "7", "8", "10", "12")))

hosp_prevented_age_long %>% ggplot() +
  geom_col(aes(x = immun_age, y = hosp_prev, fill = season_birth), position="dodge") +
  labs(title = "Number of averted RSV-associated hospitalisation by age at immunisation",
       x = "\nAge at immunisation (months)",
       y = "Averted hospitalisations\n",
       fill = "Season of birth") +
  theme_light()
