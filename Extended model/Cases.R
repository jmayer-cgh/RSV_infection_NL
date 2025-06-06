# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)

# Path to files
path_output <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"  #Where model outputs are stored
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
                                    current_season = "all",
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
                                    current_season = "all",
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

# Plot the numbers
severe_cases_u1_de %>% filter(!is.na(age_months)) %>% ggplot() +
  geom_point(aes(x = age_months, y = n_severe, colour = season_birth)) +
  labs(title = "Number of severe RSV cases", x = "Age (months)",
       y = "Cases", colour = "Season of birth") +
  theme_light()

mild_cases_u1_de %>% filter(!is.na(age_months)) %>% ggplot() +
  geom_point(aes(x = age_months, y = n_mild, colour = season_birth)) +
  labs(title = "Number of mild RSV cases", x = "Age (months)",
       y = "Cases", colour = "Season of birth") +
  theme_light()

pain_au_lait <- severe_cases_u1_de %>% filter(!is.na(age_months)) %>% 
  ggplot(aes(x = age_midpoint, y = n_severe)) +
  geom_point(aes(colour = current_season), size = 3.5) + 
  geom_smooth() +
  labs(title = "Number of severe RSV cases", x = "Age (days)",
       y = "Cases", colour = "Current season") +
  theme_light() +
  theme (axis.text.y=element_blank(), 
         axis.ticks.y=element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         title = element_text(size = 20),
         legend.text = element_text (size = 25),
         legend.title = element_text (size = 25),
         strip.text.x = element_text(size = 20, color = "black")) +
  facet_wrap(~season_birth) 

pain_au_lait %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Severe cases DE by season of birth and age.png",
         width = 14, height = 16, units = "in", 
         device='png')

plt <- mild_cases_u1_de %>% filter(!is.na(age_months)) %>%
  ggplot(aes(x = age_midpoint, y = n_mild)) +
  geom_point(aes(colour = current_season), size = 3.5) + 
  geom_smooth() +
  labs(title = "Number of mild RSV cases", x = "Age (days)",
       y = "Cases", colour = "Current season") +
  theme_light() +
  theme (axis.text.y=element_blank(), 
         axis.ticks.y=element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         title = element_text(size = 20),
         legend.text = element_text (size = 18),
         strip.text.x = element_text(size = 20, color = "black")) +
  facet_wrap(~season_birth) 

plt %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Mild cases DE by season of birth and age.png",
         width = 14, height = 16, units = "in", 
         device='png')

severe_cases_u1_de %>% filter(!is.na(age_months)) %>% 
  group_by(age_midpoint) %>%
  mutate(total = sum(n_severe_scaled)) %>%
  ggplot(aes(x = age_months, y = total)) +
  geom_point(size = 3.5) + 
  labs(title = "Number of severe RSV cases", x = "Age (months)",
       y = "Cases") +
  theme_light() +
  theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         title = element_text(size = 20),
         legend.text = element_text (size = 25),
         legend.title = element_text (size = 25),
         strip.text.x = element_text(size = 20, color = "black"))

severe_illness_u1 %>% filter(!is.na(age_months)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("spring", "summer", "autumn", "winter"))) %>%
  ggplot(aes(x = age_months, y = total_ma_severe_cases_prop)) +
  geom_point(size = 3.5) + 
  labs(title = "Proportion of severe disease among RSV cases", x = "Age (months)",
       y = "Proportion of cases") +
  theme_light() +
  theme (axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         title = element_text(size = 20),
         legend.text = element_text (size = 25),
         legend.title = element_text (size = 25),
         strip.text.x = element_text(size = 20, color = "black")) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~season_birth, scales = "free")

# ------------- Abrysvo ----------------------------------------------------
# We have estimates of VE by age so we will use those
# ------------------------------------------------------------------------------
# Estimate VE at different ages from Fabienne's S6 figure
# Hospitalisations
# hosp_prevented_vacc <- severe_cases_u1_de %>% filter(!is.na(age_months)) %>%
#   select(season_birth, age_bracket, age_months, age_midpoint, current_season,
#          n_ma_severe_scaled) %>%
#   mutate(ve = case_when(age_months == 0 ~ 0.75, # 0.87 from Ayaka // 0.62 from Hodgson
#                         age_months == 1 ~ 0.70, # 0.83 // 0.62
#                         age_months == 3 ~ 0.50, # 0.66 // 0.5
#                         age_months == 5 ~ 0.25, # 0.45 // 0.4
#                         age_months == 7 ~ 0.125, # 0.32 // 0.3
#                         age_months == 8 ~ 0.10, # 0.26 // 0.25
#                         age_months == 9 ~ 0.06, # 0.17 // 0.14
#                         age_months == 11 ~ 0.04)) # 0.10 // 0.10

# Read in numbers from Akaya's code
VE_distribution <- read.csv("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Erlang waning/Model outputs.csv")

# Draw a random number from the distribution
VE_distribution <- VE_distribution %>% group_by (group, t) %>%
  mutate(VE_value = runif(1, min = VE_t_lo, max = VE_t_hi)) %>%
  ungroup()

# Add the estimated VE to the number of cases
hosp_prevented_vacc <- severe_cases_u1_de %>% filter(!is.na(age_months)) %>%
  select(season_birth, age_bracket, age_months, age_midpoint, current_season,
         n_ma_severe_scaled) %>%
  merge(
    VE_distribution %>% filter(group == "severe") %>%
      select(t, VE_value),
    by.x = "age_midpoint", by.y = "t"
  )

# Get number of prevented cases
hosp_prevented_vacc <- hosp_prevented_vacc %>% 
  group_by(season_birth, age_bracket, age_months, age_midpoint, current_season, VE_value) %>%
  summarise(n_ma_severe_vacc = n_ma_severe_scaled * (1-VE_value), # number of cases despite vaccination
            n_ma_severe_averted = n_ma_severe_scaled * VE_value)

hosp_prevented_vacc <- hosp_prevented_vacc %>% 
  rbind(
    hosp_prevented_vacc %>% 
      group_by(season_birth) %>%
      summarise(age_bracket = 'all', 
                age_months = NA,
                age_midpoint = NA, 
                current_season = NA,
                n_ma_severe_vacc = sum(n_ma_severe_vacc),
                n_ma_severe_averted = sum(n_ma_severe_averted))
  )

# Medically assisted cases -- not sure we can actually get these
# First, we need to combine mild and severe MA cases
ma_cases_u1_de <- mild_cases_u1_de %>% filter(!is.na(age_months)) %>%
  select (season_birth, current_season, age_months, age_midpoint, n_ma_mild_scaled) %>%
  mutate(severity = "mild") %>%
  rename(n_ma_cases = "n_ma_mild_scaled") %>%
  rbind(
    severe_cases_u1_de %>% filter(!is.na(age_months)) %>%
      select (season_birth, current_season, age_months, age_midpoint, n_ma_severe_scaled) %>%
      mutate(severity = "severe") %>%
      rename(n_ma_cases = "n_ma_severe_scaled")
  )

ma_cases_u1_de <- ma_cases_u1_de %>%
  rbind(
    ma_cases_u1_de %>% group_by(season_birth, current_season, age_months, age_midpoint) %>%
      summarise (n_ma_cases = sum(n_ma_cases),
                 severity = "all")
  )

# Now we can get the # of prevented cases
# ma_prevented_vacc <- ma_cases_u1_de %>% 
#   select(season_birth, age_months, n_ma_cases, severity) %>%
#   mutate(ve = case_when(age_months == 0 ~ 0.75, # 0.65 from Ayaka // 0.9 from Hodgson
#                         age_months == 1 ~ 0.70, # 0.62 // 0.87
#                         age_months == 3 ~ 0.50, # 0.50 // 0.8
#                         age_months == 5 ~ 0.25, # 0.34 // 0.65
#                         age_months == 7 ~ 0.125, # 0.24 // 0.45
#                         age_months == 8 ~ 0.10, # 0.20 // 0.43
#                         age_months == 9 ~ 0.06, # 0.13 // 0.42
#                         age_months == 11 ~ 0.04)) # 0.08 // 0.29

ma_prevented_vacc <- ma_cases_u1_de %>% filter(!is.na(age_months)) %>%
  select(season_birth, age_months, age_midpoint, current_season,
         n_ma_cases, severity) %>%
  merge(
    VE_distribution %>%
      select(t, VE_value, group) %>%
      mutate(group = case_when (group == "LRTI" ~ "mild",
                                TRUE ~ group)),
    by.x = c("age_midpoint","severity"), by.y = c("t", "group")
  )

ma_prevented_vacc <- ma_prevented_vacc %>% 
  group_by(season_birth, age_months, severity, VE_value) %>%
  summarise(n_ma_cases_vacc = n_ma_cases * (1-VE_value), # number of cases despite vaccination
            n_ma_cases_averted = n_ma_cases * VE_value)

ma_prevented_vacc <- ma_prevented_vacc %>% 
  rbind(
    ma_prevented_vacc %>% 
      group_by(season_birth) %>%
      summarise(age_months = NA,
                severity = "all",
                n_ma_cases_vacc = sum(n_ma_cases_vacc, na.rm = T),
                n_ma_cases_averted = sum(n_ma_cases_averted, na.rm = T))
  )

# ------------- Nirvesimab ----------------------------------------------------
# From Muller et al, nirvesimab prevents 76.8% (62.3–85.2) of RSV-associated hospitalisations.
# From the SA data, MA severe illness = hospitalisation
# We want to know the added effect of nirsevimab with maternal vaccination
# Combine the two estimates
# ------------------------------------------------------------------------------
nirs_eff_hosp <- runif(1, 0.623, 0.852)

hosp_prevented <- hosp_prevented_vacc %>% filter(!is.na(age_months)) %>%
  select(season_birth, age_bracket, age_months, age_midpoint, current_season,
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

# Pivot longer for easier plotting
hosp_prevented_age_nirs_long <- hosp_prevented_age %>% 
  select(season_birth, prev_u1, prev_1, prev_3, prev_4, prev_7, prev_8, prev_10, 
         prev_11, prev_spring_10, prev_summer_4) %>%
  rename(`0` = "prev_u1",
         `1` = "prev_1",
         `3` = "prev_3",
         `4` = "prev_4",
         `7` = "prev_7",
         `8` = "prev_8",
         `10` = "prev_10",
         `11` = "prev_11",
         `10_spring` = "prev_spring_10",
         `4_summer` = "prev_summer_4") %>%
  pivot_longer(!season_birth,
               names_to = "immun_age", values_to = "hosp_prev") %>%
  mutate(immun_age = factor(immun_age, levels = c("0", "1", "3", "4", "7", "8", 
                                                  "10", "11", "10_spring", "4_summer")))

hosp_prevented_age_nirs_long  %>% ggplot() +
  geom_col(aes(x = immun_age, y = hosp_prev, fill = season_birth), position="dodge") +
  labs(title = "Number of averted RSV-associated hospitalisation by age at immunisation",
       x = "\nAge at immunisation (months)",
       y = "Averted hospitalisations\n",
       fill = "Season of birth") +
  theme_light()

# From Hammit et al, nirvesimab prevents 74.5% (49.6 to 87.1) of medically attended RSV-associated lower respiratory tract infections
# Calculate this
nirs_eff_LRTI <- runif(1, 0.496, 0.871)

# Get number of prevented cases
ma_prevented_nirs <- ma_prevented_vacc %>% filter(!is.na(age_months)) %>%
  select(season_birth, age_months, n_ma_cases_vacc, severity) %>%
  mutate(n_ma_cases_nirs = n_ma_cases_vacc * (1-nirs_eff_LRTI),
         n_ma_cases_averted_nirs = n_ma_cases_vacc * nirs_eff_LRTI)

ma_prevented_nirs <- ma_prevented_nirs %>%
  rbind(
    ma_prevented_nirs %>%
      group_by(season_birth, age_months) %>%
      mutate(n_ma_cases_vacc = sum(n_ma_cases_vacc),
             severity = "all",
             n_ma_cases_nirs = sum(n_ma_cases_nirs),
             n_ma_cases_averted_nirs = sum(n_ma_cases_averted_nirs)) %>%
      unique()
  ) %>%
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
             prev_11 = sum(n_ma_cases_averted_nirs[age_months >= 11], na.rm = T)) %>%
  ungroup() %>%
  rbind(
    data.frame(season_birth = "all",
               severity = "all",
               prev_u1 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$severity == "all"]),
               prev_1 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 1 & ma_prevented_nirs$severity == "all"]),
               prev_3 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 3 & ma_prevented_nirs$severity == "all"]),
               prev_4 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 4 & ma_prevented_nirs$severity == "all"]),
               prev_7 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 7 & ma_prevented_nirs$severity == "all"]),
               prev_8 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 8 & ma_prevented_nirs$severity == "all"]),
               prev_10 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 10 & ma_prevented_nirs$severity == "all"]),
               prev_11 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 11 & ma_prevented_nirs$severity == "all"]))
  ) %>%
  rbind(
    data.frame(season_birth = "all",
               severity = "mild",
               prev_u1 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$severity == "mild"]),
               prev_1 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 1 & ma_prevented_nirs$severity == "mild"]),
               prev_3 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 3 & ma_prevented_nirs$severity == "mild"]),
               prev_4 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 4 & ma_prevented_nirs$severity == "mild"]),
               prev_7 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 7 & ma_prevented_nirs$severity == "mild"]),
               prev_8 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 8 & ma_prevented_nirs$severity == "mild"]),
               prev_10 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 10 & ma_prevented_nirs$severity == "mild"]),
               prev_11 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 11 & ma_prevented_nirs$severity == "mild"]))
  )%>%
  rbind(
    data.frame(season_birth = "all",
               severity = "severe",
               prev_u1 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$severity == "severe"]),
               prev_1 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 1 & ma_prevented_nirs$severity == "severe"]),
               prev_3 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 3 & ma_prevented_nirs$severity == "severe"]),
               prev_4 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 4 & ma_prevented_nirs$severity == "severe"]),
               prev_7 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 7 & ma_prevented_nirs$severity == "severe"]),
               prev_8 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 8 & ma_prevented_nirs$severity == "severe"]),
               prev_10 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 10 & ma_prevented_nirs$severity == "severe"]),
               prev_11 = sum(ma_prevented_nirs$n_ma_cases_averted_nirs[ma_prevented_nirs$age_months >= 11 & ma_prevented_nirs$severity == "severe"]))
  ) %>%
  ungroup()

# Pivot longer for easier plotting
ma_prevented_age_long_nirs <- ma_prevented_age_nirs %>% 
  rename(`0` = "prev_u1",
         `1` = "prev_1",
         `3` = "prev_3",
         `4` = "prev_4",
         `7` = "prev_7",
         `8` = "prev_8",
         `10` = "prev_10",
         `11` = "prev_11") %>%
  pivot_longer(!c(season_birth, severity),
               names_to = "immun_age", values_to = "ma_prev") %>%
  mutate(immun_age = factor(immun_age, levels = c("0", "1", "3", "4", "7", "8", "10", "11")))

ma_prevented_age_long_nirs %>% ggplot() +
  geom_col(aes(x = immun_age, y = ma_prev, fill = season_birth), position="dodge") +
  labs(title = "Number of averted RSV-associated MA-cases by age at immunisation",
       x = "\nAge at immunisation (months)",
       y = "Averted MA cases\n",
       fill = "Season of birth") +
  theme_light() +
  facet_wrap(~severity)

# Plot total hospitalisation number by intervention
# Get the numbers
total_hosp_intervention <- severe_cases_u1_de %>%
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

# Check calculation of remaining hospitalisations with nirsevimab
remaining_hosp_nirs <- total_hosp_intervention$n_hospitalisations[total_hosp_intervention$intervention == "Maternal vaccination"] - total_hosp_intervention$prevented_hospitalisations[total_hosp_intervention$intervention == "Nirsevimab for spring births"]
round(remaining_hosp_nirs) == round(total_hosp_intervention %>% filter(intervention == "Nirsevimab for spring births") %>% select (n_hospitalisations))

remaining_hosp_nirs <- total_hosp_intervention$n_hospitalisations[total_hosp_intervention$intervention == "Maternal vaccination"] - total_hosp_intervention$prevented_hospitalisations[total_hosp_intervention$intervention == "Nirsevimab for summer births"]
round(remaining_hosp_nirs) == round(total_hosp_intervention %>% filter(intervention == "Nirsevimab for summer births") %>% select (n_hospitalisations))

# Plot
palette <- c("No immunisation" = "#9C964A", "Maternal vaccination" = "#6A9D96",
             "Nirsevimab for spring births" = "#88BBA0",
             "Nirsevimab for summer births" = "#85D4E3",
             "Nirsevimab for spring and summer births" = "#B39BC8")

total_hosp_intervention_plt <- total_hosp_intervention %>% ggplot() +
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

total_hosp_intervention_plt %>% ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Outputs/Hosp by intervention.png",
                                       width = 14, height = 16, units = "in", 
                                       device='png')

# Save outputs
# Estimated cases
severe_cases_u1_de %>% write.csv(paste0(path_output, "Severe cases DE.csv"), row.names = F)
mild_cases_u1_de %>% write.csv(paste0(path_output, "Mild cases DE.csv"), row.names = F)

# Vaccination
hosp_prevented_vacc %>% write.csv(paste0(path_output, "Prevented hospitalisations vaccine.csv"), row.names = F)
ma_prevented_vacc %>% write.csv(paste0(path_output, "Prevented MA cases vaccine.csv"), row.names = F)

# mAB
hosp_prevented_age %>% write.csv(paste0(path_output, "Prevented hospitalisations nirsevimab.csv"), row.names = F)
ma_prevented_age_nirs %>% write.csv(paste0(path_output, "Prevented MA cases nirsevimab.csv"), row.names = F)
s