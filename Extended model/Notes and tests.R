# Alternative LL
N_spring <- data()
n_infection_spring <- data()
N_summer <- data()
n_infection_summer <- data()
N_autumn <- data()
n_infection_autumn <- data()
N_winter <- data()
n_infection_winter <- data()

N_tot <- N_spring + N_summer + N_autumn + N_winter
n_infection <- n_infection_spring + n_infection_summer + n_infection_autumn + n_infection_winter
n_infection ~ Binomial(N_tot, R_all)

# ------------------------------------------------------------------------------
# Plotting CIs
# ------------------------------------------------------------------------------

data <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/master/data/infection_status_csv.txt",
                 sep=",")

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
    Birth_mo = birthday %>% lubridate::month(),
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
          prop_seroconv = n_infection/N,
          cum_infection = prop_seroconv * N_tot) %>%
  ungroup() %>% 
  distinct()

incidence_data_season <- data %>% select (age_grp, age_days, infection, season_birth, xMidpoint) %>%
  mutate(N_tot = n()) %>%
  group_by(age_grp, xMidpoint, season_birth) %>%
  summarise(age_mid = round(median(age_days)), 
            N = n(),
            n_infection = sum(infection),
            prop_seroconv = n_infection/N,
            cum_infection = prop_seroconv * N_tot) %>%
  ungroup() %>% 
  distinct()

# Calculate seroprevalence and binomial confidence intervals
incidence_data[,c("seroprev_mean","lower_CI","upper_CI")] <- binom::binom.confint(incidence_data$n_infection, incidence_data$N, method="exact")[,c("mean","lower","upper")]


incidence_data_season[,c("seroprev_mean","lower_CI","upper_CI")] <- binom::binom.confint(incidence_data_season$n_infection, 
                                                                                         incidence_data_season$N, 
                                                                                         method="exact")[,c("mean","lower","upper")]


incidence_data_season_wide <- incidence_data_season %>% select (!c(age_mid, age_grp, seroprev_mean)) %>%
  pivot_wider(
    names_from = season_birth,
    values_from = c(N , n_infection, prop_seroconv, lower_CI, upper_CI, cum_infection),
    values_fill = 0
  ) %>%
  rename(time = "xMidpoint")


# Build a filter and test it
filter <- dust_filter_create(msr, time_start = 0, data = incidence_data_season_wide, n_particles = 1000)
dust_likelihood_run(filter, list(spring_comp = 0.003478,
                                 summer_comp = 0.00060,
                                 autumn_comp = 0.00264,
                                 winter_comp = 0.00692, 
                                 mu = 0.09,
                                 prop = 0.3),
                    save_trajectories = T)

# Plot the fit
h <- dust_likelihood_last_trajectories(filter)
matplot(incidence_data$time, t(dust_unpack_state(filter, h)$R_all), type = "l",
        lty = 1, col = "#00000044",
        xlab = "Age", ylab = "Proportion seroconverted")
points(prop_seroconv ~ time, incidence_data, pch = 19, col = "red")
arrows(x0 = incidence_data$time, y0 = incidence_data$lower_CI, 
       x1 = incidence_data$time, y1 = incidence_data$upper_CI, 
       angle = 90, code = 3, length = 0.1, col = "blue")

# --------- Checking seroconversion predictions -------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"

seroconversion_df <- read.csv(paste0(path, "seroconversion by age.csv"))
seroconversion_df <- seroconversion_df %>%
  mutate(median_all = 0.26*median_sp + 
           0.29*median_sm + 
           0.24*median_au + 
           0.20*median_wt)

# Get incident cases
incidence_df <- data.frame(age_midpoint = seroconversion_df$age_midpoint,
                           incidence_median = c(0, diff(seroconversion_df$median_sp)),
                           incidence_low95 = c(0, diff(seroconversion_df$low95_sp)), # is this right?
                           incidence_up95 = c(0, diff(seroconversion_df$up95_sp)),
                           season_birth = "spring") %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence_median = c(0, diff(seroconversion_df$median_sm)),
               incidence_low95 = c(0, diff(seroconversion_df$low95_sm)), # is this right?
               incidence_up95 = c(0, diff(seroconversion_df$up95_sm)),
               season_birth = "summer")
  ) %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence_median = c(0, diff(seroconversion_df$median_au)),
               incidence_low95 = c(0, diff(seroconversion_df$low95_au)), # is this right?
               incidence_up95 = c(0, diff(seroconversion_df$up95_au)),
               season_birth = "autumn")
  ) %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence_median = c(0, diff(seroconversion_df$median_wt)),
               incidence_low95 = c(0, diff(seroconversion_df$low95_wt)), # is this right?
               incidence_up95 = c(0, diff(seroconversion_df$up95_wt)),
               season_birth = "winter")
  ) %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence_median = c(0, diff(seroconversion_df$median_all)),
               incidence_low95 = NA,
               incidence_up95 = NA,
               season_birth = "All")
  )

seroconversion_df %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_sp)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_sp, ymax = up95_sp)) +
  labs(title = "Proportion of seroconverted children (spring)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

seroconversion_df %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_sm)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_sm, ymax = up95_sm)) +
  labs(title = "Proportion of seroconverted children (summer)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

seroconversion_df %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_au)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_au, ymax = up95_au)) +
  labs(title = "Proportion of seroconverted children (autumn)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

seroconversion_df %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_wt)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_wt, ymax = up95_wt)) +
  labs(title = "Proportion of seroconverted children (winter)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()


seroconversion_df %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_all)) +
  labs(title = "Proportion of seroconverted children (all)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

incidence_df %>% filter(age_midpoint <= 365 & season_birth != "All") %>%
  ggplot()+
  geom_point(aes(x = age_midpoint, y = incidence_median, col = season_birth)) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = incidence_low95, ymax = incidence_up95,
                    col = season_birth)) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversion", col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth)

incidence_df %>% filter(age_midpoint <= 500 & season_birth == "All") %>%
  ggplot()+
  geom_point(aes(x = age_midpoint, y = incidence_median, col = season_birth)) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth)