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
# Read in and the data and put it in the right format
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

# Plot number of children by age and season
incidence_data_season %>% filter (xMidpoint <= 365) %>%
  mutate(season_birth = factor(season_birth, levels = c("autumn", "winter", "spring", "summer"))) %>%
  ggplot()+
  geom_bar(aes(x = xMidpoint/30, y = N, fill = season_birth), 
           stat = "identity", position = "dodge") +
  scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11, 13),
                     labels=c("1","3","5", "7", "9", "11", "13")) +
  labs(title = "Number of children by age and season of birth", x = "Age (months)",
       y = "Number of children", fill = "Season of birth") +
  theme_light()

# Plot % of children by age and season
plt_child_prop <- incidence_data_season %>% filter (xMidpoint <= 365) %>%
  mutate(season_birth = factor(season_birth, levels = c("autumn", "winter", "spring", "summer"))) %>%
  group_by(xMidpoint) %>%
  mutate(N_age = sum(N),
         prop = N/N_age) %>%
  ungroup() %>%
  ggplot()+
  geom_bar(aes(x = xMidpoint/30, y = prop, fill = season_birth), 
           stat = "identity", position = "dodge") +
  scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11, 13),
                     labels=c("1","3","5", "7", "9", "11", "13")) +
  labs(title = "Proportion of children by age and season of birth", x = "Age midpoint (months)",
       y = "Proprtion of children", fill = "Season of birth") +
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         legend.text = element_text(size = 20))

plt_child_prop

plt_child_prop  %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Prop birth cohort by age.png",
         width = 20, height = 14, units = "in", 
         device='png')

# Plot seroprevalence by age and season
plt_season <- incidence_data_season %>% 
  filter (xMidpoint <= 365) %>%
  mutate(season_birth = factor(season_birth, levels = c("autumn", "winter", "spring", "summer"))) %>%
  ggplot()+
  geom_point(aes(x = xMidpoint, y = seroprev_mean, col = season_birth), size = 3) +
  geom_errorbar(aes(x = xMidpoint, ymin = seroprev_low95, ymax = seroprev_up95,
                    col = season_birth), width = 10, linewidth = 1) +
  labs(title = "Proportion of seroconverted children by age and season of birth", x = "Age (days)",
       y = "% seroconverted\n", col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth, scales = "free") +
  theme (axis.ticks.y = element_blank(),
       legend.position = "none",
       axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
       axis.text.y = element_text(size = 20),
       axis.title.x = element_text(size = 25),
       axis.title.y = element_text(size = 25),
       title = element_text(size = 25),
       strip.text.x = element_text(size = 25, color = "black"))

plt_season

plt_season  %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Seroconversion by season of birth and age data.png",
         width = 20, height = 14, units = "in", 
         device='png')

plt_season_all <- incidence_data_season %>% 
  mutate(season_birth = factor(season_birth, levels = c("autumn", "winter", "spring", "summer"))) %>%
  ggplot()+
  geom_point(aes(x = xMidpoint, y = seroprev_mean, col = season_birth), size = 3) +
  geom_errorbar(aes(x = xMidpoint, ymin = seroprev_low95, ymax = seroprev_up95,
                    col = season_birth), width = 10, linewidth = 1) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "orange", size = 1.5) +
  annotate("text", x = 900, y = 0.01, label = "Switch from IgA to IgG", color = "orange", size = 8) +
  labs(title = "Proportion of seroconverted children by age and season of birth", x = "Age (days)",
       y = "% seroconverted\n", col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth, scales = "free") +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))

plt_season_all

plt_season_all  %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Seroconversion by season of birth and age data with Ig.png",
         width = 20, height = 14, units = "in", 
         device='png')

# Plot incidence by age and season
plt_season_inc <- incidence_data_season %>% 
  filter (xMidpoint <= 365) %>%
  mutate(season_birth = factor(season_birth, levels = c("autumn", "winter", "spring", "summer"))) %>%
  ggplot()+
  geom_point(aes(x = xMidpoint, y = incidence_mean, col = season_birth), size = 3) +
  geom_errorbar(aes(x = xMidpoint, ymin = incidence_low95, ymax = incidence_up95,
                    col = season_birth), width = 10, linewidth = 1) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversion\n", col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))

plt_season_inc

# ------------------------------------------------------------------------------
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

seroconversion_df %>% 
  ggplot()+
  geom_point(aes(x = age_midpoint, y = median_wt)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_wt, ymax = up95_wt)) +
  labs(title = "Proportion of seroconverted children (winter)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()


seroconversion_df %>% filter (age_midpoint <= 395) %>%
  ggplot()+
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
  facet_wrap(~season_birth, scales = "free")

incidence_df %>% filter(age_midpoint <= 500 & season_birth == "All") %>%
  ggplot()+
  geom_point(aes(x = age_midpoint, y = incidence_median, col = season_birth)) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth)

# -----------------------------------------------------------------------------
incidence_data %>% filter(time <= 365) %>%
  ggplot() + 
  geom_point(aes(x = time, y = incidence_mean)) +
  geom_errorbar(aes(x = time, 
                    ymin = incidence_low95, ymax = incidence_up95)) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

plt_incidence <- incidence_data_season %>% 
  filter(xMidpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = xMidpoint, y = incidence_mean, col = season_birth), size = 3) +
  geom_errorbar(aes(x = xMidpoint, 
                    ymin = incidence_low95, ymax = incidence_up95,
                    col = season_birth), width = 10, linewidth = 1) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversions\n",
       col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))

plt_incidence

plt_incidence %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Incidence by season of birth and age data.png",
         width = 22, height = 14, units = "in", 
         device='png')

# Seroconversion estimates
converted_summary <- readRDS("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/RDS files/seroconversion model outputs.rds")
converted_all <- do.call(rbind.data.frame, converted_summary[1])
converted_sp <- do.call(rbind.data.frame, converted_summary[2])
converted_sm <- do.call(rbind.data.frame, converted_summary[3])
converted_au <- do.call(rbind.data.frame, converted_summary[4])
converted_wt <- do.call(rbind.data.frame, converted_summary[5])

seroconversion <- list()

for (i in 1:length(incidence_data_season_wide$time)){
  age_midpoint <- incidence_data_season_wide$time[i]
  
  low95_sp <- quantile(t(converted_sp[i,]), 0.05, na.rm = T)
  median_sp <- quantile(t(converted_sp[i,]), 0.5, na.rm = T)
  up95_sp <- quantile(t(converted_sp[i,]), 0.95, na.rm = T)
  
  low95_sm <- quantile(t(converted_sm[i,]), 0.05, na.rm = T)
  median_sm <- quantile(t(converted_sm[i,]), 0.5, na.rm = T)
  up95_sm <- quantile(t(converted_sm[i,]), 0.95, na.rm = T)
  
  low95_au <- quantile(t(converted_au[i,]), 0.05, na.rm = T)
  median_au <- quantile(t(converted_au[i,]), 0.5, na.rm = T)
  up95_au <- quantile(t(converted_au[i,]), 0.95, na.rm = T)
  
  low95_wt <- quantile(t(converted_wt[i,]), 0.05, na.rm = T)
  median_wt <- quantile(t(converted_wt[i,]), 0.5, na.rm = T)
  up95_wt <- quantile(t(converted_wt[i,]), 0.95, na.rm = T)
  
  low95_all <- quantile(t(converted_all[i,]), 0.05, na.rm = T)
  median_all <- quantile(t(converted_all[i,]), 0.5, na.rm = T)
  up95_all <- quantile(t(converted_all[i,]), 0.95, na.rm = T)
  
  combined <- data.frame(age_midpoint = age_midpoint,
                         low95_sp, median_sp, up95_sp, 
                         low95_sm, median_sm, up95_sm, 
                         low95_au, median_au, up95_au, 
                         low95_wt, median_wt, up95_wt,
                         low95_all, median_all, up95_all)
  
  seroconversion[[i]] <- combined
}

seroconversion_df <- do.call("rbind", seroconversion)
rownames(seroconversion_df) <- NULL

seroconversion_long <- seroconversion_df %>%
  pivot_longer(
    cols = -age_midpoint,  # keep age as is
    names_to = c("measure", "season"),  # split column names into two parts
    names_sep = "_",  # separator is underscore (e.g. "low95_sp")
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = measure,  # now spread low95, median, up95 into columns
    values_from = value
  ) %>%
  mutate(season = case_when(season == "sp" ~ "Spring",
                            season == "sm" ~ "Summer",
                            season == "wt" ~ "Winter",
                            season == "au" ~ "Autumn",
                            season == "all" ~ "All"),
         season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter", "All")))

plt <- seroconversion_long %>% filter(age_midpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median, col = season), size = 3) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95, ymax = up95,
                    col = season),width = 10, linewidth = 1) +
  labs(title = "Proportion of seroconverted children", x = "Age (days)",
       y = "% seroconverted\n",
       col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))
plt

plt %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Proportion seroconversion by season of birth and age.png",
         width = 20, height = 14, units = "in", 
         device='png')

# Get incident cases
incidence <- list()

for (i in 1:length(incidence_data_season_wide$time)){
  age_midpoint <- incidence_data_season_wide$time[i]
  
  if(i>1){
    low95_spring <- quantile(
                        (t(converted_sp[i,]) - t(converted_sp[i-1,])), 0.05)
    median_spring <- quantile(
                          (t(converted_sp[i,]) - t(converted_sp[i-1,])), 0.5)
    up95_spring <- quantile(
      (t(converted_sp[i,]) - t(converted_sp[i-1,])), 0.95)
    
    low95_summer <- quantile((t(converted_sm[i,]) - t(converted_sm[i-1,])), 0.05)
    median_summer <- quantile((t(converted_sm[i,]) - t(converted_sm[i-1,])), 0.5)
    up95_summer <- quantile((t(converted_sm[i,]) - t(converted_sm[i-1,])), 0.95)
    
    low95_autumn <- quantile((t(converted_au[i,]) - t(converted_au[i-1,])), 0.05)
    median_autumn <- quantile((t(converted_au[i,]) - t(converted_au[i-1,])), 0.5)
    up95_autumn <- quantile((t(converted_au[i,]) - t(converted_au[i-1,])), 0.95)
    
    low95_winter <- quantile((t(converted_wt[i,]) - t(converted_wt[i-1,])), 0.05)
    median_winter <- quantile((t(converted_wt[i,]) - t(converted_wt[i-1,])), 0.5)
    up95_winter <- quantile((t(converted_wt[i,]) - t(converted_wt[i-1,])), 0.95)
    
    low95_all <- quantile((t(converted_all[i,]) - t(converted_all[i-1,])), 0.05)
    median_all <- quantile((t(converted_all[i,]) - t(converted_all[i-1,])), 0.5)
    up95_all <- quantile((t(converted_all[i,]) - t(converted_all[i-1,])), 0.95)
  } else{
    low95_spring <- 0
    median_spring <- 0
    up95_spring <- 0
    
    low95_summer <- 0
    median_summer <- 0
    up95_summer <- 0
    
    low95_autumn <- 0
    median_autumn <- 0
    up95_autumn <- 0
    
    low95_winter <- 0
    median_winter <- 0
    up95_winter <- 0
    
    low95_all <- 0
    median_all <- 0
    up95_all <- 0
  }
  
  combined <- data.frame(age_midpoint = age_midpoint,
                         low95_spring, median_spring, up95_spring, 
                         low95_summer, median_summer, up95_summer, 
                         low95_autumn, median_autumn, up95_autumn, 
                         low95_winter, median_winter, up95_winter,
                         low95_all, median_all, up95_all)
  
  incidence[[i]] <- combined
}

incidence_test <- do.call("rbind", incidence)
rownames(incidence_test) <- NULL

incidence_test %>% filter(age_midpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median_spring)) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95_spring, ymax = up95_spring)) +
  labs(title = "Proportion of new seroconversions in spring", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

incidence_test %>% filter(age_midpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median_summer)) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95_summer, ymax = up95_summer)) +
  labs(title = "Proportion of new seroconversions in summer", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

incidence_test %>% filter(age_midpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median_autumn)) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95_autumn, ymax = up95_autumn)) +
  labs(title = "Proportion of new seroconversions in autumn", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

incidence_test %>% filter(age_midpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median_winter)) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95_winter, ymax = up95_winter)) +
  labs(title = "Proportion of new seroconversions in winter", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

incidence_test %>% #filter(age_midpoint <= 365) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median_all)) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95_all, ymax = up95_all)) +
  labs(title = "Proportion of new seroconversions (all)", x = "Age (days)",
       y = "% additional seroconversion") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()

incidence_test_long <- incidence_test %>%
  pivot_longer(
    cols = -age_midpoint,  # keep age as is
    names_to = c("measure", "season"),  # split column names into two parts
    names_sep = "_",  # separator is underscore (e.g. "low95_sp")
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = measure,  # now spread low95, median, up95 into columns
    values_from = value
  ) %>%
  mutate(season = factor(season, levels = c("spring", "summer", "autumn", "winter", "all")))

path_model <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"  #Where model outputs are stored
incidence_test_long %>% write.csv(paste0(path_model, "incidence by age2.csv")) 


plt <- incidence_test_long %>% filter(age_midpoint <= 365 & season != "all") %>%
        ggplot() + 
        geom_point(aes(x = age_midpoint, y = median, col = season), size = 3) +
        geom_errorbar(aes(x = age_midpoint, 
                          ymin = low95, ymax = up95,
                          col = season),width = 10, linewidth = 1) +
        labs(title = "New RSV seroconversions by age", x = "Age (days)",
             y = "% new seroconversions",
             col = "Season of birth") +
        scale_y_continuous(labels = scales::percent) +
        theme_light() +
        facet_wrap(~season, scale = "free") +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))
plt

plt %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Seroconversion by season of birth and age.png",
         width = 20, height = 14, units = "in", 
         device='png')

# Plot this together with true estimates
incidence_comp <- incidence_data_season %>% 
  filter (xMidpoint <= 365) %>%
  mutate(data = "Data") %>%
  rbind(
    incidence_test_long %>% 
      filter(age_midpoint <= 365 & season != "all") %>%
      rename(xMidpoint = age_midpoint) %>%
      mutate(season_birth = season,
             incidence_mean = median,
             incidence_low95 = low95,
             incidence_up95 = up95,
             data = "Model") %>%
      select(c(xMidpoint, season_birth, incidence_mean, incidence_low95, incidence_up95, data))
  )

plt_incidence_comp <- incidence_comp %>%
  mutate(season_birth = factor(season_birth, levels = c("autumn", "winter", "spring", "summer"))) %>%
  ggplot()+
  geom_point(aes(x = xMidpoint, y = incidence_mean, col = data, shape = data), size = 3) +
  geom_errorbar(aes(x = xMidpoint, ymin = incidence_low95, ymax = incidence_up95,
                    col = data, linetype = data), width = 10, linewidth = 1) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% additional seroconversions\n", col = "Data",
       shape = "Data", linetype = "Data") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season_birth) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"),
         legend.text = element_text(size = 20))

plt_incidence_comp
plt_incidence_comp %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Incidence comparison.png",
         width = 20, height = 14, units = "in", 
         device='png')

plt_all <- incidence_test_long %>% 
  filter(season == "all" & age_midpoint > 30) %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median, col = season), size = 5) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95, ymax = up95,
                    col = season), width = 45, linewidth = 1.5) +
  labs(title = "New RSV seroconversions by age", x = "Age (days)",
       y = "% new seroconversions\n") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460, 1825)) +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
         axis.text.y = element_text(size = 30),
         axis.title.x = element_text(size = 45),
         axis.title.y = element_text(size = 45),
         title = element_text(size = 35))
  
plt_all

plt_all %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Seroconversion by age all.png",
         width = 20, height = 14, units = "in", 
         device='png')

plt_all <- incidence_test_long %>% 
  filter(age_midpoint <= 365 & age_midpoint > 30 & season == "all") %>%
  ggplot() + 
  geom_point(aes(x = age_midpoint, y = median, col = season), size = 3) +
  geom_errorbar(aes(x = age_midpoint, 
                    ymin = low95, ymax = up95,
                    col = season),width = 10, linewidth = 1) +
  labs(title = "New RSV seroconversions by age", x = "Age (days)",
       y = "% new seroconversions\n") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25))

plt_all

plt_all %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Seroconversion by age all <1.png",
         width = 22, height = 14, units = "in", 
         device='png')

# --------------------------------- Scaling hospitalisations -------------------
library(readxl)
path_paper <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/" # Where results from other studies are saved
path_model <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"  #Where model outputs are stored

# Read in RSV-illness numbers
# mild_illness <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Mild illness numbers") %>% 
#   janitor::clean_names()
severe_illness <- read_excel(paste0(path_paper,"SA estimates/RSV illness rates SA.xlsx"), sheet = "Pooled severe illness") %>% 
  janitor::clean_names() %>%
  select(age_days, total_severe_illness_rate, total_severe_illness_lower_ci, 
         total_severe_illness_upper_ci)

severe_illness_progression <- incidence_test_long %>%
  merge(severe_illness, by.x = "age_midpoint", by.y = "age_days") %>%
  mutate(cases = case_when(age_midpoint == 30 ~ 0,
                           T ~ (total_severe_illness_rate/100000) / median),
         age_midpoint_text = case_when(
           age_midpoint == 30 ~ "[0-2)",
           age_midpoint == 91 ~ "[2-4)",
           age_midpoint == 152 ~ "[4-6)",
           age_midpoint == 212 ~ "[6-8)",
           age_midpoint == 272 ~ "[8-10)",
           age_midpoint == 332 ~ "[10-12)",
           T ~ as.character(age_midpoint)
         )) %>%
  mutate(age_midpoint_text = factor(age_midpoint_text, 
                                    levels = c("[0-2)", "[2-4)", 
                                               "[4-6)", "[6-8)", 
                                               "[8-10)", "[10-12)"))) %>%
  select(age_midpoint_text, age_midpoint, season, cases)

plt_progression <- severe_illness_progression %>% ggplot() +
  geom_point(aes(x = age_midpoint_text, y = cases/100, col = season), size = 5) +
  facet_wrap(~season) +
  labs(title = "New RSV hospitalisations by age", x = "Age midpoint (months)",
       y = "% new hospitalisations\n",
       colour = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"),
         legend.text = element_text(size = 20))

plt_progression

plt_progression %>%
  ggsave(filename = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/Hosp by age and season.png",
         width = 22, height = 14, units = "in", 
         device='png')

# -----------------
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

seroconversion <- list()

for (i in 1:(5*365)){
  age_midpoint <- i
  
  low95_sp <- quantile(t(converted_sp[,i]), 0.05, na.rm = T)
  median_sp <- quantile(t(converted_sp[,i]), 0.5, na.rm = T)
  up95_sp <- quantile(t(converted_sp[,i]), 0.95, na.rm = T)
  
  low95_sm <- quantile(t(converted_sm[,i]), 0.05, na.rm = T)
  median_sm <- quantile(t(converted_sm[,i]), 0.5, na.rm = T)
  up95_sm <- quantile(t(converted_sm[,i]), 0.95, na.rm = T)
  
  low95_au <- quantile(t(converted_au[,i]), 0.05, na.rm = T)
  median_au <- quantile(t(converted_au[,i]), 0.5, na.rm = T)
  up95_au <- quantile(t(converted_au[,i]), 0.95, na.rm = T)
  
  low95_wt <- quantile(t(converted_wt[,i]), 0.05, na.rm = T)
  median_wt <- quantile(t(converted_wt[,i]), 0.5, na.rm = T)
  up95_wt <- quantile(t(converted_wt[,i]), 0.95, na.rm = T)
  
  low95_all <- quantile(t(converted_all[,i]), 0.05, na.rm = T)
  median_all <- quantile(t(converted_all[,i]), 0.5, na.rm = T)
  up95_all <- quantile(t(converted_all[,i]), 0.95, na.rm = T)
  
  combined <- data.frame(age_midpoint = age_midpoint,
                         low95_sp, median_sp, up95_sp, 
                         low95_sm, median_sm, up95_sm, 
                         low95_au, median_au, up95_au, 
                         low95_wt, median_wt, up95_wt,
                         low95_all, median_all, up95_all)
  
  seroconversion[[i]] <- combined
}

seroconversion_df <- do.call("rbind", seroconversion)
rownames(seroconversion_df) <- NULL

seroconversion_long <- seroconversion_df %>%
  pivot_longer(
    cols = -age_midpoint,  # keep age as is
    names_to = c("measure", "season"),  # split column names into two parts
    names_sep = "_",  # separator is underscore (e.g. "low95_sp")
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = measure,  # now spread low95, median, up95 into columns
    values_from = value
  ) %>%
  mutate(season = case_when(season == "sp" ~ "Spring",
                            season == "sm" ~ "Summer",
                            season == "wt" ~ "Winter",
                            season == "au" ~ "Autumn",
                            season == "all" ~ "All"),
         season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter", "All")))

seroconversion_long %>% ggplot() +
  geom_line(aes(x = age_midpoint, y = median, col = season)) + 
  geom_ribbon(aes(x = age_midpoint, ymin=low95, ymax=up95, fill=season), alpha=.2) +
  labs(title = "Proportion of seroconverted children", x = "Age (days)",
       y = "% seroconverted\n",
       col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))

# Get incident cases
incidence <- list()

for (i in 1:(5*365)){
  age_midpoint <- i
  
  if(i>1){
    low95_spring <- quantile(
      (t(converted_sp[,i]) - t(converted_sp[,i-1])), 0.05)
    median_spring <- quantile(
      (t(converted_sp[,i]) - t(converted_sp[,i-1])), 0.5)
    up95_spring <- quantile(
      (t(converted_sp[,i]) - t(converted_sp[,i-1])), 0.95)
    
    low95_summer <- quantile((t(converted_sm[,i]) - t(converted_sm[,i-1])), 0.05)
    median_summer <- quantile((t(converted_sm[,i]) - t(converted_sm[,i-1])), 0.5)
    up95_summer <- quantile((t(converted_sm[,i]) - t(converted_sm[,i-1])), 0.95)
    
    low95_autumn <- quantile((t(converted_au[,i]) - t(converted_au[,i-1])), 0.05)
    median_autumn <- quantile((t(converted_au[,i]) - t(converted_au[,i-1])), 0.5)
    up95_autumn <- quantile((t(converted_au[,i]) - t(converted_au[,i-1])), 0.95)
    
    low95_winter <- quantile((t(converted_wt[,i]) - t(converted_wt[,i-1])), 0.05)
    median_winter <- quantile((t(converted_wt[,i]) - t(converted_wt[,i-1])), 0.5)
    up95_winter <- quantile((t(converted_wt[,i]) - t(converted_wt[,i-1])), 0.95)
    
    low95_all <- quantile((t(converted_all[,i]) - t(converted_all[,i-1])), 0.05)
    median_all <- quantile((t(converted_all[,i]) - t(converted_all[,i-1])), 0.5)
    up95_all <- quantile((t(converted_all[,i]) - t(converted_all[,i-1])), 0.95)
  } else{
    low95_spring <- 0
    median_spring <- 0
    up95_spring <- 0
    
    low95_summer <- 0
    median_summer <- 0
    up95_summer <- 0
    
    low95_autumn <- 0
    median_autumn <- 0
    up95_autumn <- 0
    
    low95_winter <- 0
    median_winter <- 0
    up95_winter <- 0
    
    low95_all <- 0
    median_all <- 0
    up95_all <- 0
  }
  
  combined <- data.frame(age_midpoint = age_midpoint,
                         low95_spring, median_spring, up95_spring, 
                         low95_summer, median_summer, up95_summer, 
                         low95_autumn, median_autumn, up95_autumn, 
                         low95_winter, median_winter, up95_winter,
                         low95_all, median_all, up95_all)
  
  incidence[[i]] <- combined
}

incidence_test <- do.call("rbind", incidence)
rownames(incidence_test) <- NULL

incidence_test_long <- incidence_test %>%
  pivot_longer(
    cols = -age_midpoint,  # keep age as is
    names_to = c("measure", "season"),  # split column names into two parts
    names_sep = "_",  # separator is underscore (e.g. "low95_sp")
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = measure,  # now spread low95, median, up95 into columns
    values_from = value
  ) %>%
  mutate(season = case_when(season == "spring" ~ "Spring",
                            season == "summer" ~ "Summer",
                            season == "winter" ~ "Winter",
                            season == "autumn" ~ "Autumn",
                            season == "all" ~ "All"),
         season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter", "All")))

incidence_test_long %>% ggplot() +
  geom_line(aes(x = age_midpoint, y = median, col = season)) + 
  geom_ribbon(aes(x = age_midpoint, ymin=low95, ymax=up95, fill=season), alpha=.2) +
  labs(title = "Proportion of new seroconversions", x = "Age (days)",
       y = "% new seroconversions\n",
       col = "Season of birth") +
  scale_y_continuous(labels = scales::percent) +
  theme_light() +
  facet_wrap(~season) +
  theme (axis.ticks.y = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25),
         strip.text.x = element_text(size = 25, color = "black"))


 # -----------------------------------------------------------------------------
# Testing time values
spring_FOI_sp_func <- function(t){
  spring_FOI_sp <- if ((t <= 30.41*1.5) ||  
                     ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                     ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                     ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                     ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                     ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
  return (spring_FOI_sp)
}

summer_FOI_sp_func <- function(t){
  summer_FOI_sp <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      # FOI in summer for those born in spring
                       ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                       ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                       ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                       ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
  return(summer_FOI_sp)
  
}

autumn_FOI_sp_func <- function(t){
  autumn_FOI_sp <- if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                        ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                        ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                        ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                        ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
 return(autumn_FOI_sp)
   
}

winter_FOI_sp_func <- function(t){
  winter_FOI_sp <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                        ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                        ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                        ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                        ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
  
  return(winter_FOI_sp)
  
}

spring_FOI_sm_func <- function (t){
  spring_FOI_sm <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                        ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                        ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                        ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                        ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
  
  return (spring_FOI_sm)
}

summer_FOI_sm_func <- function (t){
  summer_FOI_sm <- if ((t <= 30.41*1.5) ||  
                       ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                       ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                       ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                       ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                       ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
  
 return (summer_FOI_sm) 
}

autumn_FOI_sm_func <- function(t){
  autumn_FOI_sm <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      # FOI in summer for those born in spring
                       ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                       ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                       ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                       ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
 return (autumn_FOI_sm) 
}

winter_FOI_sm_func <- function(t){
  winter_FOI_sm <-  if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                         ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                         ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                         ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                         ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
 return (winter_FOI_sm) 
}

spring_FOI_au_func <- function(t) {
  spring_FOI_au <-  if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                         ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                         ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                         ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                         ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
 return (spring_FOI_au) 
}

summer_FOI_au_func <- function(t){
  summer_FOI_au <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                        ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                        ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                        ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                        ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
 return (summer_FOI_au) 
}

autumn_FOI_au_func <- function(t){
  autumn_FOI_au <- if ((t <= 30.41*1.5) ||  
                       ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                       ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                       ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                       ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                       ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
 return (autumn_FOI_au) 
}

winter_FOI_au_func <- function(t){
  winter_FOI_au <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      # FOI in summer for those born in spring
                       ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                       ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                       ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                       ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
  
 return (winter_FOI_au) 
}

spring_FOI_wt_func <- function(t){
  spring_FOI_wt <- if ((t > 30.41*1.5 && t <= 30.41*4.5 ) ||                      # FOI in summer for those born in spring
                       ((t > 30.41*13.5) && (t <= 30.41*16.5)) ||
                       ((t > 30.41*25.5) && (t <= 30.41*28.5)) ||
                       ((t > 30.41*37.5) && (t <= 30.41*40.5)) ||
                       ((t > 30.41*49.5) && (t <= 30.41*52.5))) 1 else 0
  
 return (spring_FOI_wt) 
}

summer_FOI_wt_func <- function (t){
  summer_FOI_wt <-  if ( (t > 30.41*4.5 && t <= 30.41*7.5 ) ||                    
                         ((t > 30.41*16.5) && (t <= 30.41*19.5)) ||
                         ((t > 30.41*28.5) && (t <= 30.41*31.5)) ||
                         ((t > 30.41*40.5) && (t <= 30.41*43.5)) ||
                         ((t > 30.41*52.5) && (t <= 30.41*55.5))) 1 else 0
 return (summer_FOI_wt) 
}

autumn_FOI_wt_func <- function(t){
  autumn_FOI_wt <- if ( (t > 30.41*7.5 && t <= 30.41*10.5 ) ||                   
                        ((t > 30.41*19.5) && (t <= 30.41*22.5)) ||
                        ((t > 30.41*31.5) && (t <= 30.41*34.5)) ||
                        ((t > 30.41*43.5) && (t <= 30.41*46.5)) ||
                        ((t > 30.41*55.5) && (t <= 30.41*58.5))) 1 else 0
 return (autumn_FOI_wt) 
}

winter_FOI_wt_func <- function (t){
  winter_FOI_wt <- if ((t <= 30.41*1.5) ||  
                       ((t > 30.41*10.5) && (t <= (30.41*13.5))) ||
                       ((t > 30.41*22.5) && (t <= (30.41*25.5))) ||
                       ((t > 30.41*34.5) && (t <= (30.41*37.5))) || 
                       ((t > 30.41*46.5) && (t <= (30.41*49.5))) || 
                       ((t > 30.41*58.5) && (t <= (30.41*61.5)))) 1 else 0
  
}

# Get parameter estimates from run to estimates the FOI for meach birth cohort
spring_comp <- 3
summer_comp <- 1
autumn_comp <- 4
winter_comp <- 5

# Define where to store the values
lambda_sp_list <- list()
lambda_sm_list <- list()
lambda_au_list <- list()
lambda_wt_list <- list()

# Calculate the FOIs 
for (t in seq_along(1:(365*5))){
  spring_FOI_sp <- spring_FOI_sp_func(t)
  summer_FOI_sp <- summer_FOI_sp_func(t)
  autumn_FOI_sp <- autumn_FOI_sp_func(t)
  winter_FOI_sp <- winter_FOI_sp_func(t)
  
  spring_FOI_sm <- spring_FOI_sm_func(t)
  summer_FOI_sm <- summer_FOI_sm_func(t)
  autumn_FOI_sm <- autumn_FOI_sm_func(t)
  winter_FOI_sm <- winter_FOI_sm_func(t)
  
  spring_FOI_au <- spring_FOI_au_func(t)
  summer_FOI_au <- summer_FOI_au_func(t)
  autumn_FOI_au <- autumn_FOI_au_func(t)
  winter_FOI_au <- winter_FOI_au_func(t)
  
  spring_FOI_wt <- spring_FOI_wt_func(t)
  summer_FOI_wt <- summer_FOI_wt_func(t)
  autumn_FOI_wt <- autumn_FOI_wt_func(t)
  winter_FOI_wt <- winter_FOI_wt_func(t)
  
  # Putting it all together into four FOIs
  lambda_sp = (summer_comp + spring_comp) * spring_FOI_sp + 
    summer_comp * summer_FOI_sp + 
    (summer_comp + autumn_comp) * autumn_FOI_sp +
    (summer_comp + winter_comp) * winter_FOI_sp 
  
  lambda_sm = (summer_comp + spring_comp) * spring_FOI_sm + 
    summer_comp * summer_FOI_sm + 
    (summer_comp + autumn_comp) * autumn_FOI_sm + 
    (summer_comp + winter_comp) * winter_FOI_sm 
  
  lambda_au = (summer_comp + spring_comp) * spring_FOI_au + 
    summer_comp * summer_FOI_au + 
    (summer_comp + autumn_comp) * autumn_FOI_au + 
    (summer_comp + winter_comp) * winter_FOI_au 
  
  lambda_wt = (summer_comp + spring_comp) * spring_FOI_wt + 
    summer_comp * summer_FOI_wt +
    (summer_comp + autumn_comp) * autumn_FOI_wt + 
    (summer_comp + winter_comp) * winter_FOI_wt 
  
  lambda_sp_list[[t]] <- as.data.frame(lambda_sp)
  lambda_sp_list[[t]]$time <- t
  lambda_sm_list[[t]] <- as.data.frame(lambda_sm)
  lambda_sm_list[[t]]$time <- t
  lambda_au_list[[t]] <- as.data.frame(lambda_au)
  lambda_au_list[[t]]$time <- t
  lambda_wt_list[[t]] <- as.data.frame(lambda_wt)
  lambda_wt_list[[t]]$time <- t
  
}

# Combine into 4 dataframes
lambda_sp_df <- do.call("rbind", lambda_sp_list)
lambda_sm_df <- do.call("rbind", lambda_sm_list)
lambda_au_df <- do.call("rbind", lambda_au_list)
lambda_wt_df <- do.call("rbind", lambda_wt_list)

# Combine them into one dataframe
lambda_values <- lambda_sp_df %>%
  merge(lambda_sm_df, by = "time") %>%
  merge(lambda_au_df, by = "time") %>%
  merge(lambda_wt_df, by = "time") 

# Pivot longer
lambda_values_long <- lambda_values %>%
  pivot_longer(
    cols = -time,  # keep time as is
    names_to = "season_birth",  # split column names into two parts
    values_to = "FOI"
  ) %>%
  mutate(season_birth = case_when(season_birth == "lambda_sp" ~ "Spring",
                            season_birth == "lambda_sm" ~ "Summer",
                            season_birth == "lambda_au" ~ "Autumn",
                            season_birth == "lambda_wt" ~ "Winter"),
         season_birth = factor(season_birth, levels = c("Spring", "Summer", "Autumn", "Winter")))

# Plot FOI over time
lambda_values_long %>%# filter(time <= 365) %>%
  ggplot() +
  geom_line(aes(x = time/30.41, y = FOI, col = season_birth)) +
  labs(title = "Force of infection by age and season of birth", x = "Age (months)",
       y = "Force of infection\n",
       col = "Season of birth") +
  theme_light() +
  theme (axis.ticks.y = element_blank(),
         legend.position = "right",
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
         axis.text.y = element_text(size = 20),
         axis.title.x = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         title = element_text(size = 25)) +
  facet_wrap(~season_birth)


