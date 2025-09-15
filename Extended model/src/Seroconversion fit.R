# Fitting the odin model to the data

# Load libraries
library(dust2)
library(odin2)
library(monty)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(bayesplot)

# ----------- Model ------------------------------------------------------------
# Load model
msr <- odin("~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Extended model/Odin2 model.R")

# ----------- Data -------------------------------------------------------------
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

# ----------- Fitting ----------------------------------------------------------
# Build a filter and test it
filter <- dust_filter_create(msr, time_start = 0, data = incidence_data_season_wide, n_particles = 1000)

# Now running an MCMC
# List of inputs
packer <- monty_packer(c("spring_comp", "summer_comp", "autumn_comp", "winter_comp", 
                         "mu", "prop"))

# Likelihood
likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = T)
dust_likelihood_run(filter, list(spring_comp = 0.003478,
                                 summer_comp = 0.00060,
                                 autumn_comp = 0.00264,
                                 winter_comp = 0.00692, 
                                 mu = 0.09,
                                 prop = 0.3),
                    save_trajectories = T)

# Priors
prior <- monty_dsl({
  spring_comp ~ Uniform(0, 0.1)
  summer_comp ~ Uniform(0, 0.1)
  autumn_comp ~ Uniform(0, 0.1)
  winter_comp ~ Uniform(0, 0.1)
  mu ~ Uniform(0.0001, 0.1)
  prop ~ Uniform(0, 1)
})

# Posterior
posterior <- likelihood + prior

# Define a sampler (adaptive MCMC)
vcv <- diag(6) * 0.0004
sampler <- monty_sampler_adaptive(vcv)

# Parallelise the process
runner <- monty_runner_callr(2, progress = T)

# Sample
properties <- monty_model_properties(is_stochastic = F)
model <- monty_model(posterior, properties = properties)
samples <- monty_sample(model, sampler, 80000, 
                        initial = c(1e-05, 0.02002, 3e-05, 4e-05, 0.09, 0.5),
                        n_chains = 2,
                        runner = runner)

# Tune the sampler
draws <- as_draws_df(samples)
vcv_tuned <- cov(draws[1:6])

# Sample
runner <- monty_runner_callr(6, progress = T)
sampler_tuned <- monty_sampler_adaptive(vcv_tuned)
samples_tuned <- monty_sample(model, sampler_tuned, 80000, 
                              initial = c(1e-05, 0.02002, 3e-05, 4e-05, 0.09, 0.5),
                              runner = runner,
                              n_chains = 6)

# Thin
samples_thinned <- monty_samples_thin(samples_tuned, burnin = 16000)

# Trajectories
trajectories <- dust_unpack_state(filter,
                                  samples_thinned$observations$trajectories) # 17 states, 19 time points, 80000 - 16000 = 64000 samples, 6 chains

# Extract seroconversion numbers by season from all samples of one chain
converted_all <- array(trajectories$R_all, c(19, 64000))
converted_sp <-  array(trajectories$R_sp, c(19, 64000))
converted_sm <-  array(trajectories$R_sm, c(19, 64000))
converted_au <-  array(trajectories$R_au, c(19, 64000))
converted_wt <-  array(trajectories$R_wt, c(19, 64000))

# Save the outputs in case we don't want to run the whole model again
converted_summary <- list(converted_all, converted_sp, converted_sm, converted_au, converted_wt)
saveRDS(converted_summary, file = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/RDS files/seroconversion model outputs.rds")