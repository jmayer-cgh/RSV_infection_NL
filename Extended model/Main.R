# Fitting the odin model to the data

# Housekeeping
rm(list=ls())

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
          prop_seroconv = n_infection/N,
          cum_infection = prop_seroconv * N_tot) %>%
  ungroup() %>% 
  distinct()

# Calculate seroprevalence and binomial confidence intervals
incidence_data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom::binom.confint(incidence_data$n_infection, 
                                                                                             incidence_data$N,
                                                                                             method="exact")[,c("mean","lower","upper")]

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

incidence_data_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom::binom.confint(incidence_data_season$n_infection, 
                                                                                                    incidence_data_season$N, 
                                                                                                    method="exact")[,c("mean","lower","upper")]

incidence_data_season_wide <- incidence_data_season %>% 
  select (!c(age_mid, age_grp, seroprev_mean)) %>%
  pivot_wider(
    names_from = season_birth,
    values_from = c(N , n_infection, prop_seroconv, seroprev_low95, seroprev_up95, cum_infection),
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

# ----------- Diagnostics -----------------------------------------------------
# Check mixing
matplot(samples_tuned$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")

# Thin and check mixing
samples_thinned <- monty_samples_thin(samples_tuned, burnin = 16000)
matplot(samples_thinned$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Trace log posterior season mu adaptive.png");
dev.off ()

# Check trace
draws_thinned <- as_draws_df(samples_thinned)
mcmc_trace(draws_thinned, facet_args = list(ncol = 1, strip.position = "left"))
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Trace parameters season mu adaptive.png");
dev.off ()

# Summary
params_est <- summarise_draws(draws_thinned)

draws_array_tuned <- as_draws_array(draws_thinned)
# Posterior uncertainty intervals
# Plot them separately because the values are very different
mcmc_intervals(draws_array_tuned, pars = c("spring_comp", "summer_comp", "autumn_comp", "winter_comp"))
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Uncertainty FOI season mu adaptive.png");
dev.off ()

mcmc_intervals(draws_array_tuned, pars = "mu")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Uncertainty mu season mu adaptive.png");
dev.off ()

mcmc_intervals(draws_array_tuned, pars = "prop")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Uncertainty pi season mu adaptive.png");
dev.off ()

# Univariate marginal posterior distributions
mcmc_hist(draws_array_tuned)
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Marginal posterior distributions season mu adaptive.png");
dev.off ()

# Trajectories
trajectories <- dust_unpack_state(filter,
                                  samples_thinned$observations$trajectories)
sero_conv <- array(trajectories$R_all, c(19, 1000))

# Plot the fit
matplot(incidence_data$time, sero_conv, type = "l", col = "#00000044", lty = 1,
        xlab = "Age", ylab = "Sero-conversion", main = "Model fit")
points(x = incidence_data$time, y = incidence_data$prop_seroconv, pch = 19, col = "red")
arrows(x0 = incidence_data$time, y0 = incidence_data$seroprev_low95, 
       x1 = incidence_data$time, y1 = incidence_data$seroprev_up95, 
       angle = 90, code = 3, length = 0.1, col = "red")

dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Model fit season mu adaptive.png");
dev.off ()

# By season
# Spring
sero_conv_spring <- array(trajectories$R_sp, c(19, 1000))
matplot(incidence_data_season_wide$time, sero_conv_spring, type = "l", col = "#00000044", lty = 1,
        xlab = "Age", ylab = "Sero-conversion", main = "Model fit - spring")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_spring, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_spring, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_spring, 
       angle = 90, code = 3, length = 0.1, col = "red")

dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Model fit season spring mu adaptive.png");
dev.off ()

# Summer
sero_conv_summer <- array(trajectories$R_sm, c(19, 1000))
matplot(incidence_data_season_wide$time, sero_conv_summer, type = "l", col = "#00000044", lty = 1,
        xlab = "Age", ylab = "Sero-conversion", main = "Model fit - summer")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_summer, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_summer, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_summer, 
       angle = 90, code = 3, length = 0.1, col = "red")

dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Model fit season summer mu adaptive.png");
dev.off ()

# Autumn
sero_conv_autumn <- array(trajectories$R_au, c(19, 1000))
matplot(incidence_data_season_wide$time, sero_conv_autumn, type = "l", col = "#00000044", lty = 1,
        xlab = "Age", ylab = "Sero-conversion", main = "Model fit - autumn")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_autumn, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_autumn, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_autumn, 
       angle = 90, code = 3, length = 0.1, col = "red")

dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Model fit season autumn adaptive.png");
dev.off ()

# Winter
sero_conv_winter <- array(trajectories$R_wt, c(19, 1000))
matplot(incidence_data_season_wide$time, sero_conv_winter, type = "l", col = "#00000044", lty = 1,
        xlab = "Age", ylab = "Sero-conversion", main = "Model fit - winter")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_winter, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_winter, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_winter, 
       angle = 90, code = 3, length = 0.1, col = "red")

dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Model fit season winter adaptive.png");
dev.off ()

# ----------- Estimates --------------------------------------------------------
hpd <- apply(pmcmc_tuned_run$pars, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm=T) 
rownames(hpd) <- c("low95", "median", "up95")

# Get the mean highest probability distribution
mean_hpd <- apply(pmcmc_tuned_run$pars, 2, mean)
lower_95_hpd <- apply(pmcmc_tuned_run$pars, 2, quantile, probs = 0.05)
upper_95_hpd <- apply(pmcmc_tuned_run$pars, 2, quantile, probs = 0.95)
summary_hpd <- rbind(mean_hpd, lower_95_hpd, upper_95_hpd)
summary_hpd

# Seroconversion estimates
seroconversion <- list()

for (i in 1:length(incidence_data_season_wide$time)){
  age_midpoint <- incidence_data_season_wide$time[i]
  
  low95_sp <- quantile(trajectories$R_sp[i, , ], 0.05)
  median_sp <- quantile(trajectories$R_sp[i, , ], 0.5)
  up95_sp <- quantile(trajectories$R_sp[i, , ], 0.95)
  
  low95_sm <- quantile(trajectories$R_sm[i, , ], 0.05)
  median_sm <- quantile(trajectories$R_sm[i, , ], 0.5)
  up95_sm <- quantile(trajectories$R_sm[i, , ], 0.95)
  
  low95_au <- quantile(trajectories$R_au[i, , ], 0.05)
  median_au <- quantile(trajectories$R_au[i, , ], 0.5)
  up95_au <- quantile(trajectories$R_au[i, , ], 0.95)
  
  low95_wt <- quantile(trajectories$R_wt[i, , ], 0.05)
  median_wt <- quantile(trajectories$R_wt[i, , ], 0.5)
  up95_wt <- quantile(trajectories$R_wt[i, , ], 0.95)
  
  combined <- data.frame(age_midpoint = age_midpoint,
                         low95_sp, median_sp, up95_sp, 
                         low95_sm, median_sm, up95_sm, 
                         low95_au, median_au, up95_au, 
                         low95_wt, median_wt, up95_wt)
  
  seroconversion[[i]] <- combined
}

seroconversion_df <- do.call("rbind", seroconversion)
rownames(seroconversion_df) <- NULL

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
  )

# Save files
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"

write.csv(params_est, paste0(path, "highest probability distribution.csv"), row.names = F)
write.csv(seroconversion_df, paste0(path, "seroconversion by age.csv"), row.names = F)
write.csv(incidence_df, paste0(path, "incidence by age.csv"), row.names = F)
