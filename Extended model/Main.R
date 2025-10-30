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
library(magrittr)

# ----------- Model ------------------------------------------------------------
# Load model
msr <- odin("~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Extended model/Odin2 model.R")

# ----------- Data -------------------------------------------------------------
# Read in and the data and put it in the right format
data <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/refs/heads/master/data/infection_status.csv")

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
# hpd <- apply(pmcmc_tuned_run$pars, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm=T) 
hpd <- apply(samples_tuned$pars, 2, quantile, prob = c(0.025, 0.5, 0.975), na.rm=T) 
rownames(hpd) <- c("low95", "median", "up95")

# Get the mean highest probability distribution
mean_hpd <- apply(samples_tuned$pars, 2, mean)
lower_95_hpd <- apply(samples_tuned$pars, 2, quantile, probs = 0.025)
upper_95_hpd <- apply(samples_tuned$pars, 2, quantile, probs = 0.975)
summary_hpd <- rbind(mean_hpd, lower_95_hpd, upper_95_hpd)
summary_hpd

params_est <- summarise_draws(draws_thinned, 
                              mean, sd, median, mad, 
                              q5 = ~quantile(.x, 0.025), 
                              q95 = ~quantile(.x , 0.975),
                              rhat, ess_bulk, ess_tail,
                              min, max)
params_est %>% write.csv("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/Parameter estimates.csv")


# Seroconversion estimates
seroconversion <- list()

for (i in 1:length(incidence_data_season_wide$time)){
  age_midpoint <- incidence_data_season_wide$time[i]
  
  low95_sp <- quantile(trajectories$R_sp[i, , ], 0.025)
  median_sp <- quantile(trajectories$R_sp[i, , ], 0.5)
  up95_sp <- quantile(trajectories$R_sp[i, , ], 0.975)
  
  low95_sm <- quantile(trajectories$R_sm[i, , ], 0.025)
  median_sm <- quantile(trajectories$R_sm[i, , ], 0.5)
  up95_sm <- quantile(trajectories$R_sm[i, , ], 0.975)
  
  low95_au <- quantile(trajectories$R_au[i, , ], 0.025)
  median_au <- quantile(trajectories$R_au[i, , ], 0.5)
  up95_au <- quantile(trajectories$R_au[i, , ], 0.975)
  
  low95_wt <- quantile(trajectories$R_wt[i, , ], 0.025)
  median_wt <- quantile(trajectories$R_wt[i, , ], 0.5)
  up95_wt <- quantile(trajectories$R_wt[i, , ], 0.975)
  
  low95_all <- quantile(trajectories$R_all[i, , ], 0.025)
  median_all <- quantile(trajectories$R_all[i, , ], 0.5)
  up95_all <- quantile(trajectories$R_all[i, , ], 0.975)
  
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

# Get incident cases
incidence <- list()

for (i in 1:length(incidence_data_season_wide$time)){
  age_midpoint <- incidence_data_season_wide$time[i]
  
  if(i>1){
    low95_sp <- quantile(trajectories$R_sp[i, , ] - trajectories$R_sp[i-1, , ], 0.025)
    median_sp <- quantile(trajectories$R_sp[i, , ] - trajectories$R_sp[i-1, , ], 0.5)
    up95_sp <- quantile(trajectories$R_sp[i, , ] - trajectories$R_sp[i-1, , ], 0.975)
    
    low95_sm <- quantile(trajectories$R_sm[i, , ] - trajectories$R_sm[i-1, , ], 0.025)
    median_sm <- quantile(trajectories$R_sm[i, , ] - trajectories$R_sm[i-1, , ], 0.5)
    up95_sm <- quantile(trajectories$R_sm[i, , ] - trajectories$R_sm[i-1, , ], 0.975)
    
    low95_au <- quantile(trajectories$R_au[i, , ] - trajectories$R_au[i-1, , ], 0.025)
    median_au <- quantile(trajectories$R_au[i, , ] - trajectories$R_au[i-1, , ], 0.5)
    up95_au <- quantile(trajectories$R_au[i, , ] - trajectories$R_au[i-1, , ], 0.975)
    
    low95_wt <- quantile(trajectories$R_wt[i, , ] - trajectories$R_wt[i-1, , ], 0.025)
    median_wt <- quantile(trajectories$R_wt[i, , ] - trajectories$R_wt[i-1, , ], 0.5)
    up95_wt <- quantile(trajectories$R_wt[i, , ] - trajectories$R_wt[i-1, , ], 0.975)
    
    low95_all <- quantile(trajectories$R_all[i, , ] - trajectories$R_all[i-1, , ], 0.025)
    median_all <- quantile(trajectories$R_all[i, , ] - trajectories$R_all[i-1, , ], 0.5)
    up95_all <- quantile(trajectories$R_all[i, , ] - trajectories$R_all[i-1, , ], 0.975)
  } else{
    low95_sp <- 0
    median_sp <- 0
    up95_sp <- 0
    
    low95_sm <- 0
    median_sm <- 0
    up95_sm <- 0
    
    low95_au <- 0
    median_au <- 0
    up95_au <- 0
    
    low95_wt <- 0
    median_wt <- 0
    up95_wt <- 0
    
    low95_all <- 0
    median_all <- 0
    up95_all <- 0
  }
  
  combined <- data.frame(age_midpoint = age_midpoint,
                         low95_sp, median_sp, up95_sp, 
                         low95_sm, median_sm, up95_sm, 
                         low95_au, median_au, up95_au, 
                         low95_wt, median_wt, up95_wt,
                         low95_all, median_all, up95_all)
  
  incidence[[i]] <- combined
}

incidence_test <- do.call("rbind", incidence)
rownames(incidence_test) <- NULL

# Save files
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"

write.csv(params_est, paste0(path, "highest probability distribution.csv"), row.names = F)
write.csv(seroconversion_df, paste0(path, "seroconversion by age.csv"), row.names = F)
write.csv(incidence_df, paste0(path, "incidence by age.csv"), row.names = F)

incidence_test %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_sp)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_sp, ymax = up95_sp)) +
  labs(title = "Proportion of newly seroconverted children (spring)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Incidence spring.png");
dev.off ()

incidence_test %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_sm)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_sm, ymax = up95_sm)) +
  labs(title = "Proportion of newly seroconverted children (summer)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Incidence summer.png");
dev.off ()

incidence_test %>% ggplot() +
  geom_point(aes(x = age_midpoint, y = median_au)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_au, ymax = up95_au)) +
  labs(title = "Proportion of newly seroconverted children (autumn)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Incidence autumn.png");
dev.off ()

incidence_test %>% ggplot()+
  geom_point(aes(x = age_midpoint, y = median_wt)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_wt, ymax = up95_wt)) +
  labs(title = "Proportion of newly seroconverted children (winter)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Incidence winter.png");
dev.off ()


incidence_test %>% filter(age_midpoint <= 365) %>%
  ggplot()+
  geom_point(aes(x = age_midpoint, y = median_all)) +
  geom_errorbar(aes(x = age_midpoint, ymin = low95_all, ymax = up95_all)) +
  labs(title = "Proportion of seroconverted children (all)", x = "Age (days)",
       y = "% seroconverted") +
  scale_y_continuous(labels = scales::percent) +
  theme_light()
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Incidence all.png");
dev.off ()

# ------------------------------------------------------------------------------
# Simulate the model with the parameter estimates
msr_sim <- odin("~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Extended model/Simulation model.R")


# Get mean and sd for each parameter
mu_mean <- params_est %$% mean[params_est$variable == "mu"]
mu_sd <- params_est %$% sd[params_est$variable == "mu"]
prop_mean <- params_est %$% mean[params_est$variable == "prop"]
prop_sd <- params_est %$% sd[params_est$variable == "prop"]
spring_comp_mean <- params_est %$% mean[params_est$variable == "spring_comp"]
spring_comp_sd <- params_est %$% sd[params_est$variable == "spring_comp"]
summer_comp_mean <- params_est %$% mean[params_est$variable == "summer_comp"]
summer_comp_sd <- params_est %$% sd[params_est$variable == "summer_comp"]
autumn_comp_mean <- params_est %$% mean[params_est$variable == "autumn_comp"]
autumn_comp_sd <- params_est %$% sd[params_est$variable == "autumn_comp"]
winter_comp_mean <- params_est %$% mean[params_est$variable == "winter_comp"]
winter_comp_sd <- params_est %$% sd[params_est$variable == "winter_comp"]

pars <- list(spring_comp_mean = spring_comp_mean,
             spring_comp_sd = spring_comp_sd,
             summer_comp_mean = summer_comp_mean,
             summer_comp_sd = summer_comp_sd,
             autumn_comp_mean = autumn_comp_mean,
             autumn_comp_sd = autumn_comp_sd,
             winter_comp_mean = winter_comp_mean,
             winter_comp_sd = winter_comp_sd,
             mu_mean = mu_mean,
             mu_sd = mu_sd,
             prop_mean = prop_mean,
             prop_sd = prop_sd)

sys <- dust_system_create(msr_sim(), pars, n_particles = 10000, dt = 1)
dust_system_set_state_initial(sys)
time <- 0:(5*365)
y <- dust_system_simulate(sys, time)

# Save the outputs in case we don't want to run the whole model again
prevalence <- dust_unpack_state(sys, y)$R_all
prevalence_sp <- dust_unpack_state(sys, y)$R_sp
prevalence_sm <- dust_unpack_state(sys, y)$R_sm
prevalence_au <- dust_unpack_state(sys, y)$R_au
prevalence_wt <- dust_unpack_state(sys, y)$R_wt
prevalence_summary <- list(prevalence, prevalence_sp, prevalence_sm, prevalence_au, prevalence_wt)
saveRDS(prevalence_summary, file = "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/RDS files/seroconversion simulation outputs.rds")

# Plot the results
# All cohorts
matplot(time, t(prevalence), type = "l", lty = 1, col = "#00000055",
        xlab = "Time (days)", ylab = "% seroconverted", las = 1)
points(prop_seroconv ~ time, incidence_data, pch = 19, col = "red")
arrows(x0 = incidence_data$time, y0 = incidence_data$seroprev_low95, 
       x1 = incidence_data$time, y1 = incidence_data$seroprev_up95, 
       angle = 90, code = 3, length = 0.1, col = "red")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Simulated fit all.png")
dev.off ()

# Spring
matplot(time, t(prevalence_sp), type = "l", lty = 1, col = "#00000055",
        xlab = "Time (days)", ylab = "% seroconverted (spring)", las = 1)
points(prop_seroconv_spring ~ time, incidence_data_season_wide, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_spring, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_spring, 
       angle = 90, code = 3, length = 0.1, col = "red")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Simulated fit spring.png")
dev.off ()

# Summer
matplot(time, t(prevalence_sm), type = "l", lty = 1, col = "#00000055",
        xlab = "Time (days)", ylab = "% seroconverted (summer)", las = 1)
points(prop_seroconv_summer ~ time, incidence_data_season_wide, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_summer, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_summer, 
       angle = 90, code = 3, length = 0.1, col = "red")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Simulated fit summer.png")
dev.off ()

# Autumn
matplot(time, t(prevalence_au), type = "l", lty = 1, col = "#00000055",
        xlab = "Time (days)", ylab = "% seroconverted (autumn)", las = 1)
points(prop_seroconv_autumn ~ time, incidence_data_season_wide, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_autumn, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_autumn, 
       angle = 90, code = 3, length = 0.1, col = "red")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Simulated fit autumn.png")
dev.off ()

# Winter
matplot(time, t(prevalence_wt), type = "l", lty = 1, col = "#00000055",
        xlab = "Time (days)", ylab = "% seroconverted (winter)", las = 1)
points(prop_seroconv_winter ~ time, incidence_data_season_wide, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_winter, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_winter, 
       angle = 90, code = 3, length = 0.1, col = "red")
dev.copy(jpeg,filename="/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Plots/Checks/monty/Simulated fit winter.png")
dev.off ()