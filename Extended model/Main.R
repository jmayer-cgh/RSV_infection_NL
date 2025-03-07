# Fitting the odin model to the data

# Housekeeping
rm(list=ls())

library (odin.dust)
library (odin)
library (dplyr)
library (tidyr)
library (ggplot2)
library (mcstate)

# ----------- Model ------------------------------------------------------------
# Load model
gen_msr <- odin_dust("~/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Extended model/Odin model.R")

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

# Get number of cases by age
incidence_data <- data %>% select (age_grp, age_days, infection) %>%
  mutate(N_tot = n()) %>%
  group_by(age_grp) %>%
  summarise(age_mid = round(median(age_days)), 
            N = n(),
            n_infection = sum(infection),
            prop_seroconv = n_infection/N,
            cum_infection = prop_seroconv * N_tot) %>%
  ungroup() %>% 
  distinct() %>%
  mutate(incidence = c(0,diff(cum_infection))) %>%
  mutate(incidence = case_when(incidence >= 0 ~ round(incidence),
                               TRUE ~ 0))

# Put it in the right format
model_data <- particle_filter_data(data = incidence_data,
                                   initial_time = 0,
                                   time = "age_mid",
                                   rate = 1)

# ----------- Fitting ----------------------------------------------------------
# Comparison function
case_compare <- function(state, observed, pars = NULL){
  nconv <- observed$n_infection # n seroconverted at each age point  (data)
  N <- observed$N # total N at each age point  (data)
  prob <- state[18, , drop = T] #/ colSums(state[2:17, , drop = T]) # proportion seroconverted at each age point (model output) - all seroconverted divided by total populatoion
  # prob <- observed$prop_seroconv
  
  dbinom(x = nconv,
         size = N,
         prob = prob,
         log = TRUE)
}

# Set up a particle filter
n_particles <- 100
filter <- particle_filter$new(data = model_data,
                              model = gen_msr, 
                              n_particles = n_particles,
                              compare = case_compare,
                              seed = 1L)

# Use MCMC to infer parameter values
# Define parameters first
spring_comp <- pmcmc_parameter("spring_comp", initial = 0.000001, max = 0.1, min = 0, prior = function(p)
  dunif(p, max = 0.1, min = 0))
summer_comp <- pmcmc_parameter("summer_comp", initial = 0.000202, max = 0.1, min = 0, prior = function(p)
  dunif(p, max = 0.1, min = 0))
autumn_comp <- pmcmc_parameter("autumn_comp", initial = 0.000001, max = 0.1, min = 0, prior = function(p)
  dunif(p, max = 0.1, min = 0))
winter_comp <- pmcmc_parameter("winter_comp", initial = 0.000001, max = 0.1, min = 0, prior = function(p)
  dunif(p, max = 0.1, min = 0))
pi <- pmcmc_parameter("pi", initial = 0.5, max = 1, min = 0, prior = function(p)
  dunif(p, max = 1, min = 0))

proposal_matrix <- diag(0.1, 5)
mcmc_pars <- pmcmc_parameters$new(list(spring_comp = spring_comp, 
                                       summer_comp = summer_comp,
                                       autumn_comp = autumn_comp,
                                       winter_comp = winter_comp,
                                       pi = pi
                                       ),
                                  proposal_matrix)
# Chain parameters
n_steps <- 1e4
n_burnin <- 500
control <- pmcmc_control(
  n_steps, 
  save_state = T,
  save_trajectories = T,
  progress = T)

# Run the MCMC once
pmcmc_run <-pmcmc(mcmc_pars, filter, control = control)

# Tune the MCMC
proposal_matrix <- cov(pmcmc_run$pars)
mcmc_pars <- mcstate::pmcmc_parameters$new(
  list(spring_comp = spring_comp, 
       summer_comp = summer_comp,
       autumn_comp = autumn_comp,
       winter_comp = winter_comp, 
       pi = pi
       ),
  proposal_matrix)

# Run the tuned filter
control <- pmcmc_control(
  n_steps,
  save_state = T,
  save_trajectories = T,
  progress = T,
  n_chains = 4)

pmcmc_tuned_run <- pmcmc(mcmc_pars, filter, control = control)

mcmc2 <- coda::as.mcmc(cbind(
  pmcmc_tuned_run$probabilities, pmcmc_tuned_run$pars))

# ----------- Diagnostics -----------------------------------------------------
summary(mcmc2) # summary
coda::effectiveSize(mcmc2) # ESS
plot(mcmc2) # plot the trace


1 - coda::rejectionRate(mcmc2) # check the acceptance rate

# Check if the chains converged using a sample from the posterior distribution
mcmc_sample <- pmcmc_sample(pmcmc_tuned_run, n_sample = 500)

par(mfrow = c(4, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
plot(pmcmc_tuned_run$pars[, "winter_comp"], type = "l", xlab = "Iteration",
     ylab = "winter_comp")
plot(pmcmc_tuned_run$pars[, "spring_comp"], type = "l", xlab = "Iteration",
     ylab = "spring_comp")
plot(pmcmc_tuned_run$pars[, "summer_comp"], type = "l", xlab = "Iteration",
     ylab = "summer_comp")
plot(pmcmc_tuned_run$pars[, "autumn_comp"], type = "l", xlab = "Iteration",
     ylab = "autumn_comp")
plot(pmcmc_tuned_run$pars[, "pi"], type = "l", xlab = "Iteration",
     ylab = "pi")
hist(mcmc_sample$pars[, "winter_comp"], main = "", xlab = "winter_comp",
     freq = FALSE)
hist(mcmc_sample$pars[, "spring_comp"], main = "", xlab = "spring_comp",
     freq = FALSE)
hist(mcmc_sample$pars[, "summer_comp"], main = "", xlab = "summer_comp",
     freq = FALSE)
hist(mcmc_sample$pars[, "autumn_comp"], main = "", xlab = "autumn_comp",
     freq = FALSE)
hist(mcmc_sample$pars[, "pi"], main = "", xlab = "pi",
     freq = FALSE)

# Compare with the data
times <- c(0, incidence_data$age_mid)
cols <- c(`M1 spring` = "#8c8cd9", `M2 spring` = "turquoise", `S spring` = "#cc0044", `R spring` = "#999966", `Total prop \nseroconverted` = "pink")

state <- mcmc_sample$trajectories$state

prop_seroconv <-  t(state[18, , -1])

par(mfrow = c(1, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0), bty = "n")
matplot(times, t(state[2, , ]), type = "l", lty = 1, col = cols["M1 spring"],
        ylim = c(0, max(t(state[18, , ])) * 1.2), xlab = "Day", ylab = "Proportion of individuals")
matlines(times, t(state[3, , ]), lty = 1, col = cols["M2 spring"])
matlines(times, t(state[4, , ]), lty = 1, col = cols["S spring"])
matlines(times, t(state[5, , ]), lty = 1, col = cols["R spring"])
matlines(times, t(state[18, , ]), lty = 1, col = cols["Total prop \nseroconverted"])
legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")

matplot(model_data$age_mid_start, prop_seroconv, type = "l", lty = 1, col = grey(0.7, 0.2),
        xlab = "Day", ylab = "Proportion of seroconverted children", ylim = c(0, 1.2))
points(model_data$age_mid_start, model_data$prop_seroconv, pch = 20)

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

for (i in 1:length(pmcmc_tuned_run$trajectories$time)){
  age_midpoint <- pmcmc_tuned_run$trajectories$time[i]
  
  low95_sp <- quantile(pmcmc_tuned_run$trajectories$state[5, , i], 0.05)
  median_sp <- quantile(pmcmc_tuned_run$trajectories$state[5, , i], 0.5)
  up95_sp <- quantile(pmcmc_tuned_run$trajectories$state[5, , i], 0.95)
  
  low95_sm <- quantile(pmcmc_tuned_run$trajectories$state[9, , i], 0.05)
  median_sm <- quantile(pmcmc_tuned_run$trajectories$state[9, , i], 0.5)
  up95_sm <- quantile(pmcmc_tuned_run$trajectories$state[9, , i], 0.95)
  
  low95_au <- quantile(pmcmc_tuned_run$trajectories$state[13, , i], 0.05)
  median_au <- quantile(pmcmc_tuned_run$trajectories$state[13, , i], 0.5)
  up95_au <- quantile(pmcmc_tuned_run$trajectories$state[13, , i], 0.95)
  
  low95_wt <- quantile(pmcmc_tuned_run$trajectories$state[17, , i], 0.05)
  median_wt <- quantile(pmcmc_tuned_run$trajectories$state[17, , i], 0.5)
  up95_wt <- quantile(pmcmc_tuned_run$trajectories$state[17, , i], 0.95)
  
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
                           incidence = c(0, diff(seroconversion_df$median_sp)),
                           season_birth = "spring") %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence = c(0, diff(seroconversion_df$median_sm)),
               season_birth = "summer")
  ) %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence = c(0, diff(seroconversion_df$median_au)),
               season_birth = "autumn")
  ) %>%
  rbind(
    data.frame(age_midpoint = seroconversion_df$age_midpoint,
               incidence = c(0, diff(seroconversion_df$median_wt)),
               season_birth = "winter")
  )

# Save files
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/mcstate/"

write.csv(summary_hpd, paste0(path, "highest probability distribution, assumed mu.csv"), row.names = F)
write.csv(seroconversion_df, paste0(path, "seroconversion by age, assumed mu.csv"), row.names = F)
write.csv(incidence_df, paste0(path, "incidence by age, assumed mu.csv"), row.names = F)

