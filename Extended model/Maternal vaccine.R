################################################################
# RSV seroconversion MSc project
# Adding daycare attendance as a multiplication
# Author: Julia Mayer
# Last updated:5 December 2024
################################################################


# HOUSEKEEPING ------------------------------------------------------------

rm(list=ls())
set.seed(42)

library(tidyverse)
library(BayesianTools)
library(binom)
library(plyr)
library(deSolve)
library(coda)
library(lubridate)
theme_set(theme_minimal())

# DATA PREP  ------------------------------------------------------------

data <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/master/data/infection_status_csv.txt",
                 sep=",")
# use if offline
#library(readr)
#data <- read_csv("LSHTM/Project/Data/infection_status.csv")

# Group age into intervals 
# bi-monthly for 0-2 years and 6-monthly for 2-5 years
data$agegrp <- cut(data$age_days,
                   breaks = c(seq(0,730, by = 30.25*2),
                              seq(909,2000, by = 30.25*6)), 
                   include.lowest = T, right = F)
# Divide by season of birth
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)

data <- data %>% 
  mutate(Birth_mo = birthday %>% month(),
         season_birth = case_when (Birth_mo %in% spring ~ "Spring",
                                   Birth_mo %in% summer ~ "Summer",
                                   Birth_mo %in% autumn ~ "Autumn",
                                   Birth_mo %in% winter ~ "Winter"),
         visitnursery_child = case_when (visitnursery_child == 0 ~ FALSE,
                                         visitnursery_child == 1 ~ TRUE))

#different groupings to plot different things
# only grouped by age
data_no_season <- data %>% group_by(agegrp) %>% 
  dplyr::summarise(agemid = round(median(age_days)), # Age midpoint in age group
                   N = n(), # Total N in age group
                   nconv = sum(infection))

# grouped by age and season of birth
data_season <- data %>% group_by(agegrp, season_birth) %>% 
  dplyr::summarise(agemid = round(median(age_days)), # Age midpoint in age group
                   N = n(), # Total N in age group
                   nconv = sum(infection)) # n seroconverted in age group

# Calculate seroprevalence and binomial confidence intervals
data_no_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_season$nconv, data_no_season$N, method="exact")[,c("mean","lower","upper")]
data_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_season$nconv, data_season$N, method="exact")[,c("mean","lower","upper")]

# Plot the whole data frame (difficult to read)
data_season %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = season_birth)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, 
                    ymax = seroprev_up95, colour = season_birth)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(colour = "season of birth")

# Plot data grouped by age only
data_no_season %>% ggplot() + 
  geom_point(aes(x = agemid, y = seroprev_mean)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (years)") +
  scale_x_continuous(breaks = c(0,365,730,1095,1460,1825), labels = c(0, 1, 2, 3, 4 ,5))

# Plot by season
spring_df <- data_season %>% subset(season_birth == 'Spring')
spring_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in spring")

summer_df <- data_season %>% subset(season_birth == 'Summer')
summer_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in summer")

autumn_df <- data_season %>% subset(season_birth == 'Autumn')
autumn_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in autumn")

winter_df <- data_season %>% subset(season_birth == 'Winter')
winter_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in winter")

# MODEL EQUATION  ------------------------------------------------------------
# Define function for waning following the Erlang distribution
# We will use it to model the waning of maternal immunity
erlang.decay = function(immunity_init, T, k=2, t) {
  res = 0
  for(n in 0:(k-1)){
    res = res + 1/factorial(n)*exp(-T*t)*(T*t)^n } # T is the rate, t is the age
  return(immunity_init * res)
}

# Notes:
# - The states are integrated on a log-scale to avoid negative states, which is 
# why we log-transform them when in the model input and exponentiate them inside the model

model <- function(theta, age, inits, data) {
  
  catalytic <- function(age, state, param, data) {
    
    # FOI / seroconversion rate
    # Define booleans to know which season the cohort is in
    spring_FOI_sp = 0 # FOI in spring for those born in spring
    summer_FOI_sp = 0 # FOI in summer for those born in spring
    autumn_FOI_sp = 0
    winter_FOI_sp = 0
    spring_FOI_sm = 0
    summer_FOI_sm = 0
    autumn_FOI_sm = 0
    winter_FOI_sm = 0
    spring_FOI_au = 0 # FOI in spring for those born in autumn
    summer_FOI_au = 0
    autumn_FOI_au = 0
    winter_FOI_au = 0
    spring_FOI_wt = 0
    summer_FOI_wt = 0
    autumn_FOI_wt = 0
    winter_FOI_wt = 0
    
    # Baseline contact parameter for each birth cohort
    contacts_sp = 1
    contacts_sm = 1
    contacts_au = 1
    contacts_wt = 1
    
    
    # Get seasonal FOI
    if ( (age <= 30.41*3)                                      # FOI of the season in which the children were born
         || ((age >= 365) & (age <= (365+30.41*3))) 
         || ((age >= 2*365) & (age <= (2*365+30.41*3))) 
         || ((age >= 3*365) & (age <= (3*365+30.41*3))) 
         || ((age >= 4*365) & (age <= (4*365+30.41*3))) 
         || ((age >= 5*365) & (age <= (5*365+30.41*3)))){
      spring_FOI_sp = 1
      summer_FOI_sm = 1
      autumn_FOI_au = 1
      winter_FOI_wt = 1
    } 
    
    if ( (age > 30.41*3 & age <= 30.41*6 )                       # FOI of the season after the children were born
         || ((age > 365 + 30.41*3) & (age <= (365+30.41*6))) 
         || ((age > 2*365 + 30.41*3) & (age <= (2*365+30.41*6))) 
         || ((age > 3*365 + 30.41*3) & (age <= (3*365+30.41*6))) 
         || ((age > 4*365 + 30.41*3) & (age <= (4*365+30.41*6))) 
         || ((age > 5*365 + 30.41*3) & (age <= (5*365+30.41*6))) ){
      summer_FOI_sp = 1
      autumn_FOI_sm = 1
      winter_FOI_au = 1
      spring_FOI_wt = 1
    } 
    
    if ( (age > 30.41*6 & age <= 30.41*9 )                     # FOI 2 seasons after birth
         || ((age > 365+30.41*6) & (age <= (365+30.41*9))) 
         || ((age > 2*365 + 30.41*6) & (age <= (2*365+30.41*9))) 
         || ((age > 3*365 + 30.41*6) & (age <= (3*365+30.41*9))) 
         || ((age > 4*365 + 30.41*6) & (age <= (4*365+30.41*9))) 
         || ((age > 5*365 + 30.41*6) & (age <= (5*365+30.41*9))) ){
      autumn_FOI_sp = 1
      winter_FOI_sm = 1
      spring_FOI_au = 1
      summer_FOI_wt = 1
    } 
    
    if ( (age > 30.41*9 & age <= 30.41*12 )                     # FOI 3 seasons after birth
         || ((age > 365+30.41*9) & (age <= (365+30.41*12))) 
         || ((age > 2*365 + 30.41*9) & (age <= (2*365+30.41*12))) 
         || ((age > 3*365 + 30.41*9) & (age <= (3*365+30.41*12))) 
         || ((age > 4*365 + 30.41*9) & (age <= (4*365+30.41*12))) 
         || ((age > 5*365 + 30.41*9) & (age <= (5*365+30.41*12))) ){
      winter_FOI_sp = 1
      spring_FOI_sm = 1
      summer_FOI_au = 1
      autumn_FOI_wt = 1
    } 
    
    # Update contact and daycare parameters 
    '#the median number of contacts for children aged 9-12 months is x/y*(median 
    number of contacts of children aged 0-6 months) where:
    - x is the average number of contacts at that age
    - y is the average number of contacts at ages 9-12?#'
    if ((age > 30.41*6) & (age <= 30.41*12)){ 
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 3/4
      contacts_wt = 4/4
    } 
    if ((age > 30.41*12) & (age <= 30.41*18)){
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 5/4
      contacts_wt = 4/4
    }
    if ((age > 30.41*18) & (age <= 30.41*24)){
      contacts_sp = 4.5/4
      contacts_sm = 7/1
      contacts_au = 4/4
      contacts_wt = 4.5/4
    }
    if ((age > 30.41*24) & (age <= 30.41*30)){
      contacts_sp = 6/4
      contacts_sm = 7.5/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*30) & (age <= 30.41*36)){
      contacts_sp = 8/4
      contacts_sm = 6/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*36) & (age <= 30.41*42)){
      contacts_sp = 20/4
      contacts_sm = 9/1
      contacts_au = 9/4
      contacts_wt = 7.5/4
    }
    if ((age > 30.41*42) & (age <= 30.41*48)){
      contacts_sp = 14/4
      contacts_sm = 11/1
      contacts_au = 12/4
      contacts_wt = 11/4
    }
    if ((age > 30.41*48) & (age <= 30.41*54)){
      contacts_sp = 29/4
      contacts_sm = 25/1
      contacts_au = 12/4  # NB: no children born in autumn are aged 48-54 months
      contacts_wt = 21/4
    }
    if (age > 30.41*54){
      contacts_sp = 25/4
      contacts_sm = 13/1
      contacts_au = 11.5/4
      contacts_wt = 26.5/4
    }
    
    
    # FOI for each birth cohort
    # participants not attending daycare
    # FOI for each birth cohort
    # The FOI is a sum of the seasonal FOIs and the FOI due to contacts
    # the seasonal FOIs are a scaled version of the summer FOI (?)
    lambda_sp = (param[["M"]] + param[["P"]]) * spring_FOI_sp + 
      param[["M"]] * summer_FOI_sp + 
      (param[["M"]] + param[["A"]]) * autumn_FOI_sp +
      (param[["M"]]+param[["W"]]) * winter_FOI_sp +
      param[["C"]] * contacts_sp
    lambda_sm = (param[["M"]] + param[["P"]]) * spring_FOI_sm + 
      param[["M"]] * summer_FOI_sm + 
      (param[["M"]] + param[["A"]]) * autumn_FOI_sm + 
      (param[["M"]] + param[["W"]]) * winter_FOI_sm +
      param[["C"]] * contacts_sm
    lambda_au = (param[["M"]] + param[["P"]]) * spring_FOI_au + 
      param[["M"]] * summer_FOI_au + 
      (param[["M"]] + param[["A"]]) * autumn_FOI_au + 
      (param[["M"]] + param[["W"]]) * winter_FOI_au +
      param[["C"]] * contacts_au
    lambda_wt = (param[["M"]]+param[["P"]]) * spring_FOI_wt + 
      param[["M"]] * summer_FOI_wt +
      (param[["M"]] + param[["A"]]) * autumn_FOI_wt + 
      (param[["M"]] + param[["W"]]) * winter_FOI_wt +
      param[["C"]] * contacts_wt
    
    # waning maternal immunity, follows an Erlang distribution
    mu = erlang.decay(immunity_init = 1, T =  param[["V"]], k = 2, t = age) # test with rate set to 1/6 months. We probbaly want to estimate it later
    
    # states 
    # proportion with maternal immunity
    M_sp = exp(state[1]) # born in spring
    M_sm = exp(state[2]) # born in summer
    M_au = exp(state[3]) # born in autumn
    M_wt = exp(state[4]) # born in winter
    # Susceptible
    S_sp = exp(state[5]) # susceptible born in spring
    S_sm = exp(state[6]) # susceptible born in summer
    S_au = exp(state[7]) # susceptible born in autumn
    S_wt = exp(state[8]) # susceptible born in winter
    #Seroconverted
    Z_sp = exp(state[9]) # seroconverted after infection born in spring
    Z_sm = exp(state[10]) # seroconverted after infection born in summer
    Z_au = exp(state[11]) # seroconverted after infection born in autumn
    Z_wt = exp(state[12]) # seroconverted after infection born in winter
    
    # changes in states
    dM_sp = - mu*M_sp
    dM_sm = - mu*M_sm
    dM_au = - mu*M_au
    dM_wt = - mu*M_wt
    
    dS_sp = + mu*M_sp - lambda_sp*S_sp
    dS_sm = + mu*M_sm - lambda_sm*S_sm
    dS_au = + mu*M_au - lambda_au*S_au
    dS_wt = + mu*M_wt - lambda_wt*S_wt
    
    dZ_sp = + lambda_sp*S_sp
    dZ_sm = + lambda_sm*S_sm
    dZ_au = + lambda_au*S_au
    dZ_wt = + lambda_wt*S_wt
    
    
    return(list(c(dM_sp/M_sp,dM_sm/M_sm, dM_au/M_au,dM_wt/M_wt,
                  dS_sp/S_sp, dS_sm/S_sm, dS_au/S_au,dS_wt/S_wt,
                  dZ_sp/Z_sp, dZ_sm/Z_sm, dZ_au/Z_au, dZ_wt/Z_wt),
                lambda_sp = lambda_sp, lambda_sm = lambda_sm, 
                lambda_au = lambda_au, lambda_wt = lambda_wt,
                mu = mu))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sp = log(inits[["M_sp"]]),
                             M_sm = log(inits[["M_sm"]]),
                             M_au = log(inits[["M_au"]]),
                             M_wt = log(inits[["M_wt"]]),
                             S_sp = log(inits[["S_sp"]]),
                             S_sm = log(inits[["S_sm"]]),
                             S_au = log(inits[["S_au"]]),
                             S_wt = log(inits[["S_wt"]]),
                             Z_sp = log(inits[["Z_sp"]]),
                             Z_sm = log(inits[["Z_sm"]]),
                             Z_au = log(inits[["Z_au"]]),
                             Z_wt = log(inits[["Z_wt"]])),
                         times = age, 
                         func = catalytic, 
                         parms = theta, 
                         data = data,
                         method ="lsoda",
                         verbose = F))
  
  # Cumulative seroconversion
  traj$conv_spring <- exp(traj$Z_sp) # cumulative seroconversion in spring cohort (=observed state)
  traj$conv_summer <- exp(traj$Z_sm) # cumulative seroconversion in summer cohort (=observed state)
  traj$conv_autumn <- exp(traj$Z_au) # cumulative seroconversion in autumn cohort (=observed state)
  traj$conv_winter <- exp(traj$Z_wt) # cumulative seroconversion in winter cohort (=observed state)
  
  traj$Z_all <- 0.26*traj$Z_sp + # total seroconverted is the sum by birth cohort
    0.29*traj$Z_sm + # we scale by the proportion of children born in each season
    0.24*traj$Z_au + 
    0.20*traj$Z_wt
  traj$conv <- exp(traj$Z_all)
  
  # Incident seroconversion
  traj$inc_spring <- c(inits[["Z_sp"]], diff(exp(traj$Z_sp))) # incident seroconversion in spring cohort
  traj$inc_summer <- c(inits[["Z_sm"]], diff(exp(traj$Z_sm))) # incident seroconversion in summer cohort
  traj$inc_autumn <- c(inits[["Z_au"]], diff(exp(traj$Z_au))) # incident seroconversion in autumn cohort
  traj$inc_winter <- c(inits[["Z_wt"]], diff(exp(traj$Z_wt))) # incident seroconversion in winter cohort
  
  inc_all <- 0.26*diff(exp(traj$Z_sp)) + # total incident seroconversion is the sum by birth cohort
    0.29*diff(exp(traj$Z_sm)) + # we scale by the proportion of children born in each season
    0.24*diff(exp(traj$Z_au)) + 
    0.20*diff(exp(traj$Z_wt))
  inc_all_inits <- 0.26*inits[["Z_sp"]] + # initial incident seroconversion
    0.29*inits[["Z_sm"]] + 
    0.24*inits[["Z_au"]] +
    0.20*inits[["Z_wt"]]
  traj$inc_all <- c(inc_all_inits, inc_all) # total incident seroconversion at all time points
  
  
  return(traj)
  
}


# TRAJECTORY SIMULATION  ------------------------------------------------------------

# This function draws n samples from the posterior to calculate the posterior predictive uncertainty (95%)

maketrajsim <- function(trace, theta, age, model, inits, ndraw, data) {
  
  #Draw n fitted parameter vectors theta from the MCMC object
  sample <- getSample(trace, parametersOnly = TRUE, thin = 1, numSamples = ndraw) #trace is a sampler, parametersOnly = T means that likelihood, posterior and prior values are not provided in the output, thin = thinning parameter
  
  traj.rep <- adply(.data = sample, .margins = 1, .progress = "text", .parallel = F, .fun = function(x) { #split sample by 1 = rows and apply function
    
    #Define the theta as the estimated parameters from the current draw and the fixed parameters
    theta.sample <- c(x, theta[!is.na(theta)])
    
    #Simulate trajectory for selected theta
    traj <- match.fun(model)(theta.sample, age, inits, data) #match.fun extracts the underlying function
    traj <- cbind(as.data.frame(t(x)), traj)
  })
  
  colnames(traj.rep)[1] <- "replicate"
  return(traj.rep)
  
}


# THETA ---------------------------------------------------------

# P = mean FOI (proportion infected per day) in spring
# M = mean FOI for children in summer
# A = mean FOI for children in autumn
# W = mean FOI for children in winter
# V = waning rate of protection offered by maternal immunity from the vaccine
# C = contact parameter
theta <- c(P = 0.00001, M = 0.02002, A = 0.00003, W = 0.00004, C = 0.02, V = 1/180) # these are just random values, to be fitted

# INITS ---------------------------------------------------------
inits <- c(M_sp = (1-2*1e-12), M_sm = (1-2*1e-12), M_au= (1-2*1e-12), M_wt = (1-2*1e-12),
           S_sp = 1e-12, S_sm = 1e-12, S_au = 1e-12, S_wt = 1e-12,
           Z_sp = 1e-12, Z_sm = 1e-12, Z_au = 1e-12, Z_wt = 1e-12)
# initial conditions for the states (as proportions)
# --> since we integrate on a log-scale, the initial conditions cannot be 0 (not defined on a log-scale)

# SIMULATION TIME  ---------------------------------------------------------
data_season <- arrange(data_season, agemid)
agepred <- data_season$agemid

# TEST MODEL  --------------------------------------------------------
# Notes: lambda should vary over time, mu shouldn't
test <- model(theta, agepred, inits, data_season)
ggplot(test) + geom_line(aes(x=time, y = conv)) # cumulative conversion

# Incident seroconversion
ggplot(test) + geom_line(aes(x=time, y = inc_spring), colour = "red") +
  geom_line(aes(x=time, y = inc_summer), colour = "yellow") +
  geom_line(aes(x=time, y = inc_autumn), colour = "green") +
  geom_line(aes(x=time, y = inc_winter), colour = "blue")
ggplot(test) + geom_line(aes(x=time, y = inc_all))

ggplot(test) + geom_line(aes(x=time, y=lambda_sp)) + 
  xlab("Age (days)") + ylab("FOI") + labs(color = 'Daycare') + ggtitle("FOI for the spring birth cohort") 

ggplot(test) + geom_line(aes(x=time, y=lambda_sm))+ 
  xlab("Age (days)") + ylab("FOI") + ggtitle("FOI for the summer birth cohort") +
  labs(color = 'Daycare')

ggplot(test) + geom_line(aes(x=time, y=lambda_au)) +
  xlab("Age (days)") + ylab("FOI") + ggtitle("FOI for the autumn birth cohort ") +
  labs (color = 'Daycare')

ggplot(test) + geom_line(aes(x=time, y=lambda_wt)) +
  xlab("Age (days)") + ylab("FOI") + ggtitle("FOI for the winter birth cohort") + 
  labs(color = 'Daycare')

ggplot(test) + geom_line(aes(x=time, y=mu)) + xlab("Age (days)") + ylab("Waning rate of maternal immunity")

# LOG LIKELIHOOD FUNCTION ---------------------------------------------------------

loglik <- function(theta, age, data, model, inits) {
  
  traj <- match.fun(model)(theta, age, inits)
  
  # Whole data
  nconv <- data$nconv[!is.na(data$nconv)] # n seroconverted at each age point  (data)
  N <- data$N[!is.na(data$nconv)] # total N at each age point  (data) 
  prob <- traj$conv[traj$time %in% data$agemid] # proportion seroconverted at each age point (model output)
  
  ll_all <- sum(dbinom(x = nconv,
                       size = N,
                       prob = prob,
                       log = TRUE), na.rm = TRUE)
  
  return(ll_all)
  
} 

# Test function
loglik(theta, agepred, data_season, model, inits)

# Wrapper for BT: loglik can only take the fitted parameters as argument
loglik_wrapper <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(loglik(theta = parX,
                age = agepred, 
                data = data_season,
                model = match.fun(model),
                inits = inits))
} 

# FITTING -------------------------------------------

# Estimated params
estpars <- c("P", "M", "A", "W", "C", "V") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params


# Priors
lower = c(P = 0, M = 0, A = 0, W = 0, C = 0, V = 0.001)
upper = c(P = 0.1, M = 0.1, A = 0.1, W = 0.1, C = 0.1, V= 0.99)

prior <- createUniformPrior(lower=lower[estpars], 
                            upper=upper[estpars])
# MCMC settings
nchains <- 2
cpus <- 1 # or 2 if you want parallel, but it does not seem to be faster?
mcmc_settings <- list(iterations = 2*80000, 
                      nrChains = nchains)
sampler <- "DEzs"

if (cpus == 1) {
  bayesianSetup <- createBayesianSetup(prior = prior,
                                       likelihood = loglik_wrapper,
                                       names = names(theta[index]),
                                       parallel = FALSE)
  
  system.time({trace <- runMCMC(bayesianSetup = bayesianSetup, 
                                sampler = sampler, 
                                settings = mcmc_settings)})
  
}

# DIAGNOSTICS -----------------------------------------------

plot(trace) 

# burn-in
nburn <- 10000
plot(trace, parametersOnly = TRUE, start = nburn)

# check convergence and correlations
gelmanDiagnostics(trace, plot = TRUE, start = nburn)
correlationPlot(getSample(trace, parametersOnly = TRUE, coda = TRUE, start = nburn), density = "smooth", thin = 50)
marginalPlot(trace, prior = T, singlePanel = T, start = nburn, nDrawsPrior = 1000)

# remove burn-in for trajsim simulation
tracefinal <- getSample(trace, parametersOnly = TRUE, coda = TRUE, start = nburn)
plot(tracefinal)
effectiveSize(tracefinal)

# Posterior summary
summary(tracefinal)

# save the trace
saveRDS(trace, "strat_model_trace_seasonal_FOI_contacts_mat_vaccine_no_daycare_Erlang.rds")


# POSTPROCESSING AND RESULTS -----------------------------------

# Calculate simulated trajectory quantiles
trajsim <- maketrajsim(tracefinal, theta, agepred, model, inits, 1000)

# Proportion converted overall
trajquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                             function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles) <- c("agemid", "low95", "median", "up95")

# Overall incidence overtime
trajquantiles_inc <- plyr::ddply(.data = trajsim, .variables = "time", 
                             function(x) quantile(x[,"inc_all"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_inc) <- c("agemid", "low95", "median", "up95")


# FOI per birth cohort
lambda_spquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"lambda_sp"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_spquantiles) <- c("agemid", "low95", "median", "up95")

lambda_smquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"lambda_sm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_smquantiles) <- c("agemid", "low95", "median", "up95")

lambda_auquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"lambda_au"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_auquantiles) <- c("agemid", "low95", "median", "up95")

lambda_wtquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"lambda_wt"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_wtquantiles) <- c("agemid", "low95", "median", "up95")

# Parameter estimates
Pquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"P"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Pquantiles) <- c("agemid", "low95", "median", "up95")

Mquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"M"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Mquantiles) <- c("agemid", "low95", "median", "up95")

Aquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"A"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Aquantiles) <- c("agemid", "low95", "median", "up95")

Wquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"W"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Wquantiles) <- c("agemid", "low95", "median", "up95")

Cquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"C"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Cquantiles) <- c("agemid", "low95", "median", "up95")

Vquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"V"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(vquantiles) <- c("agemid", "low95", "median", "up95")

# Waning rate
mu_quantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"mu"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(mu_quantiles) <- c("agemid", "low95", "median", "up95")


# Cumulative proportion seroconverted by season of birth
spring_conv_quantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"conv_spring"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(spring_conv_quantiles) <- c("agemid", "low95", "median", "up95")

summer_conv_quantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"conv_summer"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(summer_conv_quantiles) <- c("agemid", "low95", "median", "up95")

autumn_conv_quantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"conv_autumn"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(autumn_conv_quantiles) <- c("agemid", "low95", "median", "up95")

winter_conv_quantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"conv_winter"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(winter_conv_quantiles) <- c("agemid", "low95", "median", "up95")

# Incident proportion seroconverted by season of birth
spring_conv_quantiles_inc <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"inc_spring"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(spring_conv_quantiles_inc) <- c("agemid", "low95", "median", "up95")

summer_conv_quantiles_inc <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"inc_summer"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(summer_conv_quantiles_inc) <- c("agemid", "low95", "median", "up95")

autumn_conv_quantiles_inc <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"inc_autumn"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(autumn_conv_quantiles_inc) <- c("agemid", "low95", "median", "up95")

winter_conv_quantiles_inc <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"inc_winter"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(winter_conv_quantiles_inc) <- c("agemid", "low95", "median", "up95")

# Plot fit and FOI
# Whole dataset
fit <- ggplot() + theme_bw() + ggtitle("Model fit") +
  geom_point(data = data_no_season, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = data_no_season, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = trajquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = trajquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit

# Incidence
ggplot() + theme_bw() + ggtitle("Modelled incidence") +
  geom_point(data = data_no_season, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = data_no_season, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = trajquantiles_inc, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = trajquantiles_inc, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("Proportion seroconverted") 

# Stratified by season of birth
fit_season <- ggplot() + theme_bw() + ggtitle("Model fit on data stratified by season") +
  geom_point(data = data_season, aes(x = agemid, y = seroprev_mean, color = season_birth)) +
  geom_linerange(data = data_season, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, color = season_birth)) +
  geom_ribbon(data = trajquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = trajquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") + 
  labs(color = "Season of birth")

fit_season

# Spring birth cohort
fit_sp <- ggplot() + theme_bw() + ggtitle("Model fit on spring birth cohort") +
  geom_point(data = spring_df, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = spring_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = spring_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = spring_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted")

fit_sp

# Summmer birth cohort 
fit_sm <- ggplot() + theme_bw() + ggtitle("Model fit on summer birth cohort") +
  geom_point(data = summer_df, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = summer_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = summer_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = summer_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted")

fit_sm

# Autumn birth cohort
fit_au <- ggplot() + theme_bw() + ggtitle("Model fit on autumn birth cohort") +
  geom_point(data = autumn_df, aes(x = agemid, y = seroprev_mean,)) +
  geom_linerange(data = autumn_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = autumn_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = autumn_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted")

fit_au

# Winter birth cohort
fit_wt <- ggplot() + theme_bw() + ggtitle("Model fit on winter birth cohort") +
  geom_point(data = winter_df, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = winter_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = winter_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = winter_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted")

fit_wt

# FOI for spring birth cohort
lambda_sp <- ggplot() + theme_bw() + ggtitle("FOI for the spring birth cohort (no daycare)") +
  geom_ribbon(data = lambda_spquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = lambda_spquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("FOI") + 
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))
lambda_sp

# FOI for summer birth cohort
lambda_sm <- ggplot() + theme_bw() + ggtitle("FOI for the summer birth cohort (no daycare)") +
  geom_ribbon(data = lambda_smquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = lambda_smquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("FOI") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))
lambda_sm

# FOI for autumn birth cohort
lambda_au <- ggplot() + theme_bw() + ggtitle("FOI for the autumn birth cohort (no daycare)") +
  geom_ribbon(data = lambda_auquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = lambda_auquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("FOI") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))
lambda_au

# FOI for winter birth cohort
lambda_wt <- ggplot() + theme_bw() + ggtitle("FOI for the winter birth cohort (no daycare)") +
  geom_ribbon(data = lambda_wtquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = lambda_wtquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("FOI") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))
lambda_wt

# Waning rate of maternal immunity (should not be constant)
mu <- ggplot() + theme_bw() + ggtitle("Rate of maternal immunity waning ") +
  geom_ribbon(data = mu_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = mu_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("Maternal immunity waning rate") 
mu

# Saving the useful files ------------------------------------------------------------
path <- ("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/CSV files/Vaccination/")
write.csv(lambda_spquantiles, paste0(path, "strat_lambda_sp.csv"))
write.csv(lambda_smquantiles, paste0(path, "strat_lambda_sm.csv"))
write.csv(lambda_auquantiles, paste0(path, "strat_lambda_au.csv"))
write.csv(lambda_wtquantiles, paste0(path, "strat_lambda_wt.csv"))

write.csv(Pquantiles, paste0(path, "strat_Pquantiles.csv"))
write.csv(Mquantiles, paste0(path, "strat_Mquantiles.csv"))
write.csv(Aquantiles, paste0(path, "strat_Aquantiles.csv"))
write.csv(Wquantiles, paste0(path, "strat_Wquantiles.csv"))
write.csv(Cquantiles, paste0(path, "strat_Cquantiles.csv"))
write.csv(wquantiles, paste0(path, "strat_Immunity_quantiles.csv"))

# Cumulative seroconversion
write.csv(trajquantiles, paste0(path, "strat_trajquantiles.csv"))
write.csv(trajsim, paste0(path, "strat_trajsim.csv"))
write.csv(spring_conv_quantiles, paste0(path, "strat_spring_conv_tot.csv"))
write.csv(summer_conv_quantiles, paste0(path, "strat_summer_conv_tot.csv"))
write.csv(autumn_conv_quantiles, paste0(path, "strat_autumn_conv_tot.csv"))
write.csv(winter_conv_quantiles, paste0(path, "strat_winter_conv_tot.csv"))

# Incident seroconversion
write.csv(trajquantiles_inc, paste0(path, "Incidence/strat_trajquantiles_inc.csv"))
write.csv(trajsim, paste0(path, "Incidence/strat_trajsim.csv"))
write.csv(spring_conv_quantiles_inc, paste0(path, "Incidence/strat_spring_conv_inc.csv"))
write.csv(summer_conv_quantiles_inc, paste0(path, "Incidence/strat_summer_conv_inc.csv"))
write.csv(autumn_conv_quantiles_inc, paste0(path, "Incidence/strat_autumn_conv_inc.csv"))
write.csv(winter_conv_quantiles_inc, paste0(path, "Incidence/strat_winter_conv_inc.csv"))
