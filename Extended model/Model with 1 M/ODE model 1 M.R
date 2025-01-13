# MSc model with 1 M compartment and only seasonal components for the FOI

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

# grouped by age, season of birth and day-care attendance
data <- data %>% subset(!is.na(visitnursery_child))
data <- data %>% group_by(agegrp, season_birth, visitnursery_child) %>% 
  dplyr::summarise(agemid = round(median(age_days)), # Age midpoint in age group
                   N = n(), # Total N in age group
                   nconv = sum(infection)) # n seroconverted in age group

# Calculate seroprevalence and binomial confidence intervals
data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data$nconv, data$N, method="exact")[,c("mean","lower","upper")]

# MODEL EQUATION  ------------------------------------------------------------
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
    # FOI for each birth cohort
    # The FOI is a sum of the seasonal FOIs and the FOI due to contacts
    # the seasonal FOIs are a scaled version of the summer FOI (?)
    lambda_sp = (param[["M"]] + param[["P"]]) * spring_FOI_sp + 
      param[["M"]] * summer_FOI_sp + 
      (param[["M"]] + param[["A"]]) * autumn_FOI_sp +
      (param[["M"]]+param[["W"]]) * winter_FOI_sp
    lambda_sm = (param[["M"]] + param[["P"]]) * spring_FOI_sm + 
      param[["M"]] * summer_FOI_sm + 
      (param[["M"]] + param[["A"]]) * autumn_FOI_sm + 
      (param[["M"]] + param[["W"]]) * winter_FOI_sm
    lambda_au = (param[["M"]] + param[["P"]]) * spring_FOI_au + 
      param[["M"]] * summer_FOI_au + 
      (param[["M"]] + param[["A"]]) * autumn_FOI_au + 
      (param[["M"]] + param[["W"]]) * winter_FOI_au
    lambda_wt = (param[["M"]]+param[["P"]]) * spring_FOI_wt + 
      param[["M"]] * summer_FOI_wt +
      (param[["M"]] + param[["A"]]) * autumn_FOI_wt + 
      (param[["M"]] + param[["W"]]) * winter_FOI_wt 
    
    # waning maternal immunity, follows an Erlang distribution
    mu = 1/182.5 # take half-life of 6 months
    
    # states
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
    dM_sp = - mu * M_sp
    dM_sm = - mu * M_sm
    dM_au = - mu * M_au
    dM_wt = - mu * M_wt
    
    dS_sp = + mu * M_sp - lambda_sp*S_sp
    dS_sm = + mu * M_sm - lambda_sm*S_sm
    dS_au = + mu * M_au - lambda_au*S_au
    dS_wt = + mu * M_wt - lambda_wt*S_wt
    
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
theta <- c(P = 0.00001, M = 0.02002, A = 0.00003, W = 0.00004) # these are just random values, to be fitted

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
test <- model(theta, agepred, inits, data_season)
ggplot(test) + geom_line(aes(x=time, y = conv)) # cumulative conversion

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
log <- loglik(theta, agepred, data_no_season, model, inits) # -859.61

theta1 <- c(P = 0.003478, M = 0.00060, A = 0.00264, W = 0.00692)
log1 <- loglik(theta1, agepred, data, model, inits) # -858.52

theta2 <- c(P = 0.003478, M = 0.060, A = 0.00264, W = 0.00692)
log2 <- loglik(theta2, agepred, data, model, inits) # -1212.24

theta3 <- c(P = 0.003478, M = 0.00060, A = 0.1, W = 0.00692)
log3 <- loglik(theta3, agepred, data, model, inits) # -1182.29

theta4 <- c(P = 0.1, M = 0, A = 0.1, W = 0.1)
log4 <- loglik(theta4, agepred, data, model, inits) # -1405.04

theta5 <- c(P = 0, M = 0.1, A = 0, W = 0)
log5 <- loglik(theta5, agepred, data, model, inits) # -1208.19