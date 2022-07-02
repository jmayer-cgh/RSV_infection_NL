################################################################
# RSV seroconversion MSc project
# Adding season of birth
# Author: Julia Mayer
# Last updated: 01.07.2022
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

# Group age into intervals 
# bi-monthly for 0-2 years and 6-monthly for 2-5 years
data$agegrp <- cut(data$age_days,
                   breaks=c(seq(0,730, by=30.25*2),
                            seq(909,2000, by=30.25*6)), 
                   include.lowest = T, right=F)
# Divide by season of birth
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)

data <- data %>%
  mutate(
    Birth_mo = birthday %>% month()
  )%>%
  mutate(
    season_birth = case_when(
      Birth_mo %in% spring ~ "Spring",
      Birth_mo %in% summer ~ "Summer",
      Birth_mo %in% autumn ~ "Autumn",
      Birth_mo %in% winter ~ "Winter")
  )

data <- data %>% dplyr::group_by(agegrp, season_birth) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection)) # n seroconverted in age group


# Calculate seroprevalence and binomial confidence intervals
data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data$nconv, data$N, method="exact")[,c("mean","lower","upper")]

# Plot the whole data frame (difficult to read)
ggplot(data) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = season_birth)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = season_birth)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs(colour = "season of birth")

#Plot by season
spring.df <- subset(data, season_birth == 'Spring')
ggplot(spring.df) +
  geom_point(aes(x=agemid, y=seroprev_mean), color = 'black') +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95),colour ='black') +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs(title ="Proportion seroconverted born in spring")

summer.df <- subset(data, season_birth == 'Summer')
ggplot(summer.df) +
  geom_point(aes(x=agemid, y=seroprev_mean), color = 'black') +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95),colour ='black') +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs(title ="Proportion seroconverted born in summer")

autumn.df <- subset(data, season_birth == 'Autumn')
ggplot(autumn.df) +
  geom_point(aes(x=agemid, y=seroprev_mean), color = 'black') +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95),colour ='black') +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs(title ="Proportion seroconverted born in autumn")

winter.df <- subset(data, season_birth == 'Winter')
ggplot(winter.df) +
  geom_point(aes(x=agemid, y=seroprev_mean), color = 'black') +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95),colour ='black') +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs(title ="Proportion seroconverted born in winter")

# MODEL EQUATION  ------------------------------------------------------------
# Notes:
# - The states are integrated on a log-scale to avoid negative states, which is why we log-transform them when in the model input and exponentiate them inside the model

model <- function(theta, age, inits) {
  
  catalytic <- function(age, state, param) {
    
    # FOI / seroconversion rate by season
    lambda_sp = param[["P"]] + age * 0 #constant FOI for now, change 0 to something else later
    lambda_sm = param[["M"]] + age * 0 #constant FOI for now, change 0 to something else later
    lambda_au = param[["A"]] + age * 0 #constant FOI for now, change 0 to something else later
    lambda_wt = param[["W"]] + age * 0 #constant FOI for now, change 0 to something else later
    mu = param[["B"]] # waning maternal immunity
    
    # states by season
    M_sp = exp(state[1]) # proportion with maternal immunity
    S_sp = exp(state[2]) # susceptible
    Z_sp = exp(state[3]) # seroconverted after infection
    M_sm = exp(state[4]) # proportion with maternal immunity
    S_sm = exp(state[5]) # susceptible
    Z_sm = exp(state[6]) # seroconverted after infection
    M_au = exp(state[7]) # proportion with maternal immunity
    S_au = exp(state[8]) # susceptible
    Z_au = exp(state[9]) # seroconverted after infection
    M_wt = exp(state[10]) # proportion with maternal immunity
    S_wt = exp(state[11]) # susceptible
    Z_wt = exp(state[12]) # seroconverted after infection
    
    # changes in states
    dM_sp = - mu*M_sp
    dS_sp = + mu*M_sp - lambda_sp*S_sp
    dZ_sp = + lambda_sp*S_sp
    dM_sm = - mu*M_sm
    dS_sm = + mu*M_sm - lambda_sm*S_sm
    dZ_sm = + lambda_sm*S_sm
    dM_au = - mu*M_au
    dS_au = + mu*M_au - lambda_au*S_au
    dZ_au = + lambda_au*S_au
    dM_wt = - mu*M_wt
    dS_wt = + mu*M_wt - lambda_wt*S_wt
    dZ_wt = + lambda_wt*S_wt
    
    return(list(c(dM_sp/M_sp,dS_sp/S_sp, dZ_sp/Z_sp,
                  dM_sm/M_sm,dS_sm/S_sm, dZ_sm/Z_sm,
                  dM_au/M_au,dS_au/S_au, dZ_au/Z_au,
                  dM_wt/M_wt,dS_wt/S_wt, dZ_wt/Z_wt),
                  lambda_sp=lambda_sp,lambda_sm=lambda_sm,
                  lambda_au=lambda_au,lambda_wt=lambda_wt,mu=mu))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sp=log(inits[["M_sp"]]),
                             S_sp=log(inits[["S_sp"]]),
                             Z_sp=log(inits[["Z_sp"]]),
                             M_sm=log(inits[["M_sm"]]),
                             S_sm=log(inits[["S_sm"]]),
                             Z_sm=log(inits[["Z_sm"]]),
                             M_au=log(inits[["M_au"]]),
                             S_au=log(inits[["S_au"]]),
                             Z_au=log(inits[["Z_au"]]),
                             M_wt=log(inits[["M_wt"]]),
                             S_wt=log(inits[["S_wt"]]),
                             Z_wt=log(inits[["Z_wt"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         method="lsoda",
                         verbose=F))
  
  traj$conv_sp <- exp(traj$Z_sp) # cumulative seroconversion (=observed state)
  traj$inc_sp <- c(inits[["Z_sp"]], diff(exp(traj$Z_sp))) # incident seroconversion
  traj$conv_sm <- exp(traj$Z_sm) # cumulative seroconversion (=observed state)
  traj$inc_sm <- c(inits[["Z_sm"]], diff(exp(traj$Z_sm))) # incident seroconversion
  traj$conv_au <- exp(traj$Z_au) # cumulative seroconversion (=observed state)
  traj$inc_au <- c(inits[["Z_au"]], diff(exp(traj$Z_au))) # incident seroconversion
  traj$conv_wt <- exp(traj$Z_wt) # cumulative seroconversion (=observed state)
  traj$inc_wt <- c(inits[["Z_wt"]], diff(exp(traj$Z_wt))) # incident seroconversion
  
  return(traj)
  
}


# TRAJECTORY SIMULATION  ------------------------------------------------------------

# This function draws n samples from the posterior to calculate the posterior predictive uncertainty (95%)

maketrajsim <- function(trace, theta, age, model, inits, ndraw) {
  
  #Draw n fitted parameter vectors theta from the MCMC object
  sample <- getSample(trace, parametersOnly = TRUE, thin=1, numSamples=ndraw) #trace is a sampler, parametersOnly = T means that likelihood, posterior and prior values are not provided in the output, thin = thinning parameter
  
  traj.rep <- adply(.data=sample, .margins=1, .progress="text", .parallel=F, .fun=function(x) { #split sample by 1 = rows and apply function
    
    #Define the theta as the estimated parameters from the current draw and the fixed parameters
    theta.sample <- c(x, theta[!is.na(theta)])
    
    #Simulate trajectory for selected theta
    traj <- match.fun(model)(theta.sample, age, inits) #match.fun extracts the underlying function
    traj <- cbind(as.data.frame(t(x)), traj)
  })
  
  colnames(traj.rep)[1] <- "replicate"
  return(traj.rep)
  
}


# THETA ---------------------------------------------------------

# P = mean FOI (proportion infected per day) for those born in spring
# M = mean FOI (proportion infected per day) for those born in summer
# A = mean FOI (proportion infected per day) for those born in autumn
# W = mean FOI (proportion infected per day) for those born in winter
# B = rate of waning maternal immunity
theta <- c(P = 0.02, M = 0.02, A=0.02, W = 0.02, B = 0.01) # these are just random values, to be fitted

# INITS ---------------------------------------------------------
inits <- c(M_sp=0.25-8e-12, S_sp=1e-12, Z_sp=1e-12,
           M_sm=0.25-8e-12, S_sm=1e-12, Z_sm=1e-12,
           M_au=0.25-8e-12, S_au=1e-12, Z_au=1e-12,
           M_wt=0.25-8e-12, S_wt=1e-12, Z_wt=1e-12) # initial conditions for the states (as proportions)
# --> since we integrate on a log-scale, the initial conditions cannot be 0 (not defined on a log-scale)

# SIMULATION TIME  ---------------------------------------------------------
agepred <- data$agemid

# TEST MODEL  --------------------------------------------------------
test <- model(theta, agepred, inits)
ggplot(test) + geom_line(aes(x=time, y=conv_sp))
ggplot(test) + geom_line(aes(x=time, y=lambda_sp)) #should be constant
ggplot(test) + geom_line(aes(x=time, y=w))