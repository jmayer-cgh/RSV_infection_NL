################################################################
# RSV seroconversion MSc project
# Adding daycare attendance as a multiplication
# Author: Julia Mayer
# Last updated: 10.08.2022
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
  ) %>%
  mutate(
    visitnursery_child = case_when(
      visitnursery_child == 0 ~ FALSE,
      visitnursery_child == 1 ~ TRUE
    )
  )

data_no_season <- data %>% dplyr::group_by(agegrp) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection))

data_no_daycare <- data %>% dplyr::group_by(agegrp, season_birth) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection)) # n seroconverted in age group

data <- subset(data, !is.na(visitnursery_child))
data <- data %>% dplyr::group_by(agegrp, season_birth, visitnursery_child) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection)) # n seroconverted in age group

# Calculate seroprevalence and binomial confidence intervals
data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data$nconv, data$N, method="exact")[,c("mean","lower","upper")]
data_no_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_season$nconv, data_no_season$N, method="exact")[,c("mean","lower","upper")]
data_no_daycare[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_daycare$nconv, data_no_daycare$N, method="exact")[,c("mean","lower","upper")]

# Plot the whole data frame (difficult to read)
ggplot(data) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = season_birth, 
                 shape = factor(visitnursery_child))) +
                 geom_errorbar(aes(x=agemid, ymin=seroprev_low95, 
                                   ymax=seroprev_up95, colour = season_birth,
                                   linetype = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(colour = "season of birth", shape = "Day-care attendance", linetype = "Day-care attendance")

ggplot(data_no_season) +
  geom_point(aes(x=agemid, y=seroprev_mean)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") 

ggplot(data_no_daycare) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = season_birth)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = season_birth)) +
  ylab("Proportion seroconverted") + xlab("age (days)")  + labs (colour = "Season of birth")

#Plot by season
spring.df <- subset(data, season_birth == 'Spring')
spring_no_daycare <- subset(data_no_daycare, season_birth == 'Spring')
ggplot(spring.df) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in spring", colour = "Day-care")

summer.df <- subset(data, season_birth == 'Summer')
summer_no_daycare <- subset(data_no_daycare, season_birth == 'Summer')
ggplot(summer.df) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in summer", colour = "Day-care attendance")

autumn.df <- subset(data, season_birth == 'Autumn')
autumn_no_daycare <- subset(data_no_daycare, season_birth == 'Autumn')
ggplot(autumn.df) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in autumn", colour = "Day-care attendance")

winter.df <- subset(data, season_birth == 'Winter')
winter_no_daycare <- subset(data_no_daycare, season_birth == 'Winter')
ggplot(winter.df) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in winter", colour = "Day-care attendance")

# MODEL EQUATION  ------------------------------------------------------------
# Notes:
# - The states are integrated on a log-scale to avoid negative states, which is why we log-transform them when in the model input and exponentiate them inside the model

model <- function(theta, age, inits, data) {
  
  catalytic <- function(age, state, param, data) {
    
    # FOI / seroconversion rate
    # Define booleans to know which season the cohort is in
    spring_FOI_sp = 0
    summer_FOI_sp = 0
    autumn_FOI_sp = 0
    winter_FOI_sp = 0
    spring_FOI_sm = 0
    summer_FOI_sm = 0
    autumn_FOI_sm = 0
    winter_FOI_sm = 0
    spring_FOI_au = 0
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
    
    # Baseline percentage of daycare attendance 
    daycare = 0
    if (age >= 30.41*9){ # assume that all children who go to day-care start at 9 months
      daycare = 0.4075758
    }
    
    # Get seasonal FOI
    if ( (age <= 30.41*3)                                      # FOI of the season in which the children were born
         || ((age >= 365) & (age <= (365+30.41*3))) 
         || ((age >= 2*365) & (age <= (2*365+30.41*3))) 
         || ((age >= 3*365) & (age <= (3*365+30.41*3))) 
         || ((age >= 4*365) & (age <= (4*365+30.41*3))) 
         || ((age >= 5*365) & (age <= (5*365+30.41*3))) ){
      spring_FOI_sp = 1
      summer_FOI_sm = 1
      autumn_FOI_au = 1
      winter_FOI_wt = 1
    } 
    
    if ( (age>30.41*3 & age <=30.41*6 )                       # FOI of the season after the children were born
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
    number of contacts of children aged 0-6 months)#'
    if ((age > 30.41*6) & (age <= 30.41*12)){ 
      contacts_sp = 5/4
      contacts_sm = 3.5/1
      contacts_au = 3/4
      contacts_wt = 4/4
    } 
    if ((age > 30.41*12) & (age <= 30.41*24)){
      contacts_sp = 5/4
      contacts_sm = 6/1
      contacts_au = 4.5/4
      contacts_wt = 4/4
    }
    if ((age > 30.41*24) & (age <= 30.41*36)){
      contacts_sp = 7/4
      contacts_sm = 6/1
      contacts_au = 6/4
      contacts_wt = 6/4
    }
    if ((age > 30.41*36) & (age <= 30.41*48)){
      contacts_sp = 18/4
      contacts_sm = 10/1
      contacts_au = 10/4
      contacts_wt = 9.5/4
    }
    if ((age > 30.41*48)){
      contacts_sp = 28/4
      contacts_sm = 18.5/1
      contacts_au = 9/4
      contacts_wt = 27/4
    }
    
    
    # FOI for each birth cohort
    # participants not attending daycare
    # FOI for each birth cohort
    lambda_sp = (param[["M"]]+param[["P"]])*spring_FOI_sp + param[["M"]]*summer_FOI_sp + 
                (param[["M"]]+param[["A"]])*autumn_FOI_sp + (param[["M"]]+param[["W"]])*winter_FOI_sp 
                 #param[["C"]]*contacts_sp
    lambda_sm = (param[["M"]]+param[["P"]])*spring_FOI_sm + param[["M"]]*summer_FOI_sm + 
                (param[["M"]]+param[["A"]])*autumn_FOI_sm + (param[["M"]]+param[["W"]])*winter_FOI_sm 
                #param[["C"]]*contacts_sm
    lambda_au = (param[["M"]]+param[["P"]])*spring_FOI_au + param[["M"]]*summer_FOI_au + 
                (param[["M"]]+param[["A"]])*autumn_FOI_au + (param[["M"]]+param[["W"]])*winter_FOI_au 
                #param[["C"]]*contacts_au
    lambda_wt = (param[["M"]]+param[["P"]])*spring_FOI_wt + param[["M"]]*summer_FOI_wt +
                 (param[["M"]]+param[["A"]])*autumn_FOI_wt + (param[["M"]]+param[["W"]])*winter_FOI_wt 
                #param[["C"]]*contacts_wt
    
    # waning maternal immunity, same for all children
    mu = param[["B"]] 
    
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
    dM_sp = -mu*M_sp
    dM_sm = -mu*M_sm
    dM_au = -mu*M_au
    dM_wt = -mu*M_wt
    dS_sp = + mu*M_sp - (1-daycare)*lambda_sp*S_sp - daycare*param[["D"]]*lambda_sp*S_sp
    dS_sm = + mu*M_sm - (1-daycare)*lambda_sm*S_sm - daycare*param[["D"]]*lambda_sm*S_sm
    dS_au = + mu*M_au - (1-daycare)*lambda_au*S_au - daycare*param[["D"]]*lambda_au*S_au 
    dS_wt = + mu*M_wt - (1-daycare)*lambda_wt*S_wt - daycare*param[["D"]]*lambda_wt*S_wt
    dZ_sp = + (1-daycare)*lambda_sp*S_sp + daycare*param[["D"]]*lambda_sp*S_sp
    dZ_sm = + (1-daycare)*lambda_sm*S_sm + daycare*param[["D"]]*lambda_sm*S_sm
    dZ_au = + (1-daycare)*lambda_au*S_au + daycare*param[["D"]]*lambda_au*S_au
    dZ_wt = + (1-daycare)*lambda_wt*S_wt + daycare*param[["D"]]*lambda_wt*S_wt
    
    
    return(list(c(dM_sp/M_sp,dM_sm/M_sm, dM_au/M_au,dM_wt/M_wt,
                  dS_sp/S_sp, dS_sm/S_sm, dS_au/S_au,dS_wt/S_wt,
                  dZ_sp/Z_sp, dZ_sm/Z_sm, dZ_au/Z_au, dZ_wt/Z_wt),
                lambda_sp=lambda_sp, lambda_sm = lambda_sm, 
                lambda_au = lambda_au, lambda_wt = lambda_wt,
                mu=mu))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sp=log(inits[["M_sp"]]),
                             M_sm=log(inits[["M_sm"]]),
                             M_au=log(inits[["M_au"]]),
                             M_wt=log(inits[["M_wt"]]),
                             S_sp=log(inits[["S_sp"]]),
                             S_sm=log(inits[["S_sm"]]),
                             S_au=log(inits[["S_au"]]),
                             S_wt=log(inits[["S_wt"]]),
                             Z_sp=log(inits[["Z_sp"]]),
                             Z_sm=log(inits[["Z_sm"]]),
                             Z_au=log(inits[["Z_au"]]),
                             Z_wt=log(inits[["Z_wt"]])),
                             #Z_all = log(inits[["Z_all"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         data = data,
                         method="lsoda",
                         verbose=F))
  
  traj$conv_spring <- exp(traj$Z_sp) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring <- c(inits[["Z_sp"]], diff(exp(traj$Z_sp))) # incident seroconversion
  traj$conv_summer <- exp(traj$Z_sm) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer <- c(inits[["Z_sm"]], diff(exp(traj$Z_sm))) # incident seroconversion
  traj$conv_autumn <- exp(traj$Z_au) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn <- c(inits[["Z_au"]], diff(exp(traj$Z_au))) # incident seroconversion
  traj$conv_winter <- exp(traj$Z_wt) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter <- c(inits[["Z_wt"]], diff(exp(traj$Z_wt))) # incident seroconversion
  
  traj$Z_all = 0.26*traj$Z_sp + 0.29*traj$Z_wt + 0.24*traj$Z_au + 0.20*traj$Z_wt
  traj$conv <- exp(traj$Z_all)
  
  return(traj)
  
}


# TRAJECTORY SIMULATION  ------------------------------------------------------------

# This function draws n samples from the posterior to calculate the posterior predictive uncertainty (95%)

maketrajsim <- function(trace, theta, age, model, inits, ndraw, data) {
  
  #Draw n fitted parameter vectors theta from the MCMC object
  sample <- getSample(trace, parametersOnly = TRUE, thin=1, numSamples=ndraw) #trace is a sampler, parametersOnly = T means that likelihood, posterior and prior values are not provided in the output, thin = thinning parameter
  
  traj.rep <- adply(.data=sample, .margins=1, .progress="text", .parallel=F, .fun=function(x) { #split sample by 1 = rows and apply function
    
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
# B = rate of waning maternal immunity
# C = contact parameter
# D = daycare parameter
theta <- c(P=0.00001, M=0.02002, A=0.00003, W=0.00004, B = 0.01, 
            D = 2) # these are just random values, to be fitted

# INITS ---------------------------------------------------------

#Should this add up to 1?
inits <- c(M_sp=(1-2*1e-12), M_sm = (1-2*1e-12), M_au= (1-2*1e-12), M_wt = (1-2*1e-12),
           S_sp=1e-12, S_sm=1e-12, S_au=1e-12, S_wt=1e-12, 
           Z_sp = 1e-12, Z_sm=1e-12, Z_au=1e-12, Z_wt=1e-12)
           #Z_all = 1e-12) # initial conditions for the states (as proportions)
# --> since we integrate on a log-scale, the initial conditions cannot be 0 (not defined on a log-scale)

# SIMULATION TIME  ---------------------------------------------------------
data <- arrange(data, agemid)
agepred <- data$agemid

# TEST MODEL  --------------------------------------------------------
# Notes: lambda should vary over time, mu shouldn't
test <- model(theta, agepred, inits, data)
ggplot(test) + geom_line(aes(x=time, y=conv))
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
  
  ll_all <- sum(dbinom(x=nconv,
                   size=N,
                   prob=prob,
                   log=TRUE), na.rm=TRUE)
  
  # Spring birth cohort
  spring_cohort <- subset(data, season_birth == "Spring")
  nconv_sp <- spring_cohort$nconv[!is.na(spring_cohort$nconv)] # n seroconverted at each age point  (data)
  N_sp <- spring_cohort$N[!is.na(spring_cohort$nconv)] # total N at each age point  (data) 
  prob_sp <- traj$conv_spring[traj$time %in% spring_cohort$agemid] # proportion seroconverted at each age point (model output)
  
  ll_sp <- sum(dbinom(x=nconv_sp,
                      size=N_sp,
                      prob=prob_sp,
                      log=TRUE), na.rm=TRUE)
  
  #Summer birth cohort
  summer_cohort <- subset(data, season_birth == "Summer")
  nconv_sm <- summer_cohort$nconv[!is.na(summer_cohort$nconv)] # n seroconverted at each age point  (data)
  N_sm <- summer_cohort$N[!is.na(summer_cohort$nconv)] # total N at each age point  (data) 
  prob_sm <- traj$conv_summer[traj$time %in% summer_cohort$agemid] # proportion seroconverted at each age point (model output)
  
  ll_sm <- sum(dbinom(x=nconv_sm,
                      size=N_sm,
                      prob=prob_sm,
                      log=TRUE), na.rm=TRUE)
  
  #Autumn birth cohort
  autumn_cohort <- subset(data, season_birth == "Autumn")
  nconv_au <- autumn_cohort$nconv[!is.na(autumn_cohort$nconv)] # n seroconverted at each age point  (data)
  N_au <- autumn_cohort$N[!is.na(autumn_cohort$nconv)] # total N at each age point  (data) 
  prob_au <- traj$conv_autumn[traj$time %in% autumn_cohort$agemid] # proportion seroconverted at each age point (model output)
  
  ll_au <- sum(dbinom(x=nconv_au,
                      size=N_au,
                      prob=prob_au,
                      log=TRUE), na.rm=TRUE)
  
  # Winter birth cohort
  winter_cohort <- subset(data, season_birth == "Winter")
  nconv_wt <- winter_cohort$nconv[!is.na(winter_cohort$nconv)] # n seroconverted at each age point  (data)
  N_wt <- winter_cohort$N[!is.na(winter_cohort$nconv)] # total N at each age point  (data) 
  prob_wt <- traj$conv_winter[traj$time %in% winter_cohort$agemid] # proportion seroconverted at each age point (model output)
  
  ll_wt <- sum(dbinom(x=nconv_wt,
                      size=N_wt,
                      prob=prob_wt,
                      log=TRUE), na.rm=TRUE)
  
  ll = sum(ll_sp, ll_sm, ll_au, ll_wt)
  
  return(ll_all)
  
} 

# Test function
loglik(theta, agepred, data, model, inits)

# Wrapper for BT: loglik can only take the fitted parameters as argument
loglik_wrapper <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(loglik(theta = parX,
                age = agepred, 
                data = data,
                model = match.fun(model),
                inits = inits))
} 

# FITTING -------------------------------------------

# Estimated params
estpars <- c("P", "M", "A", "W", "B", "D") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params


# Priors
lower = c(P=0, M=0, A=0, W=0, B = 0 , D = 0)
upper = c(P=0.1, M=0.1, A=0.1, W = 0.1, B = 0.2 , D = 10)

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
plot(trace, parametersOnly = TRUE, start =nburn)

# check convergence and correlations
gelmanDiagnostics(trace, plot=TRUE, start=nburn)
correlationPlot(getSample(trace, parametersOnly = TRUE, coda=TRUE, start=nburn), density="smooth", thin=50)
marginalPlot(trace, prior=T, singlePanel=T, start=nburn, nDrawsPrior = 1000)

# remove burn-in for trajsim simulation
tracefinal <- getSample(trace, parametersOnly = TRUE, coda=TRUE, start=nburn)
plot(tracefinal)
effectiveSize(tracefinal)

# Posterior summary
summary(tracefinal)

# save the trace
saveRDS(trace, "new_eq_and_log_DEzs_trace_seasonal_FOI_contacts_one_daycare_and_w.rds")


# POSTPROCESSING AND RESULTS -----------------------------------

# Calculate simulated trajectory quantiles
trajsim <- maketrajsim(tracefinal, theta, agepred, model, inits, 1000)
trajquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles) <- c("agemid", "low95", "median", "up95")

lambda_spquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"lambda_sp"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_spquantiles) <- c("agemid", "low95", "median", "up95")

lambda_smquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"lambda_sm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_smquantiles) <- c("agemid", "low95", "median", "up95")

lambda_auquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"lambda_au"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_auquantiles) <- c("agemid", "low95", "median", "up95")

lambda_wtquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"lambda_wt"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_wtquantiles) <- c("agemid", "low95", "median", "up95")

Pquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"P"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Pquantiles) <- c("agemid", "low95", "median", "up95")

Mquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"M"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Mquantiles) <- c("agemid", "low95", "median", "up95")

Aquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"A"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Aquantiles) <- c("agemid", "low95", "median", "up95")

Wquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"W"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Wquantiles) <- c("agemid", "low95", "median", "up95")

Cquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"C"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Cquantiles) <- c("agemid", "low95", "median", "up95")

Dquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"D"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Dquantiles) <- c("agemid", "low95", "median", "up95")

wquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"mu"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(wquantiles) <- c("agemid", "low95", "median", "up95")

trajquantiles_sp <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_spring"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sp) <- c("agemid", "low95", "median", "up95")

trajquantiles_sm <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_summer"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sm) <- c("agemid", "low95", "median", "up95")

trajquantiles_au <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_autumn"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_au) <- c("agemid", "low95", "median", "up95")

trajquantiles_wt <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_winter"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_wt) <- c("agemid", "low95", "median", "up95")

spring_conv_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_spring"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(spring_conv_quantiles) <- c("agemid", "low95", "median", "up95")

summer_conv_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_summer"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(summer_conv_quantiles) <- c("agemid", "low95", "median", "up95")

autumn_conv_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_autumn"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(autumn_conv_quantiles) <- c("agemid", "low95", "median", "up95")

winter_conv_quantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv_winter"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(winter_conv_quantiles) <- c("agemid", "low95", "median", "up95")

# Plot fit and FOI
fit <- ggplot() + theme_bw() + ggtitle("Model fit") +
  geom_point(data=data_no_season, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=data_no_season, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit

fit_season <- ggplot() + theme_bw() + ggtitle("Model fit on data stratified by season") +
  geom_point(data=data, aes(x=agemid, y=seroprev_mean, color = season_birth, shape = factor(visitnursery_child))) +
  geom_linerange(data=data, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, color = season_birth, linetype = factor(visitnursery_child))) +
  geom_ribbon(data=trajquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") + 
  labs(color = "Season of birth", shape = "Day-care attendance", linetype = "Day-care attendance")

fit_season

fit_sp <- ggplot() + theme_bw() + ggtitle("Model fit on spring birth cohort") +
  geom_point(data=spring.df, aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data=spring.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data=trajquantiles_sp, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_sp, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs (colour = "Day-care attendance")

fit_sp

fit_sp2 <- ggplot() + theme_bw() + ggtitle("Model fit on spring birth cohort") +
  geom_point(data=spring_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=spring_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_sp, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_sp, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit_sp2

fit_sm <- ggplot() + theme_bw() + ggtitle("Model fit on summer birth cohort") +
  geom_point(data=summer.df, aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data=summer.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data=trajquantiles_sm, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_sm, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs(colour = "Day-care attendance")

fit_sm

fit_sm2 <- ggplot() + theme_bw() + ggtitle("Model fit on summer birth cohort") +
  geom_point(data=summer_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=summer_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_sm, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_sm, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit_sm2

fit_au <- ggplot() + theme_bw() + ggtitle("Model fit on autumn birth cohort") +
  geom_point(data=autumn.df, aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data=autumn.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data=trajquantiles_au, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_au, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs(colour = "Day-care attendance")

fit_au

fit_au2 <- ggplot() + theme_bw() + ggtitle("Model fit on autumn birth cohort") +
  geom_point(data=autumn_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=autumn_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_au, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_au, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit_au2

fit_wt <- ggplot() + theme_bw() + ggtitle("Model fit on winter birth cohort") +
  geom_point(data=winter.df, aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data=winter.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data=trajquantiles_wt, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_wt, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs(colour = "Day-care attendance")

fit_wt

fit_wt2 <- ggplot() + theme_bw() + ggtitle("Model fit on winter birth cohort") +
  geom_point(data=winter_no_daycare, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=winter_no_daycare, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_wt, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_wt, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit_wt2

lambda_sp <- ggplot() + theme_bw() + ggtitle("FOI for the spring birth cohort (no daycare)") +
  geom_ribbon(data=lambda_spquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_spquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_sp

lambda_sm <- ggplot() + theme_bw() + ggtitle("FOI for the summer birth cohort (no daycare)") +
  geom_ribbon(data=lambda_smquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_smquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_sm

lambda_au <- ggplot() + theme_bw() + ggtitle("FOI for the autumn birth cohort (no daycare)") +
  geom_ribbon(data=lambda_auquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_auquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_au

lambda_wt <- ggplot() + theme_bw() + ggtitle("FOI for the winter birth cohort (no daycare)") +
  geom_ribbon(data=lambda_wtquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_wtquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_wt

w <- ggplot() + theme_bw() + ggtitle("Rate of maternal immunity waning ") +
  geom_ribbon(data=wquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=wquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("Maternal immunity waning rate") 
w

write.csv(lambda_spquantiles, "add_DEzs_lambda_sp_ll_noc.csv")
write.csv(lambda_smquantiles, "add_DEzs_lambda_sm_ll_noc.csv")
write.csv(lambda_auquantiles, "add_DEzs_lambda_au_ll_noc.csv")
write.csv(lambda_auquantiles, "add_DEzs_lambda_wt_ll_noc.csv")

write.csv(Pquantiles, "add_DEzs_Pquantiles_ll_noc.csv")
write.csv(Mquantiles, "add_DEzs_Mquantiles_ll_noc.csv")
write.csv(Aquantiles, "add_DEzs_Aquantiles_ll_noc.csv")
write.csv(Wquantiles, "add_DEzs_Wquantiles_ll_noc.csv")
write.csv(Cquantiles, "add_DEzs_Cquantiles_ll_noc.csv")
write.csv(Dquantiles, "add_DEzs_Dquantiles_ll_noc.csv")
write.csv(wquantiles, "add_DEzs_Immunity_quantiles_ll_noc.csv")

write.csv(trajquantiles, "add_DEzs_trajquantiles_ll_noc.csv")
write.csv(trajsim, "add_DEzs_trajsim_ll_noc.csv")
write.csv(spring_conv_quantiles, "add_DEzs_spring_conv_ll_noc.csv")
write.csv(summer_conv_quantiles, "add_DEzs_summer_conv_ll_noc.csv")
write.csv(autumn_conv_quantiles, "add_DEzs_autumn_conv_ll_noc.csv")
write.csv(winter_conv_quantiles, "add_DEzs_winter_conv_ll_noc.csv")
