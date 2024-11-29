################################################################
# RSV seroconversion MSc project
# Adding daycare attendance as a multiplication
# Author: Julia Mayer
# Last updated: 23.08.2022
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
  mutate(
    Birth_mo = birthday %>% month(),
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
data_no_daycare <- data %>% group_by(agegrp, season_birth) %>% 
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
data_no_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_season$nconv, data_no_season$N, method="exact")[,c("mean","lower","upper")]
data_no_daycare[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_daycare$nconv, data_no_daycare$N, method="exact")[,c("mean","lower","upper")]

# Plot the whole data frame (difficult to read)
data %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = season_birth, 
                 shape = factor(visitnursery_child))) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, 
                    ymax = seroprev_up95, colour = season_birth,
                    linetype = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(colour = "season of birth", shape = "Day-care attendance", linetype = "Day-care attendance")

# Plot data grouped by age only
data_no_season %>% ggplot() + 
  geom_point(aes(x = agemid, y = seroprev_mean)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (years)") +
  scale_x_continuous(breaks = c(0,365,730,1095,1460,1825), labels = c(0, 1, 2, 3, 4 ,5))

# Plot data grouped by age and season of birth
data_no_daycare %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = season_birth)) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = season_birth)) +
  ylab("Proportion seroconverted") + xlab("age (days)")  + labs (colour = "Season of birth")

# Plot by season
spring_df <- data %>% subset(season_birth == 'Spring')
spring_no_daycare <- data %>% subset(season_birth == 'Spring')
spring_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in spring", colour = "Day-care")

summer_df <- data %>% subset(season_birth == 'Summer')
summer_no_daycare <- data_no_daycare %>% subset(season_birth == 'Summer')
summer_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in summer", colour = "Day-care attendance")

autumn_df <- data %>% subset(season_birth == 'Autumn')
autumn_no_daycare <- data_no_daycare %>% subset(season_birth == 'Autumn')
autumn_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in autumn", colour = "Day-care attendance")

winter_df <- data %>% subset(season_birth == 'Winter')
winter_no_daycare <- data_no_daycare %>% subset(season_birth == 'Winter')
winter_df %>% ggplot() +
  geom_point(aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title = "Proportion seroconverted born in winter", colour = "Day-care attendance")

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
    
    # waning maternal immunity, same for all children
    mu = param[["B"]] 
    
    # states 
    # proportion with maternal immunity
    M_sp_d = exp(state[1]) # born in spring attending day-care
    M_sm_d = exp(state[2]) # born in summer attending day-care
    M_au_d = exp(state[3]) # born in autumn attending day-care
    M_wt_d = exp(state[4]) # born in winter attending day-care
    M_sp_n = exp(state[5]) # born in spring not attending day-care
    M_sm_n = exp(state[6]) # born in summer not attending day-care
    M_au_n = exp(state[7]) # born in autumn not attending day-care
    M_wt_n = exp(state[8]) # born in winter not attending day-care
    # Susceptible
    S_sp_d = exp(state[9]) # susceptible born in spring attending day-care
    S_sm_d = exp(state[10]) # susceptible born in summer attending day-care
    S_au_d = exp(state[11]) # susceptible born in autumn attending day-care
    S_wt_d = exp(state[12]) # susceptible born in winter attending day-care
    S_sp_n = exp(state[13]) # susceptible born in spring not attending day-care
    S_sm_n = exp(state[14]) # susceptible born in summer not attending day-care
    S_au_n = exp(state[15]) # susceptible born in autumn not attending day-care
    S_wt_n = exp(state[16]) # susceptible born in winter not attending day-care
    #Seroconverted
    Z_sp_d = exp(state[17]) # seroconverted after infection born in spring attending day-care
    Z_sm_d = exp(state[18]) # seroconverted after infection born in summer attending day-care
    Z_au_d = exp(state[19]) # seroconverted after infection born in autumn attending day-care
    Z_wt_d = exp(state[20]) # seroconverted after infection born in winter attending day-care
    Z_sp_n = exp(state[21]) # seroconverted after infection born in spring not attending day-care
    Z_sm_n = exp(state[22]) # seroconverted after infection born in summer not attending day-care
    Z_au_n = exp(state[23]) # seroconverted after infection born in autumn not attending day-care
    Z_wt_n = exp(state[24]) # seroconverted after infection born in winter not attending day-care
    
    # changes in states
    dM_sp_d = -mu*M_sp_d
    dM_sm_d = -mu*M_sm_d
    dM_au_d = -mu*M_au_d
    dM_wt_d = -mu*M_wt_d
    dM_sp_n = -mu*M_sp_n
    dM_sm_n = -mu*M_sm_n
    dM_au_n = -mu*M_au_n
    dM_wt_n = -mu*M_wt_n
    
    if (age >= 30.41*9){ # assume that all children who go to day-care start at 9 months
      dS_sp_d = + mu*M_sp_d - param[["D"]]*lambda_sp*S_sp_d
      dS_sm_d = + mu*M_sm_d - param[["D"]]*lambda_sm*S_sm_d
      dS_au_d = + mu*M_au_d - param[["D"]]*lambda_au*S_au_d 
      dS_wt_d = + mu*M_wt_d - param[["D"]]*lambda_wt*S_wt_d
      dZ_sp_d = + param[["D"]]*lambda_sp*S_sp_d
      dZ_sm_d = + param[["D"]]*lambda_sm*S_sm_d
      dZ_au_d = + param[["D"]]*lambda_au*S_au_d
      dZ_wt_d = + param[["D"]]*lambda_wt*S_wt_d
    }else{
      dS_sp_d = + mu*M_sp_d - lambda_sp*S_sp_d
      dS_sm_d = + mu*M_sm_d - lambda_sm*S_sm_d
      dS_au_d = + mu*M_au_d - lambda_au*S_au_d 
      dS_wt_d = + mu*M_wt_d - lambda_wt*S_wt_d
      dZ_sp_d = + lambda_sp*S_sp_d
      dZ_sm_d = + lambda_sm*S_sm_d
      dZ_au_d = + lambda_au*S_au_d
      dZ_wt_d = + lambda_wt*S_wt_d
    }
    
    dS_sp_n = + mu*M_sp_n - lambda_sp*S_sp_n
    dS_sm_n = + mu*M_sm_n - lambda_sm*S_sm_n
    dS_au_n = + mu*M_au_n - lambda_au*S_au_n
    dS_wt_n = + mu*M_wt_n - lambda_wt*S_wt_n
    dZ_sp_n = + lambda_sp*S_sp_n
    dZ_sm_n = + lambda_sm*S_sm_n
    dZ_au_n = + lambda_au*S_au_n
    dZ_wt_n = + lambda_wt*S_wt_n
    
    
    return(list(c(dM_sp_d/M_sp_d,dM_sm_d/M_sm_d, dM_au_d/M_au_d,dM_wt_d/M_wt_d,
                  dM_sp_n/M_sp_n,dM_sm_n/M_sm_n, dM_au_n/M_au_n,dM_wt_n/M_wt_n,
                  dS_sp_d/S_sp_d, dS_sm_d/S_sm_d, dS_au_d/S_au_d,dS_wt_d/S_wt_d,
                  dS_sp_n/S_sp_n, dS_sm_n/S_sm_n, dS_au_n/S_au_n,dS_wt_n/S_wt_n,
                  dZ_sp_d/Z_sp_d, dZ_sm_d/Z_sm_d, dZ_au_d/Z_au_d, dZ_wt_d/Z_wt_d,
                  dZ_sp_n/Z_sp_n, dZ_sm_n/Z_sm_n, dZ_au_n/Z_au_n, dZ_wt_n/Z_wt_n),
                lambda_sp = lambda_sp, lambda_sm = lambda_sm, 
                lambda_au = lambda_au, lambda_wt = lambda_wt,
                mu=mu))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sp_d=log(inits[["M_sp_d"]]),
                             M_sm_d=log(inits[["M_sm_d"]]),
                             M_au_d=log(inits[["M_au_d"]]),
                             M_wt_d=log(inits[["M_wt_d"]]),
                             M_sp_n=log(inits[["M_sp_n"]]),
                             M_sm_n=log(inits[["M_sm_n"]]),
                             M_au_n=log(inits[["M_au_n"]]),
                             M_wt_n=log(inits[["M_wt_n"]]),
                             S_sp_d=log(inits[["S_sp_d"]]),
                             S_sm_d=log(inits[["S_sm_d"]]),
                             S_au_d=log(inits[["S_au_d"]]),
                             S_wt_d=log(inits[["S_wt_d"]]),
                             S_sp_n=log(inits[["S_sp_n"]]),
                             S_sm_n=log(inits[["S_sm_n"]]),
                             S_au_n=log(inits[["S_au_n"]]),
                             S_wt_n=log(inits[["S_wt_n"]]),
                             Z_sp_d=log(inits[["Z_sp_d"]]),
                             Z_sm_d=log(inits[["Z_sm_d"]]),
                             Z_au_d=log(inits[["Z_au_d"]]),
                             Z_wt_d=log(inits[["Z_wt_d"]]),
                             Z_sp_n=log(inits[["Z_sp_n"]]),
                             Z_sm_n=log(inits[["Z_sm_n"]]),
                             Z_au_n=log(inits[["Z_au_n"]]),
                             Z_wt_n=log(inits[["Z_wt_n"]])),
                         times = age, 
                         func = catalytic, 
                         parms = theta, 
                         data = data,
                         method ="lsoda",
                         verbose = F))
  
  # Cumulative seroconversion if day-care == True
  traj$conv_spring_d <- exp(traj$Z_sp_d) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring_d <- c(inits[["Z_sp_d"]], diff(exp(traj$Z_sp_d))) # incident seroconversion
  traj$conv_summer_d <- exp(traj$Z_sm_d) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer_d <- c(inits[["Z_sm_d"]], diff(exp(traj$Z_sm_d))) # incident seroconversion
  traj$conv_autumn_d <- exp(traj$Z_au_d) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn_d <- c(inits[["Z_au_d"]], diff(exp(traj$Z_au_d))) # incident seroconversion
  traj$conv_winter_d <- exp(traj$Z_wt_d) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter_d <- c(inits[["Z_wt_d"]], diff(exp(traj$Z_wt_d))) # incident seroconversion
  
  traj$Z_all_d <- 0.26*traj$Z_sp_d + # total seroconverted is the sum by birth cohort
    0.29*traj$Z_sm_d +               # we scale by the proportion of children born in each season
    0.24*traj$Z_au_d + 
    0.20*traj$Z_wt_d
  traj$conv_d <- exp(traj$Z_all_d)
  
  # Cumulative seroconversion if day-care == False
  traj$conv_spring_n <- exp(traj$Z_sp_n) # cumulative seroconversion in spring cohort (=observed state)
  traj$inc_spring_n <- c(inits[["Z_sp_n"]], diff(exp(traj$Z_sp_d))) # incident seroconversion
  traj$conv_summer_n <- exp(traj$Z_sm_n) # cumulative seroconversion in summer cohort (=observed state)
  traj$inc_summer_n <- c(inits[["Z_sm_n"]], diff(exp(traj$Z_sm_d))) # incident seroconversion
  traj$conv_autumn_n <- exp(traj$Z_au_n) # cumulative seroconversion in autumn cohort (=observed state)
  traj$inc_autumn_n <- c(inits[["Z_au_n"]], diff(exp(traj$Z_au_d))) # incident seroconversion
  traj$conv_winter_n <- exp(traj$Z_wt_n) # cumulative seroconversion in winter cohort (=observed state)
  traj$inc_winter_n <- c(inits[["Z_wt_n"]], diff(exp(traj$Z_wt_d))) # incident seroconversion
  
  traj$Z_all_n <- 0.26*traj$Z_sp_n + 0.29*traj$Z_sm_n + 0.24*traj$Z_au_n + 0.20*traj$Z_wt_n
  traj$conv_n <- exp(traj$Z_all_n)
  
  # Overall cumulative seroconversion (0.3944 go to day-care)
  traj$Z_sp_tot <- 0.3944*traj$Z_sp_d + (1-0.3944)*traj$Z_sp_n
  traj$spring_conv <- exp(traj$Z_sp_tot)
  traj$Z_sm_tot <- 0.3944*traj$Z_sm_d + (1-0.3944)*traj$Z_sm_n
  traj$conv_summer <- exp(traj$Z_sm_tot)
  traj$Z_au_tot <- 0.3944*traj$Z_au_d + (1-0.3944)*traj$Z_au_n
  traj$conv_autumn <- exp(traj$Z_au_tot)
  traj$Z_wt_tot <- 0.3944*traj$Z_wt_d + (1-0.3944)*traj$Z_wt_n
  traj$conv_winter <- exp(traj$Z_wt_tot)
  traj$Z_all_tot <- 0.3944*traj$Z_all_d + (1-0.3944)*traj$Z_all_n
  traj$conv_tot <- exp(traj$Z_all_tot) 
  
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
# B = rate of waning maternal immunity
# C = contact parameter
# D = daycare parameter
theta <- c(P = 0.00001, M = 0.02002, A = 0.00003, W = 0.00004, B = 0.01, 
           C = 0.02, D = 2) # these are just random values, to be fitted

# INITS ---------------------------------------------------------
inits <- c(M_sp_d=(1-2*1e-12), M_sm_d = (1-2*1e-12), M_au_d = (1-2*1e-12), M_wt_d = (1-2*1e-12),
           M_sp_n=(1-2*1e-12), M_sm_n = (1-2*1e-12), M_au_n = (1-2*1e-12), M_wt_n = (1-2*1e-12),
           S_sp_d=1e-12, S_sm_d=1e-12, S_au_d=1e-12, S_wt_d=1e-12,
           S_sp_n=1e-12, S_sm_n=1e-12, S_au_n=1e-12, S_wt_n=1e-12,
           Z_sp_d = 1e-12, Z_sm_d=1e-12, Z_au_d=1e-12, Z_wt_d=1e-12,
           Z_sp_n = 1e-12, Z_sm_n=1e-12, Z_au_n=1e-12, Z_wt_n=1e-12)
# initial conditions for the states (as proportions)
# --> since we integrate on a log-scale, the initial conditions cannot be 0 (not defined on a log-scale)

# SIMULATION TIME  ---------------------------------------------------------
data <- arrange(data, agemid)
agepred <- data$agemid

# TEST MODEL  --------------------------------------------------------
# Notes: lambda should vary over time, mu shouldn't
test <- model(theta, agepred, inits, data)
ggplot(test) + geom_line(aes(x=time, y=conv_n))
ggplot(test) + geom_line(aes(x=time, y=conv_d))
ggplot(test) + geom_line(aes(x=time, y = conv_tot))

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
  prob <- traj$conv_tot[traj$time %in% data$agemid] # proportion seroconverted at each age point (model output)
  
  ll_all <- sum(dbinom(x = nconv,
                       size = N,
                       prob = prob,
                       log = TRUE), na.rm = TRUE)
  
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
estpars <- c("P", "M", "A", "W", "B", "C", "D") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params


# Priors
lower = c(P = 0, M = 0, A = 0, W = 0, B = 0 , C = 0, D = 0)
upper = c(P = 0.1, M = 0.1, A = 0.1, W = 0.1, B = 0.2 , C = 0.1, D = 10)

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
saveRDS(trace, "strat_model_trace_seasonal_FOI_contacts_one_daycare_and_w.rds")


# POSTPROCESSING AND RESULTS -----------------------------------

# Calculate simulated trajectory quantiles
trajsim <- maketrajsim(tracefinal, theta, agepred, model, inits, 1000)

# Proportion converted overall
trajquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                             function(x) quantile(x[,"conv_tot"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles) <- c("agemid", "low95", "median", "up95")

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

Dquantiles <- plyr::ddply(.data = trajsim, .variables = "time", function(x) quantile(x[,"D"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(Dquantiles) <- c("agemid", "low95", "median", "up95")

wquantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                          function(x) quantile(x[,"mu"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(wquantiles) <- c("agemid", "low95", "median", "up95")

# Proportion seroconverted by season of birth and day-care attendance
trajquantiles_sp_d <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_spring_d"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sp_d) <- c("agemid", "low95", "median", "up95")

trajquantiles_sp_n <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_spring_n"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sp_n) <- c("agemid", "low95", "median", "up95")

trajquantiles_sm_d <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_summer_d"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sm_d) <- c("agemid", "low95", "median", "up95")

trajquantiles_sm_n <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_summer_n"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sm_n) <- c("agemid", "low95", "median", "up95")

trajquantiles_au_d <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_autumn_d"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_au_d) <- c("agemid", "low95", "median", "up95")

trajquantiles_au_n <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_autumn_n"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_au_n) <- c("agemid", "low95", "median", "up95")

trajquantiles_wt_d <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_winter_d"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_wt_d) <- c("agemid", "low95", "median", "up95")

trajquantiles_wt_n <- plyr::ddply(.data = trajsim, .variables = "time", 
                                  function(x) quantile(x[,"conv_winter_n"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_wt_n) <- c("agemid", "low95", "median", "up95")

# Proportion seroconverted by season of birth not accounting for day-care attendance
spring_conv_quantiles <- plyr::ddply(.data = trajsim, .variables = "time", 
                                     function(x) quantile(x[,"spring_conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
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

# Plot fit and FOI
# Whole dataset
fit <- ggplot() + theme_bw() + ggtitle("Model fit") +
  geom_point(data = data_no_season, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = data_no_season, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = trajquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = trajquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit

# Stratified by season of birth
fit_season <- ggplot() + theme_bw() + ggtitle("Model fit on data stratified by season") +
  geom_point(data = data, aes(x = agemid, y = seroprev_mean, color = season_birth, shape = factor(visitnursery_child))) +
  geom_linerange(data = data, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, color = season_birth, linetype = factor(visitnursery_child))) +
  geom_ribbon(data = trajquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = trajquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") + 
  labs(color = "Season of birth", shape = "Day-care attendance", linetype = "Day-care attendance")

fit_season

#Spring birth cohort by day-care attendance
fit_sp <- ggplot() + theme_bw() + ggtitle("Model fit on spring birth cohort") +
  geom_point(data = spring_df, aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data = spring_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data = spring_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = spring_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs (colour = "Day-care attendance")

fit_sp

#Spring birth cohort w/o day-care attendance
fit_sp2 <- ggplot() + theme_bw() + ggtitle("Model fit on spring birth cohort") +
  geom_point(data = spring_no_daycare, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = spring_no_daycare, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = spring_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = spring_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))
fit_sp2

#Summmer birth cohort by day-care attendance
fit_sm <- ggplot() + theme_bw() + ggtitle("Model fit on summer birth cohort") +
  geom_point(data = summer_df, aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data = summer_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data = summer_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = summer_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs(colour = "Day-care attendance")

fit_sm

#Summer birth cohort w/o day-care attendance
fit_sm2 <- ggplot() + theme_bw() + ggtitle("Model fit on summer birth cohort") +
  geom_point(data = summer_no_daycare, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = summer_no_daycare, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = summer_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = summer_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit_sm2

#Autumn birth cohort by day-care attendance
fit_au <- ggplot() + theme_bw() + ggtitle("Model fit on autumn birth cohort") +
  geom_point(data = autumn_df, aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data = autumn_df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data = autumn_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = autumn_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs(colour = "Day-care attendance")

fit_au

#Autumn birth cohort w/o day-care attendance
fit_au2 <- ggplot() + theme_bw() + ggtitle("Model fit on autumn birth cohort") +
  geom_point(data = autumn_no_daycare, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = autumn_no_daycare, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = autumn_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = autumn_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit_au2

#Winter birth cohort by day-care attendance
fit_wt <- ggplot() + theme_bw() + ggtitle("Model fit on winter birth cohort") +
  geom_point(data = winter.df, aes(x = agemid, y = seroprev_mean, colour = factor(visitnursery_child))) +
  geom_linerange(data = winter.df, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95, colour = factor(visitnursery_child))) +
  geom_ribbon(data = winter_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = winter_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("proportion seroconverted") + labs(colour = "Day-care attendance")

fit_wt

#Winter birth cohort w/o day-care attendance
fit_wt2 <- ggplot() + theme_bw() + ggtitle("Model fit on winter birth cohort") +
  geom_point(data = winter_no_daycare, aes(x = agemid, y = seroprev_mean)) +
  geom_linerange(data = winter_no_daycare, aes(x = agemid, ymin = seroprev_low95, ymax = seroprev_up95)) +
  geom_ribbon(data = winter_conv_quantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = winter_conv_quantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (years)") + ylab("proportion seroconverted") +
  scale_x_continuous(breaks=c(0,182.5, 365, 547.5, 730, 912.5, 1095, 1277.5, 1460, 1642.5, 1825), 
                     labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5 ,5))

fit_wt2

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

# Waning rate of maternal immunity (should be constant)
w <- ggplot() + theme_bw() + ggtitle("Rate of maternal immunity waning ") +
  geom_ribbon(data = wquantiles, aes(x = agemid, ymin = low95, ymax = up95), fill = "red", alpha = 0.3) +
  geom_line(data = wquantiles, aes(x = agemid, y = median), color = "red") +
  xlab("age (days)") + ylab("Maternal immunity waning rate") 
w

# Saving the useful files ------------------------------------------------------------
write.csv(lambda_spquantiles, "strat_lambda_sp.csv")
write.csv(lambda_smquantiles, "strat_lambda_sm.csv")
write.csv(lambda_auquantiles, "strat_lambda_au.csv")
write.csv(lambda_wtquantiles, "strat_lambda_wt.csv")

write.csv(Pquantiles, "strat_Pquantiles.csv")
write.csv(Mquantiles, "strat_Mquantiles.csv")
write.csv(Aquantiles, "strat_Aquantiles.csv")
write.csv(Wquantiles, "strat_Wquantiles.csv")
write.csv(Cquantiles, "strat_Cquantiles.csv")
write.csv(Dquantiles, "strat_Dquantiles.csv")
write.csv(wquantiles, "strat_Immunity_quantiles.csv")

write.csv(trajquantiles, "strat_trajquantiles.csv")
write.csv(trajsim, "strat_trajsim.csv")
write.csv(spring_conv_quantiles, "strat_spring_conv_tot.csv")
write.csv(summer_conv_quantiles, "strat_ummer_conv_tot.csv")
write.csv(autumn_conv_quantiles, "strat_autumn_conv_tot.csv")
write.csv(winter_conv_quantiles, "strat_winter_conv_tot.csv")
