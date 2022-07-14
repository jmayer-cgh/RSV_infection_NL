################################################################
# RSV seroconversion MSc project
# Adding season of birth
# Author: Julia Mayer
# Last updated: 03.07.2022
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
  )

data_no_season <- data %>% dplyr::group_by(agegrp) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection))

data <- data %>% dplyr::group_by(agegrp, season_birth) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection)) # n seroconverted in age group

# Calculate seroprevalence and binomial confidence intervals
data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data$nconv, data$N, method="exact")[,c("mean","lower","upper")]
data_no_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_season$nconv, data_no_season$N, method="exact")[,c("mean","lower","upper")]


# Get proportion of children born in each season for each age group
get_proportions <- function (dframe){
  agegrp <- c(0)
  season_birth <- c("Pomme")
  agemid <- c (0)
  N <- c(0)
  nconv <- c(0)
  seroprev_mean <- c(0)
  seroprev_low95 <- c(0)
  seroprev_up95 <- c(0)
  prop <- c(0)
  df <- data.frame(agegrp, season_birth, agemid, N, nconv, seroprev_mean, seroprev_low95, seroprev_up95, prop)
  
  groups <- unique(dframe$agegrp)
  
  for (gp in groups){
    subs <- subset (dframe, agegrp == gp)
    subs$agegrp <- as.double(subs$agegrp)
    subs$prop <- subs$N/sum(subs$N)
    df_list <- list(df, subs)
    df <- df_list %>% reduce(full_join)
  }
  return (df)
}

prop.df<- get_proportions(data)
prop.df <- prop.df[-c(1), ]

# Plot the proportions of children born in each season

p <- ggplot(prop.df, aes(x = agegrp, y=prop, fill = season_birth)) + 
  geom_bar(stat="identity", position = "dodge") +
  labs(title="Proportion born in each season by age group",
                                 x ="Age group", y = "Proportion") +
  scale_fill_discrete(name="Season")
p

born_spring <- prop.df%>%filter(season_birth == "Spring")
sp <- ggplot(born_spring, aes(x = agegrp, y=prop)) + 
  geom_bar(stat="identity", fill = "lightgreen") +
  labs(title="Proportion born in spring by age group",
       x ="Age group", y = "Proportion") 
sp

born_summer <- prop.df%>%filter(season_birth == "Summer")
sm <- ggplot(born_summer, aes(x = agegrp, y=prop)) + 
  geom_bar(stat="identity", fill = "lightblue") +
  labs(title="Proportion born in summer by age group",
       x ="Age group", y = "Proportion") 
sm

born_autumn <- prop.df%>%filter(season_birth == "Autumn")
au <- ggplot(born_autumn, aes(x = agegrp, y=prop)) + 
  geom_bar(stat="identity", fill = "red") +
  labs(title="Proportion born in autumn by age group",
       x ="Age group", y = "Proportion") 
au

born_winter <- prop.df%>%filter(season_birth == "Winter")
wt <- ggplot(born_winter, aes(x = agegrp, y=prop)) + 
  geom_bar(stat="identity", fill = "purple") +
  labs(title="Proportion born in winter by age group",
       x ="Age group", y = "Proportion") 
wt

mean_sp <- median(born_spring$prop) #median proportion of children born in spring
mean_sm <- median(born_summer$prop) #median proportion of children born in summer
mean_au <- median(born_autumn$prop) #median proportion of children born in autumn
mean_wt <- median(born_winter$prop) #median proportion of children born in winter

# Plot the whole data frame (difficult to read)
ggplot(data) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = season_birth)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = season_birth)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs(colour = "season of birth")

ggplot(data_no_season) +
  geom_point(aes(x=agemid, y=seroprev_mean)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") 

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

model_sp <- function(theta, age, inits) {
  
  catalytic <- function(age, state, param) {
    
    # FOI / seroconversion rate
    lambda_sp = param[["A"]] + age * 0 #constant FOI for children born in spring, age will be added later
    # waning maternal immunity, same for all children
    mu_sp = param[["B"]] 
    
    # states 
    # proportion with maternal immunity
    M_sp = exp(state[1]) # born in spring
    # Susceptible
    S_sp = exp(state[2]) # susceptible born in spring
    #Seroconverted
    Z_sp = exp(state[3]) # seroconverted after infection born in spring
    
    # changes in states
    dM_sp = - mu_sp*M_sp
    dS_sp = + mu_sp*M_sp - lambda_sp*S_sp 
    dZ_sp = + lambda_sp*S_sp
   
    return(list(c(dM_sp/M_sp,
                  dS_sp/S_sp, 
                  dZ_sp/Z_sp), 
                lambda_sp=lambda_sp,
                mu_sp=mu_sp))
    
    
  }
  
traj <- data.frame(ode(y=c(M_sp=log(inits[["M"]]),
                             S_sp=log(inits[["S"]]),
                             Z_sp=log(inits[["Z"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         method="lsoda",
                         verbose=F))
  

  
  traj$conv <- exp(traj$Z_sp) # cumulative seroconversion (=observed state)
  traj$inc <- c(inits[["Z"]], diff(exp(traj$Z_sp))) # incident seroconversion
  
  return(traj)
  
}

model_sm <- function(theta, age, inits) {
  
  catalytic <- function(age, state, param) {
    
    # FOI / seroconversion rate
    lambda_sm = param[["A"]] + age*0 #constant FOI for children born in summer, age will be added later
    
    # waning maternal immunity, same for all children
    mu_sm = param[["B"]] 
    
    # states 
    M_sm = exp(state[1]) # born in summer
    S_sm = exp(state[2]) # susceptible born in summer
    Z_sm = exp(state[3]) # seroconverted after infection born in summer
    
    # changes in states
    dM_sm = -mu_sm*M_sm
    dS_sm = + mu_sm*M_sm - lambda_sm*S_sm 
    dZ_sm = + lambda_sm*S_sm 
    
    return(list(c(dM_sm/M_sm,dS_sm/S_sm,dZ_sm/Z_sm), lambda_sm = lambda_sm,
                mu_sm=mu_sm))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_sm=log(inits[["M"]]),
                             S_sm=log(inits[["S"]]),
                             Z_sm=log(inits[["Z"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         method="lsoda",
                         verbose=F))
  
  
  
  traj$conv <- exp(traj$Z_sm) # cumulative seroconversion (=observed state)
  traj$inc <- c(inits[["Z"]], diff(exp(traj$Z_sm))) # incident seroconversion
  
  return(traj)
  
}

model_au <- function(theta, age, inits) {
  
  catalytic <- function(age, state, param) {
    
    # FOI / seroconversion rate
    lambda_au = param[["A"]] + age*0 #constant FOI for children born in autumn, age will be added later
    
    # waning maternal immunity, same for all children
    mu_au = param[["B"]] 
    
    # states 
    # proportion with maternal immunity
    M_au = exp(state[1]) # born in autumn
    # Susceptible
    S_au = exp(state[2]) #susceptible born in autumn
    #Seroconverted
    Z_au = exp(state[3]) # seroconverted after infection born in autumn
    
    # changes in states
    dM_au = -mu_au*M_au
    dS_au = + mu_au*M_au - lambda_au*S_au 
    dZ_au = + lambda_au*S_au
    
    return(list(c(dM_au/M_au,
                  dS_au/S_au, dZ_au/Z_au), 
                lambda_au = lambda_au, 
                mu_au=mu_au))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_au=log(inits[["M"]]),
                             S_au=log(inits[["S"]]),
                             Z_au=log(inits[["Z"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         method="lsoda",
                         verbose=F))
  
  
  
  traj$conv <- exp(traj$Z_au) # cumulative seroconversion (=observed state)
  traj$inc <- c(inits[["Z"]], diff(exp(traj$Z_au)))
  return(traj)
  
}

model_wt <- function(theta, age, inits) {
  
  catalytic <- function(age, state, param) {
    
    # FOI / seroconversion rate
    lambda_wt = param [["A"]] + age*0 #constant FOI for children born in winter, age will be added later
    
    # waning maternal immunity, same for all children
    mu_wt = param[["B"]] 
    
    # states 
    # proportion with maternal immunity
    M_wt = exp(state[1]) # born in winter
    # Susceptible
    S_wt = exp(state[2]) #susceptible born in winter
    #Seroconverted
    Z_wt = exp(state[3]) # seroconverted after infection born in winter
    
    # changes in states
    dM_wt = -mu_wt*M_wt
    dS_wt = + mu_wt*M_wt - lambda_wt*S_wt
    dZ_wt = + lambda_wt*S_wt
    
    return(list(c(dM_wt/M_wt,
                 dS_wt/S_wt, 
                  dZ_wt/Z_wt), 
                 lambda_wt = lambda_wt,
                mu_wt=mu_wt))
    
    
  }
  
  traj <- data.frame(ode(y=c(M_wt=log(inits[["M"]]),
                             S_wt=log(inits[["S"]]),
                             Z_wt=log(inits[["Z"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, 
                         method="lsoda",
                         verbose=F))
  
  
  
  traj$conv <- exp(traj$Z_wt) # cumulative seroconversion (=observed state)
  traj$inc <- c(inits[["Z"]], diff(exp(traj$Z_wt))) # incident seroconversion
  
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

# P = mean FOI (proportion infected per day) for children born in Spring
# M = mean FOI for children born in summer
# A = mean FOI for children born in autumn
# W = mean FOI for children born in winter
# B = rate of waning maternal immunity
theta <- c(A=0.02,B = 0.01) # these are just random values, to be fitted

# INITS ---------------------------------------------------------

'#inits <- c(M_sp=0.26*(1-8*1e-12), M_sm = 0.29*(1-8*1e-12), M_au= 0.24*(1-8*1e-12), M_wt = 0.20*(1-8*1e-12),
           S_sp=1e-12, S_sm=1e-12, S_au=1e-12, S_wt=1e-12, 
           Z_sp = 1e-12, Z_sm=1e-12, Z_au=1e-12, Z_wt=1e-12) # initial conditions for the states (as proportions)
# --> since we integrate on a log-scale, the initial conditions cannot be 0 (not defined on a log-scale)
#'

inits <- c(M=1-1e-12-1e-12, S=1e-12, Z=1e-12)

# SIMULATION TIME  ---------------------------------------------------------
data <- arrange(data, agemid)
agepred <- data$agemid
agepred_sp <- spring.df$agemid
agepred_sm <- summer.df$agemid
agepred_au <- autumn.df$agemid
agepred_wt <- winter.df$agemid

# TEST MODEL  --------------------------------------------------------
test_sp <- model_sp(theta, agepred_sp, inits)
test_sm <- model_sm(theta, agepred_sm, inits)
test_au <- model_au(theta, agepred_au, inits)
test_wt <- model_wt(theta, agepred_wt, inits)

conv_all <- c(test_sp$conv, test_sm$conv, test_au$conv, test_wt$conv) #get all converted together
time_all <- c(test_sp$time, test_sm$time, test_au$time, test_wt$time) #get all ages together
test_all <- data.frame(time_all, conv_all)

ggplot(test_all) + geom_line(aes(x=time_all, y=conv_all))
ggplot(test_wt) + geom_line(aes(x=time, y=conv))
ggplot(test_sm) + geom_line(aes(x=time, y=lambda_sm)) #should be constant
ggplot(test_sp) + geom_line(aes(x=time, y=mu_sp))

# LOG LIKELIHOOD FUNCTION ---------------------------------------------------------

loglik <- function(theta, age, data, model, inits) {
  
  traj <- match.fun(model)(theta, age, inits)
  
  nconv <- data$nconv[!is.na(data$nconv)] # n seroconverted at each age point  (data)
  N <- data$N[!is.na(data$nconv)] # total N at each age point  (data) 
  prob <- traj$conv[traj$time %in% data$agemid] # proportion seroconverted at each age point (model output)
  
  ll <- sum(dbinom(x=nconv,
                   size=N,
                   prob=prob,
                   log=TRUE), na.rm=TRUE)
  
  return(ll)
  
} 

# Test function for each season
loglik(theta, agepred_sp, spring.df, model_sp, inits) # spring
loglik(theta, agepred_sm, summer.df, model_sm, inits) # summer
loglik(theta, agepred_au, autumn.df, model_au, inits) # autumn
loglik(theta, agepred_wt, winter.df, model_wt, inits) # winter

# Wrapper for BT: loglik can only take the fitted parameters as argument
loglik_wrapper_sp <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(loglik(theta = parX,
                age = agepred_sp, 
                data = spring.df,
                model = match.fun(model_sp),
                inits = inits))
} 

loglik_wrapper_sm <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(loglik(theta = parX,
                age = agepred_sm, 
                data = summer.df,
                model = match.fun(model_sm),
                inits = inits))
}

loglik_wrapper_au <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(loglik(theta = parX,
                age = agepred_au, 
                data = autumn.df,
                model = match.fun(model_au),
                inits = inits))
}

loglik_wrapper_wt <- function(par) {
  
  parX = theta
  parX[index] = par
  
  return(loglik(theta = parX,
                age = agepred_wt, 
                data = winter.df,
                model = match.fun(model_wt),
                inits = inits))
}
# FITTING -------------------------------------------

# Estimated params
estpars <- c("A", "B") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params


# Priors
lower = c(A=0,B = 0)
upper = c(A=0.1,B = 0.2)

prior <- createUniformPrior(lower=lower[estpars], 
                            upper=upper[estpars])
# MCMC settings
nchains <- 2
cpus <- 1 # or 2 if you want parallel, but it does not seem to be faster?
mcmc_settings <- list(iterations = 2*80000, 
                      nrChains = nchains)
sampler <- "Metropolis"

if (cpus == 1) {
  #spring
  bayesianSetup_sp <- createBayesianSetup(prior = prior,
                                       likelihood = loglik_wrapper_sp,
                                       names = names(theta[index]),
                                       parallel = FALSE)
  
  system.time({trace_sp <- runMCMC(bayesianSetup = bayesianSetup_sp, 
                                sampler = sampler, 
                                settings = mcmc_settings)})
  
  #summer
  bayesianSetup_sm <- createBayesianSetup(prior = prior,
                                          likelihood = loglik_wrapper_sm,
                                          names = names(theta[index]),
                                          parallel = FALSE)
  
  system.time({trace_sm <- runMCMC(bayesianSetup = bayesianSetup_sm, 
                                   sampler = sampler, 
                                   settings = mcmc_settings)})
  
  #autumn
  bayesianSetup_au <- createBayesianSetup(prior = prior,
                                          likelihood = loglik_wrapper_au,
                                          names = names(theta[index]),
                                          parallel = FALSE)
  
  system.time({trace_au <- runMCMC(bayesianSetup = bayesianSetup_au, 
                                   sampler = sampler, 
                                   settings = mcmc_settings)})
  
  # winter
  bayesianSetup_wt <- createBayesianSetup(prior = prior,
                                          likelihood = loglik_wrapper_wt,
                                          names = names(theta[index]),
                                          parallel = FALSE)
  
  system.time({trace_wt <- runMCMC(bayesianSetup = bayesianSetup_wt, 
                                   sampler = sampler, 
                                   settings = mcmc_settings)})
  
}

# DIAGNOSTICS -----------------------------------------------

plot(trace_sp)
plot(trace_sm) 
plot(trace_au)
plot(trace_wt) 

# burn-in
nburn <- 10000
plot(trace_sp, parametersOnly = TRUE, start=nburn)
plot(trace_sm, parametersOnly = TRUE, start=nburn)
plot(trace_au, parametersOnly = TRUE, start=nburn)
plot(trace_wt, parametersOnly = TRUE, start=nburn)

# check convergence and correlations
gelmanDiagnostics(trace_sp, plot=TRUE, start=nburn)
gelmanDiagnostics(trace_sm, plot=TRUE, start=nburn)
gelmanDiagnostics(trace_au, plot=TRUE, start=nburn)
gelmanDiagnostics(trace_wt, plot=TRUE, start=nburn)
correlationPlot(getSample(trace_sp, parametersOnly = TRUE, coda=TRUE, start=nburn), density="smooth", thin=50)
marginalPlot(trace_sp, prior=T, singlePanel=T, start=nburn, nDrawsPrior = 1000)

# remove burn-in for trajsim simulation
tracefinal_sp <- getSample(trace_sp, parametersOnly = TRUE, coda=TRUE, start=nburn)
plot(tracefinal_sp)
effectiveSize(tracefinal_sp)

tracefinal_sm <- getSample(trace_sm, parametersOnly = TRUE, coda=TRUE, start=nburn)
plot(tracefinal_sm)
effectiveSize(tracefinal_sm)

tracefinal_au <- getSample(trace_au, parametersOnly = TRUE, coda=TRUE, start=nburn)
tracefinal_wt <- getSample(trace_wt, parametersOnly = TRUE, coda=TRUE, start=nburn)

# Posterior summary
summary(tracefinal_sp)
summary(tracefinal_sm)

# save the trace
#saveRDS(trace_sp, "trace_FOI_w_all_season_together_prop.rds")


# POSTPROCESSING AND RESULTS -----------------------------------

# Calculate simulated trajectory quantiles
trajsim <- maketrajsim(tracefinal_sp, theta, agepred_sp, model_sp, inits, 1000)
trajquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles) <- c("agemid", "low95", "median", "up95")

lambda_spquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"lambda_sp"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_spquantiles) <- c("agemid", "low95", "median", "up95")

trajsim_sm <- maketrajsim(tracefinal_sm, theta, agepred_sm, model_sm, inits, 1000)
trajquantiles_sm <- plyr::ddply(.data=trajsim_sm, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_sm) <- c("agemid", "low95", "median", "up95")

lambda_smquantiles <- plyr::ddply(.data=trajsim_sm, .variables="time", function(x) quantile(x[,"lambda_sm"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_smquantiles) <- c("agemid", "low95", "median", "up95")

trajsim_au <- maketrajsim(tracefinal_au, theta, agepred_au, model_au, inits, 1000)
trajquantiles_au <- plyr::ddply(.data=trajsim_au, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_au) <- c("agemid", "low95", "median", "up95")

lambda_auquantiles <- plyr::ddply(.data=trajsim_au, .variables="time", function(x) quantile(x[,"lambda_au"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_auquantiles) <- c("agemid", "low95", "median", "up95")

trajsim_wt <- maketrajsim(tracefinal_wt, theta, agepred_wt, model_wt, inits, 1000)
trajquantiles_wt <- plyr::ddply(.data=trajsim_wt, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles_wt) <- c("agemid", "low95", "median", "up95")

lambda_wtquantiles <- plyr::ddply(.data=trajsim_wt, .variables="time", function(x) quantile(x[,"lambda_wt"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambda_wtquantiles) <- c("agemid", "low95", "median", "up95")

wquantiles <- plyr::ddply(.data=trajsim_wt, .variables="time", function(x) quantile(x[,"mu_wt"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(wquantiles) <- c("agemid", "low95", "median", "up95")



# Plot fit and FOI
fit_sp <- ggplot() + theme_bw() + ggtitle("model fit in spring") +
  geom_point(data=spring.df, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=spring.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 


fit_sp

lambda_sp <- ggplot() + theme_bw() + ggtitle("FOI in spring") +
  geom_ribbon(data=lambda_spquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_spquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_sp

fit_sm <- ggplot() + theme_bw() + ggtitle("model fit in summer") +
  geom_point(data=summer.df, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=summer.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_sm, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_sm, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 


fit_sm

lambda_sm <- ggplot() + theme_bw() + ggtitle("FOI in summer") +
  geom_ribbon(data=lambda_smquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_smquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_sm

fit_au <- ggplot() + theme_bw() + ggtitle("model fit in autumn") +
  geom_point(data=autumn.df, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=autumn.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_au, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_au, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 


fit_au

lambda_au <- ggplot() + theme_bw() + ggtitle("FOI in autumn") +
  geom_ribbon(data=lambda_auquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_auquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_au


fit_wt <- ggplot() + theme_bw() + ggtitle("model fit in winter") +
  geom_point(data=winter.df, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=winter.df, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles_wt, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles_wt, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 


fit_wt


lambda_wt <- ggplot() + theme_bw() + ggtitle("FOI in winter") +
  geom_ribbon(data=lambda_wtquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambda_wtquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 
lambda_wt

w <- ggplot() + theme_bw() + ggtitle("Waning maternal immunity in winter") +
  geom_ribbon(data=wquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=wquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("Maternal immunity") 


w



