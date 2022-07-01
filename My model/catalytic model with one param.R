################################################################
# RSV seroconversion MSc project
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

data <- data %>% dplyr::group_by(agegrp) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection)) # n seroconverted in age group

# Calculate seroprevalence and binomial confidence intervals
data[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data$nconv, data$N, method="exact")[,c("mean","lower","upper")]

# Plot
ggplot(data) +
  geom_point(aes(x=agemid, y=seroprev_mean), color="black") +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95), colour="black") +
  ylab("Proportion seroconverted") + xlab("age (days)")


# MODEL EQUATION  ------------------------------------------------------------

# Notes:
# - This model does not contain a maternal immunity compartment, because of how "infection" is defined in the data, see Andeweg et al
# - The states are integrated on a log-scale to avoid negative states, which is why we log-transform them when in the model input and exponentiate them inside the model

model <- function(theta, age, inits) {
  
  catalytic <- function(age, state, param) {
    
    # FOI / seroconversion rate
    lambda = param[["A"]] + age * 0 #constant FOI for now, change 0 to something else later
    
    # states 
    S = exp(state[1]) # susceptible
    Z = exp(state[2]) # seroconverted after infection
    
    # changes in states
    dS = - lambda*S
    dZ = + lambda*S
    
    return(list(c(dS/S, dZ/Z), lambda=lambda))
    
    
  }
  
  traj <- data.frame(ode(y=c(S=log(inits[["S"]]),
                             Z=log(inits[["Z"]])),
                         times=age, 
                         func=catalytic, 
                         parms=theta, #what is theta? the set of parameters
                         method="lsoda",
                         verbose=F))
  
  traj$conv <- exp(traj$Z) # cumulative seroconversion (=observed state)
  traj$inc <- c(inits[["Z"]], diff(exp(traj$Z))) # incident seroconversion
  
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

# A = mean FOI (proportion infected per day)
theta <- c(A=0.05) # these are just random values, to be fitted


# INITS ---------------------------------------------------------

inits <- c(S=1-1e-12, Z=1e-12) # initial conditions for the states (as proportions)
# --> since we integrate on a log-scale, the initial conditions cannot be 0 (not defined on a log-scale)

# SIMULATION TIME  ---------------------------------------------------------

#agepred <- seq(min(data$agemid), max(data$agemid), by=1) # if you want a higher resolution of the solution, but is computationally more expensive
agepred <- data$agemid

# TEST MODEL  ---------------------------------------------------------

test <- model(theta, agepred, inits)
ggplot(test) + geom_line(aes(x=time, y=conv))
ggplot(test) + geom_line(aes(x=time, y=lambda)) #should be constant

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
estpars <- c("A") # parameters to estimate, can be modified
index <- which(names(theta) %in% estpars) # index of estimated params


# Priors
lower = c(A=0)
upper = c(A=0.1)

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

plot(trace) #put plot in full screen to work

# burn-in
nburn <- 10000
plot(trace, parametersOnly = TRUE, start =nburn)

# check convergence and correlations
gelmanDiagnostics(trace, plot=TRUE, start=nburn)
#correlationPlot(getSample(trace, parametersOnly = TRUE, coda=TRUE, start=nburn), density="smooth", thin=50)
#marginalPlot(trace, prior=T, singlePanel=T, start=nburn, nDrawsPrior = 1000)

# remove burn-in for trajsim simulation
tracefinal <- getSample(trace, parametersOnly = TRUE, coda=TRUE, start=nburn)
plot(tracefinal) #put in full screen to work
effectiveSize(tracefinal)

# Posterior summary
summary(tracefinal)

# save the trace
saveRDS(trace, "trace_A.rds")

# POSTPROCESSING AND RESULTS -----------------------------------

# Calculate simulated trajectory quantiles
trajsim <- maketrajsim(tracefinal, theta, agepred, model, inits, 1000)
trajquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"conv"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(trajquantiles) <- c("agemid", "low95", "median", "up95")

lambdaquantiles <- plyr::ddply(.data=trajsim, .variables="time", function(x) quantile(x[,"lambda"], prob = c(0.025, 0.5, 0.975), na.rm=T)) 
colnames(lambdaquantiles) <- c("agemid", "low95", "median", "up95")


# Plot fit and FOI
fit <- ggplot() + theme_bw() + ggtitle("model fit") +
  geom_point(data=data, aes(x=agemid, y=seroprev_mean)) +
  geom_linerange(data=data, aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  geom_ribbon(data=trajquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=trajquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("proportion seroconverted") 

fit #full screen

lambda <- ggplot() + theme_bw() + ggtitle("FOI") +
  geom_ribbon(data=lambdaquantiles, aes(x=agemid, ymin=low95, ymax=up95), fill="red", alpha=0.3) +
  geom_line(data=lambdaquantiles, aes(x=agemid, y=median), color="red") +
  xlab("age (days)") + ylab("FOI") 


lambda #full screen





