# mAB efficacy, taken from Ayaka

set.seed(1234)

### dependencies
if (!require(pacman)){ #load packages
  install.packages("pacman")
}

library(pacman)
pacman::p_load(char = c("tidyverse", "here","epitools",
                        "BayesianTools", "knitr", "ggplot2",  "MASS",
                        "dplyr", "conflicted", "ggtext", "MCMCpack",
                        "readxl", "fitR", "coda", "lattice", "mosaic"
))

conflicts_prefer(dplyr::filter)
conflicts_prefer(mosaic::mean)
conflicts_prefer(mosaic::sum)
conflicts_prefer(dplyr::select)

#### Data preparation ####
# Efficacy results from phase 2b and 3 trial of nirsevimab 
# Simões EAF, et al. Efficacy of nirsevimab against respiratory syncytial virus lower respiratory tract infections in preterm and term infants, and pharmacokinetic extrapolation to infants with congenital heart disease and chronic lung disease: a pooled analysis of randomised controlled trials. Lancet Child Adolesc Health. 2023 Mar;7(3):180-189. doi: 10.1016/S2352-4642(22)00321-2. Epub 2023 Jan 9. PMID: 36634694; PMCID: PMC9940918.)
# Figure 1


# calculate incremental efficacy for each time interval from trial data

#Efficacy vs first medically attended RSV LRTI
N_placebo <-  786 # number at risk at time 0 in placebo arm
N_intervention <-  1564 # number at risk at time 0 in intervention arm

# intervention arm
x <-   c(1564, 1553, 1546, 1538, 1527, 1519) 
NumAtRisk_int <- (head(x, -1) + tail(x, -1)) / 2 # number at risk as average of numbers at risk at beginning and end of that interval
NumCen_int <-  c(0, 8, 11, 17, 21, 1545) # number censored
big_int <- c(1564, 1553, 1546, 1538, 1527, 0) #number at risk at beginning of that interval (except the last number)
# number censored at end of study is considered to include number at risk

# placebo arm
y <- c(786, 772, 756, 737, 729, 724)
NumAtRisk_pla <- (head(y, -1) + tail(y, -1)) / 2
NumCen_pla  <-  c(0, 6, 7, 8, 9, 735)
big_pla <-  c(786, 772, 756, 737, 729, 0)

# calculate incremental efficacy during each time interval
dat_mab_ma <- tibble(Tmin = c(0,30,60,90,120),
                     Tmax = c(30,60,90,120,150),
                     Tmid = (Tmin+Tmax)/2,
                     inc_cases_int = diff(N_intervention - big_int - NumCen_int), # incremental cases during each interval
                     inc_cases_pla = diff(N_placebo - big_pla - NumCen_pla),
                     num_at_risk_int = NumAtRisk_int,
                     num_at_risk_pla = NumAtRisk_pla) %>%
  rowwise() %>%
  # incremental efficacy during each interval
  mutate(VEmid = 1-riskratio.wald(rbind(c(num_at_risk_pla,inc_cases_pla), c(num_at_risk_int,inc_cases_int)))$measure[2,1], 
         VEhi = 1-riskratio.wald(rbind(c(num_at_risk_pla,inc_cases_pla), c(num_at_risk_int,inc_cases_int)))$measure[2,2],
         VElo = 1-riskratio.wald(rbind(c(num_at_risk_pla,inc_cases_pla), c(num_at_risk_int,inc_cases_int)))$measure[2,3],
         group = 'medically-attended') 

# Efficacy vs first hospital admission for medically attended RSV LRTI
# intervention arm
x_hosp <-   c(1564, 1554, 1547, 1540, 1535, 1529)
NumAtRisk_int_hosp <- (head(x_hosp, -1) + tail(x_hosp, -1)) / 2
NumCen_int_hosp <-  c(0, 8, 11, 17, 21, 1555)
big_int_hosp <-c(1564, 1554, 1547, 1540, 1535, 0)

# placebo arm
y_hosp <-  c(786, 778, 769, 761, 757, 753)
NumAtRisk_pla_hosp <- (head(y_hosp, -1) + tail(y_hosp, -1)) / 2
NumCen_pla_hosp <-  c(0, 6, 7, 8, 10, 765)
big_pla_hosp <-  c(786, 778, 769, 761, 757, 0)

dat_mab_hosp <- tibble(Tmin = c(0,30,60,90,120),
                       Tmax = c(30,60,90,120,150),
                       Tmid = (Tmin+Tmax)/2,
                       inc_cases_int = diff(N_intervention - big_int_hosp - NumCen_int_hosp),
                       inc_cases_pla = diff(N_placebo - big_pla_hosp - NumCen_pla_hosp),
                       num_at_risk_int = NumAtRisk_int_hosp[1:5],
                       num_at_risk_pla = NumAtRisk_pla_hosp[1:5]) %>%
  rowwise() %>%
  mutate(VEmid = 1-riskratio.wald(rbind(c(num_at_risk_pla,inc_cases_pla), c(num_at_risk_int,inc_cases_int)))$measure[2,1],
         VEhi = 1-riskratio.wald(rbind(c(num_at_risk_pla,inc_cases_pla), c(num_at_risk_int,inc_cases_int)))$measure[2,2],
         VElo = 1-riskratio.wald(rbind(c(num_at_risk_pla,inc_cases_pla), c(num_at_risk_int,inc_cases_int)))$measure[2,3],
         group = 'hospitalisation')

#### Waning function following Erlang distribution ####
erlang.decay = function(VE0, T, k=2, t=1:150) {
  res = 0
  for(n in 0:(k-1)){ #counter cdf of Erlang-k distribution
    res = res + 1/factorial(n)*exp(-T*t)*(T*t)^n }　
  return( VE0 * res)
}


#### define LL ####
LL = function(x){
  VE0_m = x[1]
  VE0_h = x[2]
  T = x[3]
  dat_mab_ma %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_m, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(inc_cases_int,N_intervention,(1-VE_est)*inc_cases_pla/N_placebo))) -> datatmp_m
  dat_mab_hosp %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_h, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(inc_cases_int,N_intervention,(1-VE_est)*inc_cases_pla/N_placebo))) -> datatmp_h
  return(sum(datatmp_m$LL) + sum(datatmp_h$LL))
}
LL_freq = function(x) -LL(x)

#### run fitting ####
# frequentist check
out_freq <- optim(c(.5,.5,1/365), LL_freq, lower = c(0.001,0.001,0.001), upper = c(1,1,1), method= 'L-BFGS-B')

# Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(0.001,0.001,0.001), upper = c(.99,.99,.99))

settings <- list(
  sampler = "DEzs",      
  #  iterations = 10000,    
  iterations = 100000,  
  burnin =500,            
  message = FALSE, 
  nrChains = 1,  
  thin =2
) 

out_bay <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

#### convergence, ESS check ####
getSample(out_bay, start = settings$burnin, coda = T, parametersOnly = FALSE) # include log posterior and log likelihood and log prior
chains_mcmc <- getSample(out_bay, start = settings$burnin , coda = T) 

# convergence check
plot(chains_mcmc)
gelman.diag(chains_mcmc, autoburnin = FALSE)
effectiveSize(chains_mcmc) 

#  make dataframe for MCMC samples
posterior_df <- bind_rows(lapply(seq_along(chains_mcmc), function(i) {
  as.data.frame(chains_mcmc[[i]]) %>% mutate(chain = i)
}))

# randomise posterior samples
posterior_df <- posterior_df %>% slice_sample(prop = 1)  # shuffle

# posterior samples
resp_mab <- posterior_df %>%
  select('par 1', 'par 2', 'par 3') %>% 
  pivot_longer(-'par 3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='par 3') %>%
  mutate(group = case_match(group, 'par 1'~'medically-attended', 'par 2'~'hospitalisation')) 

# Outputs by severity of disease
# This gives the distribution of IE over time based on the initial values and T fitted above
model_mab_runs <- resp_mab %>%
  rowwise() %>%
  mutate(IE_t = list(ie = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(IE_t) %>%
  mutate(t = rep(1:365, times=dim(resp_mab)[1]))

# Identify which numbers were generated from the same initial values
model_mab_runs_t <- model_mab_runs %>%
  mutate(iter = rep(1:ceiling(nrow(model_mab_runs) / 365), each = 365)[1:nrow(model_mab_runs)])

# Hospitalisations
model_ie_hosp_sum <- model_mab_runs %>%
  group_by(group,t) %>%
  summarise(IE_t_mid = median(IE_t),
            IE_t_lo = quantile(IE_t, probs = .025),
            IE_t_hi = quantile(IE_t, probs = .975))%>%
  filter(group == "hospitalisation")

# Mild
model_ie_med_sum <- model_mab_runs %>%
  group_by(group,t) %>%
  summarise(IE_t_mid = median(IE_t),
            IE_t_lo = quantile(IE_t, probs = .025),
            IE_t_hi = quantile(IE_t, probs = .975))%>%
  filter(group == "medically-attended")

# Combine the two into one dataframe
model_outputs_sum <- rbind(model_ie_hosp_sum, model_ie_med_sum)
