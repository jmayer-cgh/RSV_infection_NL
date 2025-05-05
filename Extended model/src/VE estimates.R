# Vaccine efficacy, taken from Ayaka

#### libraries ####
library(tidyverse)
library(epitools)
library(BayesianTools)
library(knitr)

#### Data preparation ####
#### Phase 3 trial of Pfizer's RSVpreF (Munjai et al. 2024/2. RSVVW'24)
#### VE vs RSV-positive severe MA-LRTI
N_placebo = 3563
N_intervention = 3585
tibble(Tmin = c(0,30,60,90,120,150), # Lower bound of time period (days after birth)
       Tmax = c(30,60,90,120,150,180), # Upper bound of time period (days after birth)
       Intervention = diff(c(0,1,4,6,13,18,21)), # Number of RSV-positive severe MA-LRT in intervention group at a given time point (OG paper gives cumulative cases so we need to back-calculate this)
       Placebo = diff(c(0,10,28,34,49,61,70))) %>% # Number of RSV-positive severe MA-LRT in placebo group at a given time point
  rowwise() %>%
  mutate(Tmid = (Tmin+Tmax)/2, # Mid-point of time interval
         # VE against RSV-positive severe MA-LRTI
         VEmid = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,1],
         VEhi = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,2],
         VElo = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,3],
         group = 'severe') -> dat_mat_sev


##### VE vs RSV-positive less severe MA-LRTI
tibble(Tmin = c(0,30,60,90,120,150),
       Tmax = c(30,60,90,120,150,180),
       Intervention = diff(c(0,2,14,25,40,55,67)),
       Placebo = diff(c(0,15,38,59,88,110,132)) ) %>%
  rowwise() %>%
  mutate(Tmid = (Tmin+Tmax)/2,
         VEmid = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,1],
         VEhi = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,2],
         VElo = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,3],
         group = 'LRTI') -> dat_mat_LRTI


#### get corresponding VE with Erlang decay ####
erlang.decay = function(VE0, T, k=2, t=1:150) {
  res = 0
  for(n in 0:(k-1)){
    res = res + 1/factorial(n)*exp(-T*t)*(T*t)^n } # T is the rate
  return( VE0 * res)
}


#### define LL ####
LL = function(x){
  VE0_s = x[1] # severe
  VE0_l = x[2] # less-severe
  T = x[3]
  dat_mat_sev %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_s, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_s
  dat_mat_LRTI %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_l, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_l
  return(sum(datatmp_s$LL) + sum(datatmp_l$LL))
}
LL_freq = function(x) -LL(x)

#### run fitting - estimating VE0_s, VE0_l, and T####
# Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(0.001,0.001,0.001), upper = c(.99,.99,.99)) # Limits for priors for VE0_s, VE0_l, T are 0.001 and 0.99
iter = 100000
settings = list(iterations = iter, message = FALSE, nrChains = 2)
out_ba <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

#### get model predictions ####
out_ba[[1]]$chain[seq(1,iter,by=10),1:3] %>% # 10000 samples from chain 1?
  as_tibble() %>%
  pivot_longer(-'3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='3') %>% # T is the waning rate
  mutate(group = case_match(group, '1'~'severe', '2'~'LRTI')) -> resp

# Outputs by severity of disease
# This gives the distribution of VE over time based on the initial values and T fitted above
model_ve_runs <- resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1]))

# Identify which numbers were generated from the same initial values
model_ve_runs_t <- model_ve_runs %>%
  mutate(iter = rep(1:ceiling(nrow(model_ve_runs) / 730), each = 730)[1:nrow(model_ve_runs)])

# Severe
model_ve_severe_sum <- model_ve_runs %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "severe")

# Mild
model_ve_less_sum <- model_ve_runs %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "LRTI")

# Combine the two into one dataframe
model_outputs_sum <- rbind(model_ve_severe_sum, model_ve_less_sum)


