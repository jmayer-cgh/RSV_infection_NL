#
### Generalized Additive Models, Logistic regression for RSV serology.
### Authors: Stijn Andeweg, Jan van de Kassteele and Michiel van Boven
#

# Load packages
library(mgcv)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(grid)
library(writexl)

#set WD
folder <- "RSV_serology"
PATH <- file.path("C:/Users/julia/OneDrive - London School of Hygiene and Tropical Medicine/Documents/LSHTM/Project") #"/home/andewegs/1_RSV_scripts/"
setwd(PATH)
getwd()

# set figure save path
figure_folder <- '/home/andewegs/1_RSV_scripts/tensor_splines/figures/final_models/'

# set colors in figures
color_1 = '#0072B2' #Blue
color_2 = '#E69F00' #Orange
memory.size()
#
# Data ----
#

# Import data 
rsv.data_org <- read_csv(file = "data/infection_status_csv.txt")#"data/Riskfactors2_csv.txt")
season_border = "10-01" #MM-DD
# Some modifications
rsv.data <- rsv.data_org %>%
  # Get day of the year of birthday (= number between 1 and 365)
  mutate(
    Birth_doy = birthday %>% yday()
  ) %>%
  # Alter number of housefold members 1 == no siblings, >1 == having siblings
  mutate(
    household04_counts = case_when(
      age_days/365 >= 5 ~ (household04 + 1),  #because the child itself is also included in household size
      age_days/365 < 5 ~ household04)       #the household04 of children of age 5 should also include the child itself
  ) %>%
  # Make number of siblings 0-4 years binary factor variable
  mutate(
    Siblings04 = case_when(
      household04_counts <= 1 ~ 'False',
      household04_counts > 1 ~ 'True'
    ) %>% 
      factor()
  ) %>%
  # Make number of siblings 5-9 years binary factor variable
  mutate(
    Siblings59 = case_when(
      household59 <= 0 ~ 'False',
      household59 > 0 ~ 'True'
    ) %>% 
      factor()
  ) %>% 
  # Set nursery 0 1 to False True
  mutate(
    Nursery = case_when(
      visitnursery_child == 0 ~ 'False',
      visitnursery_child == 1 ~ 'True'
    ) %>% 
      factor()
  ) %>% 
  # Set variable for the different seasons.
  mutate(seasons = case_when(
    consultdate < paste("2006-", season_border, sep = "")  ~ "2005/2006", 
    (consultdate >= paste("2006-", season_border, sep = "") &  consultdate < "2010-01-01") ~ "2006/2007",
    (consultdate >= "2010-01-01" &  consultdate < paste("2016-", season_border, sep = "")) ~ "2015/2016",
    consultdate >= paste("2016-", season_border, sep = "")  ~ "2016/2017")
  ) %>% 
  # Variable for the two cohorts
  mutate(
    cohort = if_else(seasons == "2005/2006" | seasons == "2006/2007", "2006/2007", "2016/2017"))

# Only keep data with the variables of intrest n = 616
completeVec <- complete.cases(rsv.data[, c("Siblings04", "Nursery")])
rsv.data <- rsv.data[completeVec,]

#
# Analysis ----
#

# Model 1 with age and birth doy
model1 <- gam(
  formula = infection ~ 
    ti(age_days, bs = "ps", k = 25) + 
    ti(Birth_doy, bs = "cp", k = 11),# +
  #ti(age_days, Birth_doy, bs = c("ps", "cp")),
  knots = list(Birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model1)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = as.matrix(c(30,60,90,120,150,180,270,365,547,730,900,1095)),
  Birth_doy = seq(from = 1, to = 365, by = 1))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Figure 3

doy_age <- ggplot(
  data = rsv.doy_age.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr)) +
  geom_ribbon(
    alpha = 0.25) +
  geom_line() +
  ylab("Probability of prior infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "60" = "Age: 60 days",
                                                 "90" = "Age: 90 days",
                                                 "120" = "Age: 120 days",
                                                 "150" = "Age: 150 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age
ggsave(file = file.path(PATH, 'age_doy.svg'), plot = doy_age)

# Model1 for age time steps of 73 days
rsv.preddata <- expand.grid(
  age_days = as.matrix(c(73,146,219,292,365,
                         73+365,146+365,219+365,292+365,365+365,
                         73+365+365,146+365+365,219+365+365,292+365+365,365+365+365)),
  Birth_doy = seq(from = 1, to = 365, by = 5))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age.preddata_eqstep <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement figure
doy_age_73 <- ggplot(
  data = rsv.doy_age.preddata_eqstep,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr)) +
  geom_ribbon(
    alpha = 0.25) +
  geom_line() +
  ylab("Probability of prior infection") +
  xlab("Birth day of year") +
  facet_wrap(facets = ~ age_days, ncol = 5, labeller = labeller(age_days = 
                                                                  c("73" = "Age: 73 days",
                                                                    "146" = "Age: 146 days",
                                                                    "219" = "Age: 219 days",
                                                                    "292" = "Age: 292 days",
                                                                    "365" = "Age: 365 days",
                                                                    "438" = "Age: 438 days",
                                                                    "511" = "Age: 511 days",
                                                                    "584" = "Age: 584 days",
                                                                    "657" = "Age: 657 days",
                                                                    "730" = "Age: 730 days",
                                                                    "803" = "Age: 803 days",
                                                                    "876" = "Age: 876 days",
                                                                    "949" = "Age: 949 days",
                                                                    "1022" = "Age: 1022 days",
                                                                    "1095" = "Age: 1095 days"))) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_73
ggsave(file = file.path(PATH, 'age_doy_15grid.svg'), plot = doy_age_73)

## Model 2 with age, birth doy and Siblings04
model2 <- gam(
  formula = infection ~ 
    seasons + 
    ti(Birth_doy, bs = "cp", k = 11) + 
    ti(age_days, bs = "ps", k = 25),# +
  #ti(age_days, bs = "ps", k = 15, by = Siblings04),
  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model2)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Siblings04 = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model2,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement Figure
doy_age_sib <- ggplot(
  data = rsv.doy_age_sib.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr,
    group = Siblings04, fill = Siblings04, color = Siblings04
  )) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of prior infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "90" = "Age: 90 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  )  + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_sib
ggsave(file = file.path(PATH, 'age_doy_sib.svg'), plot = doy_age_sib)

# Additional model 2, with age, birth doy and Siblings59
model2.1 <- gam(
  formula = infection ~ 
    Siblings59 + 
    ti(Birth_doy, bs = "cp", k = 11) + 
    #ti(household04n, bs = "re", k = 10) +
    ti(age_days, bs = "ps", k = 25),# +
  #ti(age_days, bs = "ps", k = 15, by = siblings04),
  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model2.1)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Siblings59 = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model2.1,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib59.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement Figure
doy_age_sib59 <- ggplot(
  data = rsv.doy_age_sib59.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group= Siblings59, col= Siblings59, fill = Siblings59)) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of prior infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "90" = "Age: 90 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  )  + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  )

doy_age_sib59
ggsave(file = file.path(PATH, 'age_doy_sib59.svg'), plot = doy_age_sib59)

# Model3 with age, birth doy and nursery
model3 <- gam(
  formula = infection ~ 
    Nursery + 
    ti(age_days, bs = "ps", k = 25) + 
    ti(Birth_doy, bs = "cp", k = 11),# +
  #ti(contact04, bs = "re", k = 20) +
  #ti(age_days, bs = "ps", k = 25, by = visitnursery_child, m = 1) + 
  #ti(birth_doy, bs = "cp", k = 10, by = visitnursery_child, m = 1),
  # ti(contact04, bs = "re", k = 20, by = siblings04, m = 1), 
  
  knots = list(birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model3)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,90,180,270,365,547,730,900,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Nursery = factor(c('False', 'True')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model3,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_nur.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Supplement Figure
doy_age_nur <- ggplot(
  data = rsv.doy_age_nur.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group= Nursery, col= Nursery, fill = Nursery)) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of prior infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days, labeller = labeller(age_days = 
                                               c("30" = "Age: 30 days",
                                                 "90" = "Age: 90 days",
                                                 "180" = "Age: 180 days",
                                                 "270" = "Age: 270 days",
                                                 "365" = "Age: 365 days",
                                                 "547" = "Age: 547 days",
                                                 "730" = "Age: 730 days",
                                                 "900" = "Age: 900 days",
                                                 "1095" = "Age: 1095 days"))
  ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "grey65", size = 0.2),
        strip.background = element_rect(
          color="black", fill="white", linetype = "blank"),
        strip.text = element_text(size=8)
  ) + guides(fill=guide_legend(title="Day-care"), col=guide_legend(title="Day-care"))

doy_age_nur   
ggsave(file = file.path(PATH, 'age_doy_nur.svg'), plot = doy_age_nur)

# Model4 with age, birth doy, Siblings04 and nursery
model4 <- gam(
  formula = infection ~ 
    Nursery + 
    Siblings04 +
    ti(age_days, bs = "ps", k = 25) + 
    ti(Birth_doy, bs = "cp", k = 11),# +
  #ti(contact04, bs = "re", k = 10) +
  #ti(age_days, bs = "ps", k = 25, by = siblings04, m = 1) + 
  #ti(birth_doy, bs = "cp", k = 10, by = siblings04, m = 1) +
  # ti(contact04, bs = "re", k = 10, by = siblings04, m = 1) + 
  #ti(age_days, bs = "ps", k = 25, by = visitnursery_child, m = 1) + 
  #ti(birth_doy, bs = "cp", k = 10, by = visitnursery_child, m = 1),# +
  # ti(contact04, bs = "re", k = 10, by = visitnursery_child, m = 1),# 
  
  knots = list(Birth_doy = c(1, 365)),
  family = binomial,
  method = "REML",  #restricted maximum likelihood 
  data = rsv.data)
summary(model4)

# Generate some fake data for the model prediction
rsv.preddata <- expand.grid(
  age_days = c(30,180,365,547,730,1095),
  Birth_doy = seq(from = 1, to = 365, by = 5),
  Siblings04 = factor(c('True', 'False')),
  Nursery = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model4,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib_nur.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

#
# Plot results ----
#

# Plot Figure 4
doy_age_sib_nur <- ggplot(
  data = rsv.doy_age_sib_nur.preddata,
  mapping = aes(
    x = Birth_doy, y = p_inf, ymin = p_inf_lwr, ymax = p_inf_upr, 
    group=Nursery , col= Nursery, fill = Nursery)) +
  geom_ribbon(
    alpha = 0.25, colour = NA) +
  geom_line() +
  scale_fill_manual(values=c(color_1, color_2)) + 
  scale_color_manual(values=c(color_1, color_2)) + 
  ylab("Probability of prior infection") +
  xlab("Birth day of year") +
  facet_wrap(
    facets = ~ age_days + Siblings04, labeller = labeller(age_days = 
                                                            c("30" = "Age: 30 days",
                                                              "90" = "Age: 90 days",
                                                              "180" = "Age: 180 days",
                                                              "270" = "Age: 270 days",
                                                              "365" = "Age: 365 days",
                                                              "547" = "Age: 547 days",
                                                              "730" = "Age: 730 days",
                                                              "900" = "Age: 900 days",
                                                              "1095" = "Age: 1095 days"),
                                                          Siblings04 = c("False" = "Siblings04: False",
                                                                         "True" = "Siblings04: True"), .multi_line = FALSE)
  ) + theme_bw() + 
  theme(panel.grid.minor = element_blank(), #remove grid lines between ticks
        panel.grid.major = element_line(colour = "grey65", size = 0.2), #set size and color of grid lines
        strip.background = element_rect(color="black", fill="white", linetype = "blank"), #set the background
        legend.position = "bottom", #remove legend
        strip.text = element_text(size=7.4) #set text size for text above the plots
  ) + guides(fill=guide_legend(title="Day-care"), col=guide_legend(title="Day-care"))

doy_age_sib_nur
ggsave(file = file.path(PATH, 'age_doy_sib_nur.svg'), plot = doy_age_sib_nur)

# Table 1
rsv.preddata <- expand.grid(
  age_days = c(30,365,730),
  Birth_doy = c(1,182), # 1 Jan, 1 July
  Siblings04 = factor(c('True', 'False')),
  Nursery = factor(c('True', 'False')))

# Make predictions for each combination using the GAM model
# Do this in a temporary tibble, because two columns are produced: fit and se.fit
tmp <- predict(
  object = model4,
  newdata = rsv.preddata,
  type = "link",
  se.fit = TRUE) %>%
  as_tibble

# Bind tmp to rsv.preddata and calculate the p_inf including 95% lower and upper bound
rsv.doy_age_sib_nur.preddata <- bind_cols(rsv.preddata, tmp) %>%
  mutate(
    fit_lwr = fit + qnorm(0.025)*se.fit,   
    fit_upr = fit + qnorm(0.975)*se.fit,
    p_inf = fit %>% plogis,                #create vector of probabilities from log odds
    p_inf_lwr = fit_lwr %>% plogis,
    p_inf_upr = fit_upr %>% plogis)

table_model4 <- rsv.doy_age_sib_nur.preddata %>% 
  mutate(
    p_inf = p_inf,# %>% round(digits = 3),
    p_inf_lwr = p_inf_lwr,# %>% round(digits = 3),
    p_inf_upr = p_inf_upr,# %>% round(digits = 3),
    n = rsv.data %>% nrow
  ) %>% 
  select(`Age (days)` = age_days, 
         DOY = Birth_doy,
         Siblings04,
         `Day-care` = Nursery,
         n,
         `Probability of prior infection` = p_inf,
         `Lower` = p_inf_lwr,
         `Upper` = p_inf_upr) %>% 
  arrange(`Age (days)`, DOY, Siblings04) %>% 
  write_xlsx(path = "tabel_model4.xlsx")

# table 1
table_complete <- table_model4 %>% 
  left_join(table_model1, by = c("Age (days)", "DOY")) %>% 
  left_join(table_model2, by = c("Age (days)", "DOY", "Siblings04")) %>% 
  left_join(table_model3, by = c("Age (days)", "DOY", "Day-care")) %>% 
  select(`Age (days)`, DOY, `Probability of prior infection m1`, `Lower m1`, `Upper m1`,
         `Siblings04`, `Probability of prior infection m2`, `Lower m2`, `Upper m2`,
         `Day-care`, `Probability of prior infection m3`, `Lower m3`, `Upper m3`,
         `Probability of prior infection m4`, `Lower m4`, `Upper m4`) %>% 
  write_xlsx(path = "voorbeeld_tabel_versie2.xlsx")




rm(age_groups)

#The examples/numbers given in the manuscript
# Abstract summer vs winter 
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 1,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 182,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# Abstract siblings04 
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 181 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 181 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# Abstract nursery
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 181 & 
                           rsv.doy_age_nur.preddata$Nursery == 'True',]
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 181 & 
                           rsv.doy_age_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# Results July vs Januari
#half a year
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 182 & rsv.doy_age.preddata$Birth_doy == 15,c('p_inf', 'p_inf_lwr', 'p_inf_upr')] #jan. 15
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 182 & rsv.doy_age.preddata$Birth_doy == 196,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]# juli 15
# one year
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 15,c('p_inf', 'p_inf_lwr', 'p_inf_upr')] #jan. 15
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 365 & rsv.doy_age.preddata$Birth_doy == 196,c('p_inf', 'p_inf_lwr', 'p_inf_upr')]# juli 15

# Results siblings0-4
# age half a year
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 180 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 180 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# age 1 year
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_sib.preddata[rsv.doy_age_sib.preddata$age_days == 365 & 
                           rsv.doy_age_sib.preddata$Birth_doy == 196 & 
                           rsv.doy_age_sib.preddata$Siblings04 == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

# Results day-care
# age half a year
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 180 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 180 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# age a year
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
rsv.doy_age_nur.preddata[rsv.doy_age_nur.preddata$age_days == 365 & 
                           rsv.doy_age_nur.preddata$Birth_doy == 196 & 
                           rsv.doy_age_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

#HR vs LR
# at half a year
rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 180 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'True' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 180 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'False' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

# at a year
rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 365 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'True' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'True',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

rsv.doy_age_sib_nur.preddata[rsv.doy_age_sib_nur.preddata$age_days == 365 & 
                               rsv.doy_age_sib_nur.preddata$Birth_doy == 196 & 
                               rsv.doy_age_sib_nur.preddata$Siblings04 == 'False' &
                               rsv.doy_age_sib_nur.preddata$Nursery == 'False',c('p_inf', 'p_inf_lwr', 'p_inf_upr')]

#Discussion
#month one sumer
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 30 & rsv.doy_age.preddata$Birth_doy == 196, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
#month 24 sumer
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 395 & rsv.doy_age.preddata$Birth_doy == 196, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# month one winter
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 30 & rsv.doy_age.preddata$Birth_doy == 15, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]
# month 24 winter
rsv.doy_age.preddata[rsv.doy_age.preddata$age_days == 395 & rsv.doy_age.preddata$Birth_doy == 15, c('p_inf', 'p_inf_lwr', 'p_inf_upr')]



#
### Model fit measurements ----
#
AIC.data <- data.frame(
  model = c('model1', 'model2', 'model3', 'model4'),
  score = c(model1$aic, model2$aic, model3$aic, model4$aic)
)

BIC.data <- data.frame(
  model = c('model1', 'model2','model3', 'model4'),
  score = c(BIC(model1), BIC(model2), BIC(model3),BIC(model4))
)

r.sq.data <- data.frame(
  model = c('model1', 'model2', 'model3', 'model4'),
  score = c(summary(model1)$r.sq, summary(model2)$r.sq, summary(model3)$r.sq, summary(model4)$r.sq)
)
AIC.data$method <- "AIC"
BIC.data$method <- "BIC"
r.sq.data$method <- "r.sq"
model_scores <- rbind(AIC.data, BIC.data, r.sq.data)