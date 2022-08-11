################################################################
# RSV seroconversion MSc project
# Data plots
# Author: Julia Mayer
# Last updated: 11.08.2022
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

data.df <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/master/data/infection_status_csv.txt",
                 sep=",")
# use if offline
#library(readr)
#data <- read_csv("LSHTM/Project/Data/infection_status.csv")

# Group age into intervals 
# bi-monthly for 0-2 years and 6-monthly for 2-5 years
data.df$agegrp <- cut(data.df$age_days,
                   breaks=c(seq(0,730, by=30.25*2),
                            seq(909,2000, by=30.25*6)), 
                   include.lowest = T, right=F)
# Divide by season of birth
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)

data.df <- data.df %>%
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
data.day_care <- subset(data.df, !is.na(visitnursery_child))
data.day_care <- data.day_care%>% dplyr::group_by(agegrp, season_birth, visitnursery_child) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection))

data_no_season.df <- data.df %>% dplyr::group_by(agegrp) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection))

data.df <- data.df %>% dplyr::group_by(agegrp, season_birth) %>% 
  dplyr::summarise(agemid=round(median(age_days)), # Age midpoint in age group
                   N=n(), # Total N in age group
                   nconv=sum(infection)) # n seroconverted in age group



# Calculate seroprevalence and binomial confidence intervals
data.df[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data.df$nconv, data.df$N, method="exact")[,c("mean","lower","upper")]
data_no_season.df[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data_no_season$nconv, data_no_season$N, method="exact")[,c("mean","lower","upper")]
data.day_care[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom.confint(data.df$nconv, data.df$N, method="exact")[,c("mean","lower","upper")]

# Plot the whole data frame (difficult to read)
ggplot(data.day_care) +
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

ggplot(data.df) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = season_birth)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = season_birth)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + labs (colour = "Day-care attendance") 



#Plot by season
spring.df2 <- subset(data.day_care, season_birth == 'Spring')
ggplot(spring.df2) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in spring", colour = "Day-care")

spring.df3 <- subset(data.df, season_birth == 'Spring')
ggplot(spring.df3) +
  geom_point(aes(x=agemid, y=seroprev_mean)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in spring")


summer.df2 <- subset(data.day_care, season_birth == 'Summer')
ggplot(summer.df2) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in summer", colour = "Day-care attendance")

summer.df3 <- subset(data.df, season_birth == 'Summer')
ggplot(summer.df3) +
  geom_point(aes(x=agemid, y=seroprev_mean)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in summer")

autumn.df2 <- subset(data.day_care, season_birth == 'Autumn')
ggplot(autumn.df2) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in autumn", colour = "Day-care attendance")

autumn.df3 <- subset(data.df, season_birth == 'Autumn')
ggplot(autumn.df3) +
  geom_point(aes(x=agemid, y=seroprev_mean)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in autumn")


winter.df2 <- subset(data.day_care, season_birth == 'Winter')
ggplot(winter.df2) +
  geom_point(aes(x=agemid, y=seroprev_mean, colour = factor(visitnursery_child))) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95, colour = factor(visitnursery_child))) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in winter", colour = "Day-care attendance")

winter.df3 <- subset(data.df, season_birth == 'Winter')
ggplot(winter.df3) +
  geom_point(aes(x=agemid, y=seroprev_mean)) +
  geom_errorbar(aes(x=agemid, ymin=seroprev_low95, ymax=seroprev_up95)) +
  ylab("Proportion seroconverted") + xlab("age (days)") + 
  labs(title ="Proportion seroconverted born in winter")

