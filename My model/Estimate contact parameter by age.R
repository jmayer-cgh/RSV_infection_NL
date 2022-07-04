################################################################
# RSV seroconversion MSc project
# Contact patterns by age
# Author: Julia Mayer
# Last updated: 04.07.2022
################################################################

rm(list=ls())
set.seed(42)

library(tidyverse)
library(plyr)
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

season_border = "10-01" #MM-DD
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)

# Define categories like in the original model
data <- data %>%
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
  )
data <- data %>%
  mutate(
    Birth_mo = birthday %>% month()
  )%>%
  mutate(
    seasons = case_when(
      consultdate < paste("2006-", season_border, sep = "")  ~ "2005/2006", 
      (consultdate >= paste("2006-", season_border, sep = "") &  consultdate < "2010-01-01") ~ "2006/2007",
      (consultdate >= "2010-01-01" &  consultdate < paste("2016-", season_border, sep = "")) ~ "2015/2016",
      consultdate >= paste("2016-", season_border, sep = "")  ~ "2016/2017")
  ) %>% 
  # Variable for the two cohorts
  mutate(
    cohort = if_else(seasons == "2005/2006" | seasons == "2006/2007", "2006/2007", "2016/2017")
  )%>%
  mutate(
    season_birth = case_when(
      Birth_mo %in% spring ~ "Spring",
      Birth_mo %in% summer ~ "Summer",
      Birth_mo %in% autumn ~ "Autumn",
      Birth_mo %in% winter ~ "Winter")
  ) %>%
  mutate(
    nursery_house = case_when(
      visitnursery_house <= 0 ~ 'False',
      visitnursery_house > 0 ~ 'True'
    ) %>% 
      factor()
  )%>%
  mutate(
    pre_term = case_when(
      pregnancytime < 37.0 ~ 'True',
      pregnancytime >= 37.0 ~ 'False'
    ) %>%
      factor()
  )

# Group number of contacts
data <- data %>%
  mutate(
    total_cont = case_when(
      contacttotal == 0 ~ '0',
      contacttotal >= 1 & contacttotal <= 2 ~ '1-2',
      contacttotal >= 3 & contacttotal <= 5 ~ '3-5',
      contacttotal >= 6 & contacttotal <= 10 ~ '6-10',
      contacttotal >= 11 & contacttotal <= 20 ~ '11-20',
      contacttotal >= 21 & contacttotal <= 50 ~ '21-50',
      contacttotal > 50  ~ '50+'
    )
  )%>%
  mutate(
    total_cont59 = case_when(
      contact59 == 0 ~ '0',
      contact59 >= 1 & contact59 <= 2 ~ '1-2',
      contact59 >= 3 & contact59 <= 5 ~ '3-5',
      contact59 >= 6 & contact59 <= 10 ~ '6-10',
      contact59 >= 11 & contact59 <= 20 ~ '11-20',
      contact59 >= 21 & contact59 <= 50 ~ '21-50'
    )
  )%>%
  mutate(
    total_cont04 = case_when(
      contact04 == 0 ~ '0',
      contact04 >= 1 & contact04 <= 2 ~ '1-2',
      contact04 >= 3 & contact04 <= 5 ~ '3-5',
      contact04 >= 6 & contact04 <= 10 ~ '6-10',
      contact04 >= 11 & contact04 <= 20 ~ '11-20',
      contact04 >= 21 & contact04 <= 40 ~ '21-40'
    )
  )

# Infection by age
tab <- table(data$agegrp, data$infection)
ptab <- prop.table(tab, margin=1)
dframe <- data.frame(values=rownames(tab), infected=ptab[,2])
p<-ggplot(data=dframe, aes(x=values, y=infected, fill = values)) +
  geom_bar(stat="identity")+
  theme_minimal()+
  labs(title="Proportion infected by age",
       x ="Age in days", y = "Proportion infected") + scale_fill_discrete(name="Age in days")
p


# Contacts by age
data %>% filter(complete.cases(contacttotal)) %>% 
  group_by(agegrp) %>% 
  ggplot(aes(x=agegrp,y=contacttotal,fill=factor(infection)))+geom_boxplot()+
  labs(title="Total number of contacts by age and infection status",
       x ="Age in days", y = "Total contacts", fill = "Infection status") 
