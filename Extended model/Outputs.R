# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)

# Path to files
path_output <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/CSV files/2 M odin/monty/"  #Where model outputs are stored
path_pop <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/Helpful papers/German population estimates/"
path_code <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/Extended model/src/"

# Read in files
# Cases in Germany
severe_illness_de <- read.csv(paste0(path_output, "Severe cases DE.csv"))
mild_illness_de <- read.csv(paste0(path_output, "Mild cases DE.csv"))

# Run model 100 times
for (i in 1:100){
  source(paste0(path_code, "VE estimates.R"))
  
  VE_distribution <- model_ve_runs %>% group_by (group, t) %>%
    summarise(VE_value = sample(VE_t, 1)) %>%
    ungroup()
}


VE_distribution %>% ggplot(aes(x = t, y = VE_value, col = group)) +
  geom_line()
