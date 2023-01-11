rm(list=ls())

################
## Directories
################

NZdata_dir <- file.path(getwd(), "NZ_data")

fig_dir <- file.path(NZdata_dir, "Figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

library(tidyverse)
library(proj4)
library(akima)

# ################
## Load data
################
obsfull <- readRDS(file.path(NZdata_dir, "NZ_observations.rds"))
netfull <- readRDS(file.path(NZdata_dir, "NZ_network.rds"))


# NZFFD_Waitaki_Waikato <- obsfull %>% filter(data_type=="encounter", c(grepl("aitaki", CatName) | grepl("aikato", CatName)))
# 
# WRC_Waikato_count <- obsfull %>% filter(source=="waikato_study")


