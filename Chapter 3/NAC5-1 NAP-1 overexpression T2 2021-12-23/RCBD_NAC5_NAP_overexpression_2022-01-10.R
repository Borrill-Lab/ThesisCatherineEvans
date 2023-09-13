#'Experimental Design and Analysis
#'10/01/2022
#'
#'RANDOMISED COMPLETE BLOCK DESIGN for NAC5 NAP overexpression CER plants
#'
#'From
#'Practical 4
#'Experimental Design using R
#'
#'Made in R 4.0.4
#'
#'

#Load packages
require(agricolae)
require(tidyverse)

#NAC5
#This has 8 genotypes and 10 reps
set.seed(seed = 20220110)
treat<-c("5.4_hom", "5.4_null", "8.2_multi", "8.23_hom", "8.23_null", "8.3_multi", "8.9_hom", "8.9_null")
out<-design.rcbd(treat,r=10)
t(out$sketch)
out$book

NAC5_design <- out$book %>%
  group_by(treat) %>%
  mutate(plant_num = sample(10, 10, replace = FALSE), plant_name = paste(treat, plant_num, sep = "_"), gene = "NAC5")

set.seed(seed = 20220111)
treat<-c("7.19_hom", "7.19_null", "7.9_multi", "10.5_hom", "10.5_null", "10.1_multi", "10.9_hom", "10.9_null")
out<-design.rcbd(treat,r=10)
t(out$sketch)
out$book

NAP_design <- out$book %>%
  group_by(treat) %>%
  mutate(plant_num = sample(10, 10, replace = FALSE), plant_name = paste(treat, plant_num, sep = "_"), gene = "NAP")


setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER")
write_csv(rbind(NAC5_design, NAP_design), file = "RCBD_CER_2022-01-10.csv", col_names = TRUE)