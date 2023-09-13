#'Experimental Design and Analysis
#'29/11/2021
#'
#'RANDOMISED COMPLETE BLOCK DESIGN for Kronos Greenhouse plants
#'
#'22/11/2021
#'Separate NAC-5A and NAP-A into two separate designs with block size 8.
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
require(lubridate)

#What about my own Kronos experiment
#This has 12 genotypes and 12 reps

#NAC-5A, 8 genotypes and 12 reps
set.seed(seed = 20211222)
treat<-paste(rep(c("X388","X392"),each=4), c("WT","aa","bb","aabb"), sep = "_")
#If the plant number for a given rep is updated, along with sowing date, we can add here.
#OR make another script that does it better
unrandomise <- data.frame(treat=rep(treat, each=12), rep=rep(1:12, 8), plant_num=1:96)

#RCBD
out<-design.rcbd(treat,r=12)
t(out$sketch)
out$book

Kronos_design <- out$book %>%
  group_by(treat) %>%
  mutate(rep = sample(12, 12, replace = FALSE)) %>%
  inner_join(unrandomise, by = c("treat", "rep")) %>%
  mutate(plant_name = paste("K", treat, plant_num, sep = "_"))

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse")
write_csv(Kronos_design, file = paste("RCBD_Kronos_NAC5_", today(), ".csv"), col_names = TRUE)

#NAP-A, 4 genotypes and 12 reps
set.seed(seed = 20211222)
treat<-paste(rep(c("X585", "X585"),each=4), c("WT","aa","bb","aabb"), sep = "_") #put identical treatments in twice for >1 treat per block
#If the plant number for a given rep is updated, along with sowing date, we can add here.
#OR make another script that does it better
unrandomise <- data.frame(treat=rep(treat, each=12), rep=rep(1:12, 4), plant_num=97:144)

#RCBD
out<-design.rcbd(treat,r=6)
t(out$sketch)
out$book

Kronos_design <- out$book %>%
  group_by(treat) %>%
  mutate(rep = sample(12, 12, replace = FALSE)) %>%
  inner_join(unrandomise, by = c("treat", "rep")) %>%
  mutate(plant_name = paste("K", treat, plant_num, sep = "_"))

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse")
write_csv(Kronos_design, file = paste("RCBD_Kronos_NAP_", today(), ".csv"), col_names = TRUE)
