#'Experimental Design and Analysis
#'14/09/22
#'
#'RANDOMISED COMPLETE BLOCK DESIGN for NAM2 Greenhouse Experiment
#'
#'21/10/22
#'How can we account for light in a balanced way?
#'I am using the same seed for everything 
#'
#'16/03/23
#'Repeat with different seed for Spring 2023 experiment
#'
#'From
#'Practical 4
#'Experimental Design using R
#'
#'Made in R 4.0.4
#'
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Load packages
require(agricolae)
require(tidyverse)
require(lubridate)

#OPEN
#Open file with plant numbers
library(readxl)
KASP_NAM2_Greenhouse_calls_set1 <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/KASP NAM2 Exp2 2023-03-03/2023-03-10 _KASP_NAM2_Exp2_calls.xlsx",
                                              sheet = "set1")
plants_to_pick_set1 <- KASP_NAM2_Greenhouse_calls_set1 %>%
  tidyr::extract(Sample, into = "treat", regex = "(X[0-9]{2}-G[0-9])", remove = FALSE) %>%
  mutate(treat = gsub("-", "_", treat)) %>%
  mutate(Rep = rep(1:10, 8))

KASP_NAM2_Greenhouse_calls_set2 <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/KASP NAM2 Exp2 2023-03-03/2023-03-10 _KASP_NAM2_Exp2_calls.xlsx",
                                              sheet = "set2")
plants_to_pick_set2 <- KASP_NAM2_Greenhouse_calls_set2 %>%
  tidyr::extract(Sample, into = "treat", regex = "(X[0-9]{2}-G[0-9])", remove = FALSE) %>%
  mutate(treat = gsub("-", "_", treat)) %>%
  mutate(Rep = rep(1:10, 8))

#NAM2, 16 genotypes and 10 reps
set.seed(seed = 20230316) #Different seed
treat <- c(paste(rep(c("X50_G","X61_G"),each=8), c(1:8), sep = ""))
treat_set1 <- treat[c(1:8)]
treat_set1[7] <- "X54_G7"
treat_set2 <- treat[c(9:16)]

# #OR make another script that does it better
n_reps <- 10
n_treat <- length(treat)
# unrandomise <- data.frame(treat=rep(treat, each=n_reps), rep=rep(1:n_reps, length(treat)), plant_num=1:(length(treat)*n_reps))



# # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Could a Latin Square Based Design help with multiple parameters??
str(design.lsd) #this seems to list the function parameters

latin_square <- design.lsd(treat_set1, seed = 20230316) #auto reps = n_treat = 8

View(t(latin_square$sketch))

Bonus_RCBD <- design.rcbd(treat_set1,r=2, seed = 20230316)

Bonus_RCBD_table <- Bonus_RCBD$book %>%
  mutate(plots = plots + 800,
         block = as.factor(as.numeric(block) + 8))

Latin_square_plus_set1 <- full_join(latin_square$book, Bonus_RCBD_table) %>%
  group_by(treat_set1) %>%
  mutate(Rep = sample(n_reps, n_reps, replace = FALSE)) %>%
  inner_join(plants_to_pick_set1, by = c("treat_set1"="treat", "Rep")) %>%
  select(plots, row, col, Plant_num, everything())

latin_square <- design.lsd(treat_set2, seed = 20230317) #auto reps = n_treat = 8

View(t(latin_square$sketch))

Bonus_RCBD <- design.rcbd(treat_set2,r=2, seed = 20230317)

Bonus_RCBD_table <- Bonus_RCBD$book %>%
  mutate(plots = plots + 800,
         block = as.factor(as.numeric(block) + 8))

Latin_square_plus_set2 <- full_join(latin_square$book, Bonus_RCBD_table) %>%
  group_by(treat_set2) %>%
  mutate(Rep = sample(n_reps, n_reps, replace = FALSE)) %>%
  inner_join(plants_to_pick_set2, by = c("treat_set2"="treat", "Rep")) %>%
  select(plots, row, col, Plant_num, everything())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#check if the conditions are met - considering 10 blocks of 8
p<-choose(8,2); p #28 unique pairwise comparisons between treatments
b<-choose(8,2); b #28 pairwise comparisons within a block
b*10/p #lambda = 10 of each unique pairwise comparison across the study. This must be an integer!

#check if the conditions are met, ish - considering 4 rows of 20
p<-choose(8,2); p #28 unique pairwise comparisons between treatments
b<-choose(20,2); b #190 pairwise comparisons within a block
b*4/p #lambda = 27.14 of each unique pairwise comparison across the study. Not an integer so not quite balanced re light

p<-choose(2,2); p #1 unique pairwise comparison for a single "factor"
b<-choose(20,2); b #190 pairwise comparisons within a block
b*4/p #lambda = 760 of each unique pairwise comparison across the study. Fine when I consider a single homoeolog at a time

#check if the conditions are met, ish - considering 40 row*block combos of 2
p<-choose(8,2); p #28 unique pairwise comparisons between treatments
b<-choose(2,2); b #1 pairwise comparison within a block
b*40/p #lambda = 1.43 of each unique pairwise comparison across the study. Not an integer so not quite balanced re light


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #RCBD
# out<-design.rcbd(treat,r=n_reps)
# t(out$sketch)
# out$book
# 
# Full_RCBD <- out$book
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #Half RCBD
# set.seed(seed = 20220915)
# out<-design.rcbd(treat_set1,r=n_reps)
# View(t(out$sketch))
# out$book
# 
# Half_RCBD_set1 <- out$book %>%
#   group_by(treat_set1) %>%
#   mutate(Rep = sample(n_reps, n_reps, replace = FALSE)) %>%
#   inner_join(plants_to_pick_set1, by = c("treat_set1"="treat", "Rep"))
# 
# set.seed(seed = 20220915)
# out<-design.rcbd(treat_set2,r=n_reps)
# View(t(out$sketch))
# out$book
# 
# Half_RCBD_set2 <- out$book

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#BIBD
#Using agricolae package, which helpfully gives lambda and efficiency stats
# #BIBD  blocks of 8
# set.seed(seed = 20220915)
# treat <- c(paste(rep(c("X50_G","X61_G"),each=8), c(1:8), sep = ""), "Cadenza_a", "Cadenza_b")
# bibd<-design.bib(treat_set1, k = 8) #Is invalid! This shows this design would not have a balanced number of pairwise comparisons.
# t(bibd$sketch)

# #NAM2, 18 genotypes and 10 reps
# set.seed(seed = 20220915) #Different seed to field experiment
# treat <- c(paste(rep(c("X50_G","X61_G"),each=8), c(1:8), sep = ""), "Cadenza_a", "Cadenza_b")
# treat_set1 <- treat[c(1:8, 17)]
# treat_set2 <- treat[c(9:16, 18)]
#
# #Trays with plants in sets of 9
# set.seed(seed = 20220916) #Different seed
# treat <- c(paste(rep(c("X50_G","X61_G"),each=8), c(1:8), sep = ""))
# treat_set1 <- treat[c(1:8)]
# treat_set2 <- treat[c(9:16)]
# n_reps = 2
# out<-design.rcbd(treat_set1,r=n_reps)
# View(t(out$sketch))
# View(out$book[3])
# out<-design.rcbd(treat_set2,r=n_reps)
# View(t(out$sketch))
# View(out$book[3])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#WRITE
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING")
write_csv(Latin_square_plus_set1, file = paste("Latinplus_NAM2_set1", today(), ".csv"), col_names = TRUE)
write_csv(Latin_square_plus_set2, file = paste("Latinplus_NAM2_set2", today(), ".csv"), col_names = TRUE)
