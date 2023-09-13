#'04/03/21
#Analyse the NAC-5A senescence data
#Work out how to deal with dates
#'lubridate is a tidyverse package to deal with dates and times.
#'https://r4ds.had.co.nz/dates-and-times.html 
#'19/03/21
#'Edit: use file NAC5A_phenotype_clean_CE_190321.csv

library(readr)
library(ggplot2)
library(lubridate)
library(reshape2)
theme_set(theme_bw())
setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics/NAC phenotyping")
Table <- read_csv(file = "NAC5A_phenotype_clean_CE_190321.csv")

#Parse strings of dates as dmy() mdy() or ymd()
#yday gives the day of the year, from 2021-01-01 = 1
head(dmy(Table$Heading_date))
Table_v2 <- Table
#Change all date columns to day of the year
Table_v2$Heading_date <- yday(dmy(Table$Heading_date))
Table_v2$Flag_leaf_25 <- yday(dmy(Table$Flag_leaf_25))
Table_v2$Peduncle_100 <- yday(dmy(Table$Peduncle_100))
Table_v2$First_score_date <- yday(dmy(Table$First_score_date))

Table_v2$First_score_date-Table_v2$Heading_date==Table_v2$DAH #TRUE, as DAH is the number of days between Heading and First Score
#Add columns for flag leaf senescence and peduncle senescence in DAH by subtracting heading date
Table_v2 <- cbind(Table_v2, Flag_leaf_25_DAH = Table_v2$Flag_leaf_25-Table_v2$Heading_date, Peduncle_100_DAH = Table_v2$Peduncle_100-Table_v2$Heading_date)

write_csv(Table_v2, file="NAC5A_phenotype_yday_190321.csv")
