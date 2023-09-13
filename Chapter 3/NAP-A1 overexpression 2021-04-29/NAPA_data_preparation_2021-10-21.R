#NAP-A phenotyping data preparation
#21/10/2021
#Open the Excel file, fill in blank columns
#'calculate dates as "days after sowing" or "days after heading"
#'and save as a .csv file.
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)


#'SOURCE section
#'Packages and source files
require(readxl); require(readr); require(tidyverse); require(lubridate)

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics needs tidying/NAC phenotyping")
NAPA_phenotyping_results_20210913 <- read_excel("NAPA_phenotyping_results_20210913.xlsx", 
                                                sheet = "NAPA_master_sheet")

#'EDIT section
#'Key parameters

sowing_date <- lubridate::ymd('2021-04-29')

#'RUN section
#'Explore the file and reformat data

summary(NAPA_phenotyping_results_20210913)

#Replace blanks from speedy table entry with appropriate values
NAPA_phenotyping_clean <- NAPA_phenotyping_results_20210913 %>%
  mutate(Height_total = if_else(is.na(Height_total), true=Height, false=Height_total),
         Snapped_3007 = if_else(Snapped_3007 %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
         Snapped_pre_peduncle_sen = if_else(Snapped_pre_peduncle_sen %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
         Tiller_0406 = replace_na(0),
         TKW = Grain_mass/Grain_number*1000)

#Calculate all dates as days since SOWING on 29/4/21 (NB plant numbers 11 & 12 sown on 30/4/21 but all calculated from 29/4)
#Using lubridate
#Also calculate Leaf and Peduncle Senescence as Days after Heading (DAH)
#Done slightly differently to NAC-5A data where heading was set as 'day of the year'
days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

NAPA_phenotyping_yday <- NAPA_phenotyping_clean %>%
  mutate(Leaf_senescence_DAH = days.after.sowing(Leaf_senescence, sowing_date = Heading_date),
         Peduncle_senescence_DAH = days.after.sowing(Peduncle_senescence, sowing_date = Heading_date),
         across(c(Heading_date, Leaf_senescence, Peduncle_senescence, Date_1, Date_2, Date_3, Date_4, Date_5, Heading_plus_twenty),
                ~ days.after.sowing(.x, sowing_date = sowing_date)))

write_csv(NAPA_phenotyping_clean, file = paste('NAPA_phenotyping_clean_', today(), '.csv', sep = ""))
write_csv(NAPA_phenotyping_yday, file = paste('NAPA_phenotyping_yday_', today(), '.csv', sep = ""))
