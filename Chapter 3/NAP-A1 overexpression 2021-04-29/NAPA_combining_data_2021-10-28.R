#NAP-A phenotyping data preparation
#28/10/2021
#'After NAPA_data_preparation_2021-10-21.R
#'Combine ALL the data we have about each plant into one table.
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)


#'SOURCE section
#'Packages and source files
require(readxl); require(readr); require(tidyverse); require(lubridate); require(reshape2)

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics needs tidying/NAC phenotyping")
NAPA_phenotyping_yday_2021_10_25 <- read_csv("NAPA_phenotyping_yday_2021-10-29.csv")
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics needs tidying/NAC expression R")
NAPA_genotype_key <- read_excel("NAPA_genotype_key.xlsx")
qPCR_NAPA_individual_analysis_pfaffl_2021_10_28 <- read_csv("qPCR_NAPA_individual_analysis_pfaffl_2021-10-28.csv")

#Split Sample into Line name and Plant num
#Pivot so that there is one row per plant
qPCR_table_pivot <- qPCR_NAPA_individual_analysis_pfaffl_2021_10_28 %>%
  mutate(Sample2 = Sample) %>%
  separate(Sample2, into=c("Line_name", "Plant_num"), sep = "-") %>%
  dplyr::select(Sample, Line_name, Plant_num, Primer.y, fold_change) %>%
  filter(!is.na(Primer.y)) %>%
  pivot_wider(names_from = Primer.y, values_from = fold_change)
  

#Add everything together!! 144 plants * 36 variables = 5184 pieces of data
everything_table <- NAPA_phenotyping_yday_2021_10_25 %>%
  mutate(Plant_num = as.character(Plant_num)) %>%
  left_join(NAPA_genotype_key, by = "Line_name") %>%
  left_join(qPCR_table_pivot, by = c("Line_name", "Plant_num"))

#Let's also create a SPAD table with one row per SPAD value
SPAD_table <- everything_table %>%
  filter(Priority == TRUE) %>%
  rename(SPAD_0 = SPAD_heading, Date_0 = Heading_date) %>%
  pivot_longer(cols = c(5:6,12:21), names_pattern = "(.*)_([0-9])$", names_to = c("Category","Week")) %>%
  pivot_wider(names_from = Category, values_from = value)


setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics needs tidying/NAC phenotyping")
write_csv(everything_table, file = paste('NAPA_EVERYTHING_', today(), '.csv', sep = ""))
write_csv(SPAD_table, file = paste('NAPA_SPAD_', today(), '.csv', sep = ""))
