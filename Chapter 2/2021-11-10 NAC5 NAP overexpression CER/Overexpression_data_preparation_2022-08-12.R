#'Overexpression phenotyping data preparation
#'21/10/2021
#'Open the Excel file, fill in blank columns
#'calculate dates as "days after sowing" or "days after heading"
#'and save as a .csv file.
#'
#'07/03/2021
#'Adapt to use with Kronos NAC5 NAP data (incomplete). Remove unnecessary lines.
#'Using sowing date as column rather than fixed.
#'
#'29/03/2021
#'Use with Kronos_phenotyping_2022-03-28.xlsx
#'I have manually added the genotypes to the spreadsheet (so will need to sort this again when I get the final datafile)
#'
#'08/04/2022
#'Calculate metrics
#'Linear interpolation method Using the r function stats::approx
#'Based on NAM_calculate_metrics_2022-03-17
#'
#'27/06/2022
#'Use full dataset from 2022-04-28
#'
#'12/08/2022
#'Copy Kronos_data_preparation_2022-06-27
#'Use NAC5 NAP overexpression CER data
#'
#'05/09/2022
#'Join Marvin data
#'
#'29/09/2022
#'Use properly unique and adjusted Marvin data
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'SOURCE
#'Packages and source files
require(readxl); require(readr); require(lubridate); require(MESS) #for auc()
require(tidyverse)

#Run twice, once for each gene
gene = "NAC5"
# gene= "NAP"

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER")

#34 columns vs 33 columns
if(gene == "NAC5"){
Overexpression_phenotyping <- read_excel("2022-02-14_phenotyping_data/NAC5_NAP_CER_phenotyping_2022-04-28.xlsx",
                                 sheet = paste(gene),
                                 col_types = c("numeric", "numeric", "text", 
                                 "numeric", "numeric", "numeric", 
                                 "date", "numeric", "date", "date", 
                                 "numeric", "numeric", "date", "numeric", 
                                 "date", "numeric", "date", "numeric", 
                                 "date", "numeric", "date", "numeric", 
                                 "date", "numeric", "date", "numeric", 
                                  "numeric", "numeric", "skip", "numeric", 
                                 "numeric", "text", "date", "numeric"))
} else {
  Overexpression_phenotyping <- read_excel("2022-02-14_phenotyping_data/NAC5_NAP_CER_phenotyping_2022-04-28.xlsx",
                                           sheet = paste(gene),
                                           col_types = c("numeric", "numeric", "text", 
                                                         "numeric", "numeric", "numeric", 
                                                         "date", "numeric", "date", "date", 
                                                         "numeric", "numeric", "date", "numeric", 
                                                         "date", "numeric", "date", "numeric", 
                                                         "date", "numeric", "date", "numeric", 
                                                         "date", "numeric", "date", "numeric", 
                                                         "numeric", "skip", "numeric", 
                                                         "numeric", "text", "date", "numeric"))
}

Marvin_data <- read_csv("2022-08-12_clean_data/Overexpression_Marvin_unique_2022-09-29.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'EDIT
#'Key parameters

sowing_date <- lubridate::ymd('2021-12-22')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'TIDY
#'Explore the file and reformat data

summary(Overexpression_phenotyping)

# #Replace blanks from speedy table entry with appropriate values
# Overexpression_phenotyping_clean <- Overexpression_phenotyping %>%
#   mutate(Height_total = if_else(is.na(Height_total), true=Height, false=Height_total),
#          Snapped_3007 = if_else(Snapped_3007 %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
#          Snapped_pre_peduncle_sen = if_else(Snapped_pre_peduncle_sen %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
#          Tiller_0406 = replace_na(0),
#          TKW = Grain_mass/Grain_number*1000)

#Neaten various columns, add new columns, and remove photo columns
if(gene == "NAC5"){
  Overexpression_phenotyping_clean <- Overexpression_phenotyping %>%
    mutate(
      SPAD_8 = replace_na(Date_8, 0)
    )
} else if(gene == "NAP"){
  Overexpression_phenotyping_clean <- Overexpression_phenotyping %>%
    mutate(SPAD_8 = replace_na(SPAD_8, 0)) %>%
    rename(Date_6 = SPAD_52)
}

Overexpression_phenotyping_clean <- Overexpression_phenotyping_clean %>%
  arrange(Index) %>%
  mutate(Heading_date = lubridate::ymd(Heading_date),
         Tiller_2 = replace_na(Tiller_2, 0),
         Date_8 = Date_7 + days(7),
         Gene_name = gene,
         snapped_pre_ped_sen = snapped < Peduncle_senescence | !is.na(snapped) & is.na(Peduncle_senescence)) %>%
  tidyr::extract(Line_name, "Genotype", regex = "_([a-z]+)$", remove = FALSE)

not_needed <- c("photo", "Photo", "Column") #remove photo comments and unlabelled columns
Overexpression_phenotyping_clean <- Overexpression_phenotyping_clean %>%
  select(!contains(not_needed))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#CALCULATE
#Calculate all dates as days since SOWING
#Using lubridate
#Also calculate Leaf and Peduncle Senescence as Days after Heading (DAH)
#Done slightly differently to NAC-5A data where heading was set as 'day of the year'
days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

Overexpression_phenotyping_yday <- Overexpression_phenotyping_clean %>%
  mutate(Leaf_senescence_DAH = days.after.sowing(Leaf_senescence, sowing_date = Heading_date),
         Peduncle_senescence_DAH = days.after.sowing(Peduncle_senescence, sowing_date = Heading_date),
         Date_0 = Heading_date,
         across(c(starts_with("Date")),
                ~ days.after.sowing(.x, sowing_date = Heading_date)),
         across(c(Heading_date, Leaf_senescence, Peduncle_senescence),
                ~ days.after.sowing(.x, sowing_date = sowing_date)))



Long_data <- Overexpression_phenotyping_yday %>%
  rename(SPAD_0 = SPAD_heading) %>%
  pivot_longer(cols = starts_with("Date")|starts_with("SPAD"), names_pattern = "(.*)_([0-9])$", names_to = c("Category","Week")) %>%
  pivot_wider(names_from = Category, values_from = value)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#AUC NORMALISATION TEST

#Mean senescence values are dataset-specific as they depend on how many data are collected.
#Therefore if I want to compare between sowing dates I'll need to normalise to a range in days after sowing or heading

#normalise to a range in days after sowing or heading
Long_data %>% group_by(Week) %>%
  summarise(max(Date),
            min(Date),
            mean(Date)
  )
#End AUC at minimum scoring endpoint - ensures interpolation only, no extrapolation.
#Start point is less important as additional score points are 0 so do not add to total auc
auc_limit <- Long_data %>% filter(Week == max(Week)) %>%
  summarise(
    min_DAH = floor(min(Date))
  )

# Exampele <- Long_data %>%
#   group_by(Line_name, Plant_num, Index) %>%
#   summarise(
#     AUC_Leaf_Sen = auc(Date, SPAD, type = "linear")
#   ) %>%
#   pivot_wider(names_from = Plant_num, id_cols = Line_name, values_from = AUC_Leaf_Sen)
# 
# Exampele2 <- Long_data %>%
#   group_by(Line_name, Plant_num, Index) %>%
#   summarise(
#     AUC_Leaf_Sen = auc(Date, SPAD, type = "linear", to = auc_limit$min_DAH) #Add endpoint to auc()
#   ) %>%
#   pivot_wider(names_from = Plant_num, id_cols = Line_name, values_from = AUC_Leaf_Sen)

#Now all AUC values represent the same time interval, between 0% senescence and min endpoint
#Here 55 days
#Differences are subtle as SPAD close to 0 at endpoint anyway



#Calculate metrics from SPAD table
NA_count <- Long_data %>%
  filter(!is.na(Plant_num)) %>%
  group_by(Line_name, Plant_num) %>%
  arrange(-(Date)) %>%
  summarise(
    NA_count = sum(is.na(SPAD))
  )

Summary_Days_after_heading <- inner_join(Long_data, NA_count, by = c("Line_name", "Plant_num")) %>%
  filter(!is.na(Plant_num) & NA_count < 4) %>%
  group_by(Line_name, Plant_num) %>%
  arrange(-(Date)) %>%
  summarise(
    TT_SPAD10 = approx(SPAD, Date, c(10), ties = "ordered")$y,
    TT_SPAD30 = approx(SPAD, Date, c(30), ties = "ordered")$y,
    TT_SPAD40 = approx(SPAD, Date, c(40), ties = "ordered")$y,
    Dur_Leaf_Sen = TT_SPAD10 - TT_SPAD40,
    AUC_Leaf_Sen = auc(Date, SPAD, type = "linear")
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#JOIN
#
#Join calculated values AND Marvin data

Overexpression_phenotyping_yday_v2 <- left_join(Overexpression_phenotyping_yday, NA_count, by = c("Line_name", "Plant_num"))

Overexpression_phenotyping_yday_v3 <- left_join(Overexpression_phenotyping_yday_v2, Summary_Days_after_heading, by = c("Line_name", "Plant_num"))

Overexpression_phenotyping_yday_v4 <- Overexpression_phenotyping_yday_v3 %>%
  mutate(ID = paste(Line_name, Plant_num, sep = "_")) %>%
  left_join(Marvin_data, by = c("ID"))

# Example <- SPAD_table %>%
#   filter(Plant_num == 4) %>%
#   select(Plant_num, Date, SPAD) %>%
#   arrange(-(Date))
# approx(Example$SPAD, Example$Date, c(40), ties = "ordered")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#WRITE

setwd("./2022-08-12_clean_data")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Overexpression greenhouse/2022-03-07_clean_data")
write_csv(Overexpression_phenotyping_clean, file = paste('Overexpression_phenotyping_clean_', gene, "_", today(), '.csv', sep = ""))
write_csv(Overexpression_phenotyping_yday_v4, file = paste('Overexpression_phenotyping_yday_', gene, "_", today(), '.csv', sep = ""))
write_csv(Long_data, file = paste('Overexpression_SPAD_yday_', gene, "_", today(), '.csv', sep = ""))

