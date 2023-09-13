#'Kronos phenotyping data preparation
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
#'#07/09/2022
#'Oh dear this wasn't the most recent one, never mind
#'Add Marvin and NIR data
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'SOURCE
#'Packages and source files
require(readxl); require(readr); require(lubridate); require(MESS) #for auc()
require(tidyverse)

#Run twice, once for each gene
gene = "NAC5"
# gene= "NAP"

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse")
Kronos_phenotyping <- Kronos_phenotyping <- read_excel("2022-02-09_phenotyping_data/Kronos_phenotyping_2022-04-28.xlsx",
                                                       sheet = paste("Kronos", gene),
                                                       col_types = c("numeric", "numeric", "text", 
                                                                     "numeric", "date", "numeric", "date", 
                                                                     "date", "numeric", "numeric", "date", 
                                                                     "numeric", "date", "numeric", "date", 
                                                                     "numeric", "date", "numeric", "date", 
                                                                     "numeric", "date", "numeric", "date", 
                                                                     "numeric", "date", "numeric", "numeric", 
                                                                     "numeric", "numeric", "numeric", 
                                                                     "numeric", "date", "text", "text", 
                                                                     "numeric", "text", "numeric", "text", 
                                                                     "numeric", "numeric", "text", "text"))

Kronos_Marvin_clean <- read_csv("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data/Kronos_Marvin_clean.csv")

Kronos_NIR_Extensive_Report_Wheat_clean <- read_csv("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data/Kronos_NIR_Extensive Report_Wheat_clean_2022-09-07.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'EDIT
#'Key parameters

# sowing_date <- lubridate::ymd('2021-12-01') #Sowing date not all the same for Kronos data

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'TIDY
#'Explore the file and reformat data

summary(Kronos_phenotyping)

# #Replace blanks from speedy table entry with appropriate values
# Kronos_phenotyping_clean <- Kronos_phenotyping %>%
#   mutate(Height_total = if_else(is.na(Height_total), true=Height, false=Height_total),
#          Snapped_3007 = if_else(Snapped_3007 %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
#          Snapped_pre_peduncle_sen = if_else(Snapped_pre_peduncle_sen %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
#          Tiller_0406 = replace_na(0),
#          TKW = Grain_mass/Grain_number*1000)

#Neaten various columns, add new columns, and remove photo columns
Kronos_phenotyping_clean <- Kronos_phenotyping %>%
  arrange(Index) %>%
  mutate(Heading_date = lubridate::ymd(Heading_date),
         Height_total = if_else(is.na(Height_total), true=Height, false=Height_total),
         Tiller_2 = replace_na(Tiller_2, 0),
         Date_9 = Date_8 + days(7),
         SPAD_9 = replace_na(SPAD_9, 0),
         A = factor(sub("\\?", NA, A)),
         B = factor(sub("\\?", NA, B)),
         Genotype = interaction(A, B)) %>%
  tidyr::extract(Line_name, "Cross", regex = "^(X[0-9]{3})", remove = FALSE)

if(gene == "NAP"){
  Kronos_phenotyping_clean <- Kronos_phenotyping_clean %>%
  tidyr::unite(col = "Comments", Column1, Column2, sep = ", ", remove = TRUE, na.rm = TRUE)
}

not_needed <- c("photo", "Photo", "Column") #remove photo comments and unlabelled columns
Kronos_phenotyping_clean <- Kronos_phenotyping_clean %>%
  dplyr::select(!contains(not_needed))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#CALCULATE
#Calculate all dates as days since SOWING
#Using lubridate
#Also calculate Leaf and Peduncle Senescence as Days after Heading (DAH)
#Done slightly differently to NAC-5A data where heading was set as 'day of the year'
days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

Kronos_phenotyping_yday <- Kronos_phenotyping_clean %>%
  mutate(Leaf_senescence_DAH = days.after.sowing(Leaf_senescence, sowing_date = Heading_date),
         Peduncle_senescence_DAH = days.after.sowing(Peduncle_senescence, sowing_date = Heading_date),
         Date_0 = Heading_date,
         across(c(starts_with("Date")),
                ~ days.after.sowing(.x, sowing_date = Heading_date)),
         across(c(Heading_date, Leaf_senescence, Peduncle_senescence),
                ~ days.after.sowing(.x, sowing_date = Sowing_date)))



SPAD_table <- Kronos_phenotyping_yday %>%
  rename(SPAD_0 = SPAD_heading) %>%
  pivot_longer(cols = starts_with("Date")|starts_with("SPAD"), names_pattern = "(.*)_([0-9])$", names_to = c("Category","Week")) %>%
  pivot_wider(names_from = Category, values_from = value)


#Calculate metrics from SPAD table
Summary_Days_after_heading <- SPAD_table %>%
  filter(!is.na(Plant_num)) %>%
  filter(!(Week == 0 & SPAD < 40)) %>% #fudge factor to sort the two plants with a SPAD<40 at heading, which increases after
  group_by(Plant_num) %>%
  arrange(-(Date)) %>%
  summarise(
    TT_SPAD10 = approx(SPAD, Date, c(10), ties = "ordered")$y,
    TT_SPAD30 = approx(SPAD, Date, c(30), ties = "ordered")$y,
    TT_SPAD40 = approx(SPAD, Date, c(40), ties = "ordered")$y,
    Dur_Leaf_Sen = TT_SPAD10 - TT_SPAD40,
    Rate_Leaf_Sen = 30/Dur_Leaf_Sen,
    AUC_Leaf_Sen = auc(Date, SPAD, type = "linear")
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#JOIN
Kronos_phenotyping_yday_v2 <- left_join(Kronos_phenotyping_yday, Summary_Days_after_heading, by = "Plant_num")

#Choose the most recent Marvin value where there are two
Kronos_Marvin_sorted <- Kronos_Marvin_clean %>%
  mutate(Date_Marvin = lubridate::dmy(Date_Marvin)) %>%
  arrange(desc(Date_Marvin))

Kronos_pairs <- Kronos_Marvin_sorted[duplicated(Kronos_Marvin_sorted$ID),]
Kronos_Marvin_unique <- Kronos_Marvin_sorted[!duplicated(Kronos_Marvin_sorted$ID),]


Kronos_phenotyping_yday_v3 <- left_join(Kronos_phenotyping_yday_v2, Kronos_Marvin_unique, by = c("Plant_name" = "ID"))

Kronos_phenotyping_yday_v4 <- left_join(Kronos_phenotyping_yday_v3, Kronos_NIR_Extensive_Report_Wheat_clean, by = c("Plant_name" = "Sample ID"))

# Example <- SPAD_table %>%
#   filter(Plant_num == 4) %>%
#   select(Plant_num, Date, SPAD) %>%
#   arrange(-(Date))
# approx(Example$SPAD, Example$Date, c(40), ties = "ordered")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#WRITE

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data")
write_csv(Kronos_phenotyping_clean, file = paste('Kronos_phenotyping_clean_', gene, "_", today(), '.csv', sep = ""))
write_csv(Kronos_phenotyping_yday_v4, file = paste('Kronos_phenotyping_yday_', gene, "_", today(), '.csv', sep = ""))
write_csv(SPAD_table, file = paste('Kronos_SPAD_yday_', gene, "_", today(), '.csv', sep = ""))
