# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'NAM2 phenotyping data preparation
#'Base on NAM2 phenotyping data preparation
#'
#'29/03/2021
#'Use with NAM2_phenotyping_2022-03-28.xlsx
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
#'12/01/2023
#'Adapt to use with NAM-2 greenhouse data. Adapt to run all together once for NAM2 X50 and NAM2 X61.
#'Sheet formatted differently to NAC5/NAP datasets.
#'
#'15/06/2023
#'Use with NAM-2 Exp2: NAM2_Exp2_phenotyping_results_2023-06-12.xlsx
#'Fixed A / B / D genotypes
#'
#'30/06/2023
#'Added Week_approx and Dosage as variables
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
get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-04-13_phenotyping_data", sep = ""))

#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date the scripts were run through NAM2_data_preparation.R
gene = "NAM2"
date = "2023-06-12"

NAM2_phenotyping_X50 <- read_excel(paste("NAM2_Exp2_phenotyping_results", "_", date, ".xlsx", sep= ""), sheet = "X50",
                                   col_types = c("numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "text",  "numeric", 
                                                 "text", "text", "text", 
                                                 "text", "text", "numeric", "date", "text",
                                                 "numeric", "date", "date", "date", 
                                                 "numeric", "numeric", "text", "date", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "text", "text", "text", "text"))

NAM2_phenotyping_X61 <- read_excel(paste("NAM2_Exp2_phenotyping_results", "_", date, ".xlsx", sep= ""), sheet = "X61",
                                   col_types = c("numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "text",  "numeric", 
                                                 "text", "text", "text", 
                                                 "text", "text", "numeric", "date", "text",
                                                 "numeric", "date", "date", "date", 
                                                 "numeric", "numeric", "text", "date", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "text", "text", "text", "text"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'EDIT
#'Key parameters

Sowing_date <- lubridate::ymd('2023-02-15') #Sowing date the same for NAM2 data

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'TIDY
#'Explore the file and reformat data

NAM2_phenotyping <- full_join(NAM2_phenotyping_X50, NAM2_phenotyping_X61)

summary(NAM2_phenotyping)

#Replace blanks from speedy table entry with appropriate values
NAM2_phenotyping_clean <- NAM2_phenotyping %>%
  mutate(Snapped_pre_peduncle_sen = if_else(Snapped_pre_peduncle_sen %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
         Heading_backdated = if_else(Heading_backdated %in% c('T','t',TRUE,'TRUE'), true = TRUE, false = FALSE),
         Chlorosis = replace_na(Chlorosis, "none"),
         across(c(`SPAD_2023-06-02`, `SPAD_2023-06-06`, `SPAD_2023-06-09`, `SPAD_2023-06-13`), ~ replace_na(.x, 0))) %>%
  select(where(~ !all(is.na(.x)))) #Remove columns that are totally blank

#Neaten various columns, add new columns, and remove photo columns
NAM2_phenotyping_clean <- NAM2_phenotyping_clean %>%
  arrange(Index) %>%
  mutate(Heading_date = lubridate::ymd(Heading_date),
         Tiller_2 = replace_na(Tiller_2, 0),
         #THESE LINES ARE COMPLETELY WRONG AND NEED TO BE SORTED IN THE OTHER ANALYSIS TOO
         # A = na.omit(c(`CE1656-COM-2`, `CE1773-2`)),
         # B = na.omit(c(`CE1436-COM-3`, `CE0895-COM-2`)),
         # D = na.omit(c(`CE0231-COM-2`, `CE0237`)),
         A = if_else(is.na(`CE1656-COM-2`), true = `CE1773-2`, false = `CE1656-COM-2`),
         B = if_else(is.na(`CE1436-COM-3`), true = `CE0895-COM-2`, false = `CE1436-COM-3`),
         D = if_else(is.na(`CE0231-COM-2`), true = `CE0237`, false = `CE0231-COM-2`),
         Column = (Block*2-1) + (Position - 1) %/% 4,
         # Tray_Row = (Plant_num - 1) %/% 12 + 1, #NOT sure if these trays are correct
         #  Tray_Column = (Plant_num - 1) %% 12 + 1,
         #  Tray = (Plant_num - 1) %/% 96 + 1
         ) %>%
  tidyr::extract(Line_name, c("Cross", "Genotype"), regex = "^(X[0-9]{2})_(G[0-9])", remove = FALSE) %>%
  mutate(Cross = ifelse(Cross == "X54", "X50", Cross),
         Dosage = fct_recode(Genotype,
                                    "Wild type" = "G1",
                                    "single" = "G5",
                                    "single" = "G3",
                                    "single" = "G2",
                                    "double" = "G7",
                                    "double" = "G6",
                                    "double" = "G4",
                                    "triple" = "G8")
         )

not_needed <- c("photo", "Photo", "CE0", "CE1", "Pick?", "Hom") #remove photo comments, unlabelled columns, and surplus columns from KASP data
NAM2_phenotyping_clean <- NAM2_phenotyping_clean %>%
  dplyr::select(!contains(not_needed))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#CALCULATE
#Calculate all dates as days since SOWING
#Using lubridate
#Also calculate Leaf and Peduncle Senescence as Days after Heading (DAH)

days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

NAM2_phenotyping_yday <- NAM2_phenotyping_clean %>%
  mutate(Leaf_senescence_DAH = days.after.sowing(Leaf_senescence, sowing_date = Heading_date),
         Leaf_senescence_2nd_DAH = days.after.sowing(Leaf_senescence_2nd, sowing_date = Heading_date),
         Peduncle_senescence_DAH = days.after.sowing(Peduncle_senescence, sowing_date = Heading_date),
         Date_0 = Heading_date,
         across(c(starts_with("Date")),
                 ~ days.after.sowing(.x, sowing_date = Heading_date))
  )


# #Remove these cols from SPAD table
# NAM2_phenotyping_v3 <- NAM2_phenotyping_yday %>%
#   dplyr::select(!`SPAD_2022-11-26 and 23` & !`SPAD_2022-12-16_repeat` & !`SPAD_2022-10-28`)
NAM2_phenotyping_v3 <- NAM2_phenotyping_yday

SPAD_table <- NAM2_phenotyping_v3 %>%
  dplyr::select(
    Index,
    Block,
    Position,
    Row,
    Line_name,
    Cross, Genotype, Dosage,
    Plant_num, Rep, Sample,
    A, B, D,
    Heading_date,
    Scoring_day,
    Date_0,
    Date_1,
    all_of(starts_with("SPAD"))
  ) %>%
  pivot_longer(cols = starts_with("SPAD"), names_pattern = "SPAD_(.*)$", names_to = c("Date"), values_to = "SPAD")

SPAD_table <- SPAD_table %>%
  mutate(Date = ifelse(Date == "heading", as.character(Heading_date), Date),
         Date = ymd(Date),
         DAH = days.after.sowing(Date, sowing_date = Heading_date),
         Week_approx = DAH %/% 7
         ) %>%
  filter(!is.na(SPAD))


#Calculate metrics from SPAD table

#I think the issue here is it's taking the closest values of SPAD for the interpolation
#When there's more than one Zero, it may take any of the SPAD = 0 values to compare
#Removing all but the earliest Zero value could help with this
#In the meantime just ignore TTSPAD10 in this exploratory analysis

Summary_Days_after_heading <- SPAD_table %>%
  filter(!is.na(Plant_num)) %>%
  #filter(!(Week == 0 & SPAD < 40)) %>% #fudge factor to sort the two plants with a SPAD<40 at heading, which increases after
  group_by(Plant_num) %>%
  arrange(-(DAH)) %>%
  summarise(
    TT_SPAD10 = approx(SPAD, DAH, c(10), ties = list("ordered", min))$y,
    TT_SPAD30 = approx(SPAD, DAH, c(30), ties = list("ordered", min))$y,
    TT_SPAD40 = approx(SPAD, DAH, c(40), ties = list("ordered", min))$y,
    Dur_Leaf_Sen = TT_SPAD10 - TT_SPAD40,
    Rate_Leaf_Sen = 30/Dur_Leaf_Sen,
    AUC_Leaf_Sen = auc(DAH, SPAD, type = "linear")
  )

NAM2_phenotyping_yday <- inner_join(NAM2_phenotyping_yday, Summary_Days_after_heading, by = "Plant_num") %>%
  mutate(
    across(c(Heading_date, Leaf_senescence, Peduncle_senescence),
           ~ days.after.sowing(.x, sowing_date = Sowing_date))
  )

Example <- SPAD_table %>%
  filter(Plant_num == 91) %>%
  dplyr::select(Plant_num, DAH, SPAD) %>%
  arrange(-(DAH))
approx(Example$SPAD, Example$DAH, c(10), ties = c("ordered", min))
approx(Example$SPAD, Example$DAH, c(10), ties = c("ordered", max))

Example <- SPAD_table %>%
  filter(Plant_num == 253) %>%
  dplyr::select(Plant_num, DAH, SPAD) %>%
  arrange(-(DAH))
approx(Example$SPAD, Example$DAH, c(40), ties = c("ordered", min))
approx(Example$SPAD, Example$DAH, c(40), ties = c("ordered", max))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#WRITE

setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_clean_data", sep = ""))
write_csv(NAM2_phenotyping_clean, file = paste('NAM2_phenotyping_clean_', today(), '.csv', sep = ""))
write_csv(NAM2_phenotyping_yday, file = paste('NAM2_phenotyping_yday_', today(), '.csv', sep = ""))
write_csv(SPAD_table, file = paste('NAM2_SPAD_yday_', today(), '.csv', sep = ""))
