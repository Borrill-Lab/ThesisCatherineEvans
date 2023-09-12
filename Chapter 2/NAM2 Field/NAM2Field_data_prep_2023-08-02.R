# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'21/02/2022
#'Data preparation
#'
#'I need to
# 1.	Extract the raw data from the spreadsheets
# 3.	Combine together
# 4.	Make a “wide format” table with 1 plot per row and a “long format” table with 1 plot * 1 timepoint per row.
#'5. Calculate time in "days since ear emergence" and "thermal time since ear emergence"
#'
#'Based on NAM_data_prep_2022-07-19
#'
#'02/08/2023
#'Adapt for NAM2 field data from 2023
#'
#'03/08/2023
#'Add temperature Degree_days_after_heading
#'
#'04/08/2023
#'Add Dosage and A.B.D for formatting and factorial models
#'
#'#'Built with R 4.0.4 and tidyverse
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(tidyverse)
require(readxl)
require(lubridate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN file
get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()

setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/2023-03-07_raw_data", sep = ""))

table_name <- "Field_results_NAM2_"

Field_results <- read_excel("NAM2_Field23_phenotyping_results_neat_2023-08-02.xlsx", skip = 4)

#Temperature data from Church Farm 2023
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23", sep = ""))
temperature_by_day <- read_csv("Weather_data/daily_mean_temperature_2023.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'Calculates the difference between 2 dates in degree days
#'EXCLUDES the temperature on the start date
#'So that if you put in the same date for start and end, the answer is 0
#'Half days are currently rounded UP
#'EDIT: Capitals for Column Headings
degree_day_interval <- function(start_date, end_date){
  require(lubridate); require(tidyverse)
  if(length(start_date)==length(end_date)){
    degree_days = c()
    for(i in 1:length(start_date)){
      focus_interval <- interval(start = start_date[i]+days(1), end = end_date[i], tzone = "GMT")
      degree_days[i] <- temperature_by_day %>%
        filter(Date %within% focus_interval) %>%
        summarise(sum(Mean_temp)) %>%
        unlist()
    }
    return(degree_days)
  }else{
    print("Error: Start date and end date need to be vectors of the same length")
  }
  
}

days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY DATA
summary(Field_results)

#Completeness
Wide_data <- Field_results %>%
  mutate(
    `Leaf_senescence.2023-06-09` = 0,
    `Peduncle_senescence.2023-07-27` = 10,
    across(starts_with("Leaf_senescence"),
           ~ signif(.x * 2, 2)/2
    )
  ) %>%
  unite(Comment, all_of(starts_with("Comment")), sep = ",", na.rm = TRUE)

#Regex generalized for Cadenza inclusion
Wide_data_v2 <- Wide_data %>%
  rename(Line_name = "Genotype") %>%
  extract(Line_name, into = c("Cross", "Genotype"), "([[:alnum:]]+)_([[:alnum:]]+)", remove = FALSE) %>%
  mutate(Cross = ifelse(Cross == "X54", "X50", Cross),
         Dosage = fct_recode(Genotype,
                      "Wild type" = "G1",
                      "single" = "G5",
                      "single" = "G3",
                      "single" = "G2",
                      "double" = "G7",
                      "double" = "G6",
                      "double" = "G4",
                      "triple" = "G8",
                      "Wild type" = "WT"),
         A = fct_recode(Genotype,
                             "AA" = "G1",
                             "AA" = "G2",
                             "AA" = "G3",
                             "AA" = "G4",
                             "aa" = "G5",
                             "aa" = "G6",
                             "aa" = "G7",
                             "aa" = "G8",
                        "AA" = "WT"),
         B = fct_recode(Genotype,
                        "BB" = "G1",
                        "BB" = "G2",
                        "bb" = "G3",
                        "bb" = "G4",
                        "BB" = "G5",
                        "BB" = "G6",
                        "bb" = "G7",
                        "bb" = "G8",
                        "BB" = "WT"),
         D = fct_recode(Genotype,
                        "DD" = "G1",
                        "dd" = "G2",
                        "DD" = "G3",
                        "dd" = "G4",
                        "DD" = "G5",
                        "dd" = "G6",
                        "DD" = "G7",
                        "dd" = "G8",
                        "DD" = "WT")
  )

Long_data <- Wide_data_v2 %>%
  pivot_longer(cols = starts_with("Leaf_senescence")|starts_with("Peduncle_senescence"),
               names_to = c("Organ", "Date"),
               names_pattern = "(.*)\\.(.*)")
Long_data_v2 <- Long_data %>%
  pivot_wider(names_from = Organ, values_from = value)

Long_data_v3 <- Long_data_v2 %>%
  mutate(Date = ymd(Date, tz = "GMT"),
         Days_after_heading = days.after.sowing(Date, EE_date),
         Degree_days_after_heading = degree_day_interval(EE_date, Date)
         )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#WRITE
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/2023-08-02_clean_data", sep = ""))
write_csv(Wide_data_v2, file = paste(table_name, "wide_", lubridate::today(), ".csv", sep = ""))
write_csv(Long_data_v3, file = paste(table_name, "long_", lubridate::today(), ".csv", sep = ""))
