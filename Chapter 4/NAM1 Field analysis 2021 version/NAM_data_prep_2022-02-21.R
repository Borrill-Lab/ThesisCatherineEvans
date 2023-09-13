# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'21/02/2022
#'Data preparation
#'
#'I need to
# 1.	Extract the raw data from the spreadsheets
# 3.	Combine together
# 4.	Make a “wide format” table with 1 plot per row and a “long format” table with 1 plot * 1 timepoint per row.
#'
#'11/03/2022
#'5. Calculate time in "days since ear emergence" and "thermal time since ear emergence"
#'
#'23/08/2023
#'New temperature data
#'
#'#'Built with R 4.0.4 and tidyverse
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(tidyverse)
require(reshape2)
require(readxl)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN file
setwd("U:/RAGT Field Senescence 2021/Data")

table_name <- "Field_results_NAM_"
Field_results_neat_2021 <- read_excel("Field_results_neat_2021.xlsx", 
                                      sheet = "All results NAM NILs", skip = 5)

GridScore_data <- read_delim("data-little-heath-nam-nils-2021-07-28.txt", delim = "\t")
View(GridScore_data)

#GridScore_data is not useable in its current state because the 4 digit trial codes are NOT unique - they're replicated!

#Temperature data from Andrewsfield
setwd("U:/Field senescence re-analysis")
temperature_by_day <- read_csv("weather_data/daily_mean_temperature_2021.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY DATA
summary(Field_results_neat_2021)

first_of_may = ymd("2021-05-01", tz = "GMT")

Wide_data <- Field_results_neat_2021 %>%
  mutate(
    Protein = na_if(Protein, "not done"),
    #let's make those into NAs
    EE = na_if(EE, "not done"),
    #let's make those into NAs
    `Ear_senescence.2021-07-02` = 0
    #there was no ear senescence at this time point and it didn't get scored
  )


#Ear emergence
#Infer that Rep 2 will have the same EE dates as Rep 1
Wide_data <- Wide_data %>%
  mutate(EE_inferred = EE) %>%
  arrange(Trial_code, Rep) %>%
  tidyr::fill(EE_inferred, .direction = "down")

#and then put back the ones which are supposed to stay NA and shouldn't have been filled
#yes this is a fudge
Wide_data$EE_inferred[which(Wide_data$Line_name %in% c("SKYFALL", "RGT BOWEN"))] <- NA

Wide_data <- Wide_data %>%
  mutate(
EE_date = first_of_may + days(EE_inferred %/% 1) + hours((EE_inferred %% 1) * 24)
)

#Convert EE into a date in such a way as to preserve the half-days.
Wide_data$EE_date


#LONG DATA for plotting the timecourse
Long_data <- Wide_data %>%
  pivot_longer(cols = starts_with("Leaf_senescence")|starts_with("Ear_senescence"),
               names_to = c("Organ", "Date"),
               names_pattern = "(.*)\\.(.*)")
Long_data_v2 <- Long_data %>%
  pivot_wider(names_from = Organ, values_from = value)

Long_data_v3 <- Long_data_v2 %>%
  mutate(Date = ymd(Date, tz = "GMT"),
        Days_after_heading = Date - EE_date,
         Degree_days_after_heading = degree_day_interval(EE_date, Date))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'Calculates the difference between 2 dates in degree days
#'EXCLUDES the temperature on the start date
#'So that if you put in the same date for start and end, the answer is 0
#'Half days are currently rounded UP
degree_day_interval <- function(start_date, end_date){
  require(lubridate); require(tidyverse)
  if(length(start_date)==length(end_date)){
    degree_days = c()
    for(i in 1:length(start_date)){
      focus_interval <- interval(start = start_date[i]+days(1), end = end_date[i], tzone = "GMT")
      degree_days[i] <- temperature_by_day %>%
        filter(date %within% focus_interval) %>%
        summarise(sum(mean_temp)) %>%
        unlist()
    }
    return(degree_days)
  }else{
    print("Error: Start date and end date need to be vectors of the same length")
  }
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#WRITE
setwd("U:/RAGT Field Senescence 2021/Data")
write_csv(Wide_data, file = paste(table_name, "wide_", lubridate::today(), ".csv", sep = ""))
write_csv(Long_data_v3, file = paste(table_name, "long_", lubridate::today(), ".csv", sep = ""))
