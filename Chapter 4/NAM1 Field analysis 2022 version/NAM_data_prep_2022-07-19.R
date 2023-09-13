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
#'Andrewsfield temperature data
#'
#'#'Built with R 4.0.4 and tidyverse
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(tidyverse)
require(readxl)
require(lubridate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN file
setwd("U:/Field Senescence 2022")


table_name <- "Field_results_NAM_"

Field_results_neat_2022 <- read_excel("Raw_data/Catherine_Gpc_NILs_senescence_neat_2022-07-14.xlsx", skip = 5)

Last_year_data <- read_csv("U:/RAGT Field Senescence 2021/Data/Field_results_NAM_wide_2022-03-11.csv",
                           col_types = cols("NAMA1"= col_factor(),
                                            "NAMB1"= col_factor(),
                                            "Rep"= col_factor()))

Height_Ear_emergence_2022_anonymised <- read_excel("Raw_data/Height_Ear_emergence_2022_anonymised.xlsx", 
                                                        sheet = "Gpc NILs Project Block")

#Temperature data from Andrewsfield
setwd("U:/Field senescence re-analysis")
temperature_by_day <- read_csv("weather_data/daily_mean_temperature_2022.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY DATA
summary(Field_results_neat_2022)

first_of_may = ymd("2022-05-01", tz = "GMT")

Wide_data <- Field_results_neat_2022 %>%
  mutate(
    `Leaf_senescence.2022-07-15` = 100,
    `Ear_senescence.2022-06-16` = 0,
    `Ear_senescence.2022-06-22` = 0,
    `Ear_senescence.2022-07-15` = 100
  )

#Add NAM-A1 and NAM-B1 genotypes
#Add NAMA1 and NAMB1 genotypes
Genotype_key <- Last_year_data %>%
  select(Line_name, NAMA1, NAMB1) %>%
  unique()

Wide_data_v2 <- Wide_data %>%
  extract(Line_name, "Line_name_2", "(.*)\\.", remove = FALSE) %>%
  mutate(Line_name = ifelse(is.na(Line_name_2), Line_name, Line_name_2)) %>%
  left_join(Genotype_key, by = c("Line_name"))

#Add Ear emergence etc.
EE_data <- Height_Ear_emergence_2022_anonymised %>%
  select(`Data identifier (DATA_ID)`, `Hab 1`, EE_REL, PLH) %>%
  rename(DATA_ID = "Data identifier (DATA_ID)", Habit = "Hab 1", EE_REL = "EE_REL", Height = "PLH")

#Assuming EE_REL is EE relative to 1st May, as with previous years' data
Wide_data_v2 <- Wide_data_v2 %>%
  left_join(EE_data, by = "DATA_ID") %>%
  mutate(
    EE_date = first_of_may + days(EE_REL %/% 1) + hours((EE_REL %% 1) * 24)
  )



#LONG DATA for plotting the timecourse
Long_data <- Wide_data_v2 %>%
  pivot_longer(cols = starts_with("Leaf_senescence")|starts_with("Ear_senescence"),
               names_to = c("Organ", "Date"),
               names_pattern = "(.*)\\.(.*)")
Long_data_v2 <- Long_data %>%
  pivot_wider(names_from = Organ, values_from = value)

# Long_data_v3 <- Long_data_v2 %>%
#   mutate(Date = ymd(Date, tz = "GMT"),
#          days_after_1may = Date - first_of_may)

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
setwd("U:/Field senescence re-analysis/Processed_data")
write_csv(Wide_data_v2, file = paste(table_name, "2022_", "wide_", lubridate::today(), ".csv", sep = ""))
write_csv(Long_data_v3, file = paste(table_name, "2022_", "long_", lubridate::today(), ".csv", sep = ""))
