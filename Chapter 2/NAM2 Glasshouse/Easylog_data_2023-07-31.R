# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'Temperature data
#'Get something out of EasyLog glasshouse data
#'
#'31/07/2023
#'
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
library(tidyverse)
library(lubridate) #sort out the dates

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#OPEN
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-06-14 EasyLog NAC-5 greenhouse")
temp_2023 <- read_delim("Easylog_NAC5_2023-06-14.txt", delim = "," , skip =1,
                        col_names = c("Index", "Datetime", "Temp", "High_alarm", "Serial_num"))

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2022-10-06 Easylog NAM-2 greenhouse")
temp_2022_pt1 <- read_delim("Easylog_NAM2_part1_2022-10-06_2022-12-06.txt", delim = "," , skip =1,
                        col_names = c("Index", "Datetime", "Temp", "Serial_num"))
temp_2022_pt2 <- read_delim("Easylog_NAM2_part2_2022-12-06_2023-01-23.txt", delim = "," , skip =1,
                            col_names = c("Index", "Datetime", "Temp", "Serial_num"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#THEME

theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 24), plot.title = element_text(size = 24), legend.title = element_text(size=24),
                  axis.text = element_text(size = 20), legend.text = element_text(size = 20), strip.text = element_text(size = 20),
                  panel.grid = element_blank(), legend.position = "right")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
obs_per_day = 24*4

#Set daytime only 8am - 8pm (based on graph image)

daylight_hours <- function(temp_table){
  temp_table <- temp_table %>%
    select(Index, Datetime, Temp) %>%
    mutate(Hour = hour(Datetime), Temp = as.numeric(Temp)) %>%
    filter(Hour %in% 8:19)
}

temp_2023_v2 <- daylight_hours(temp_2023)
temp_2022_pt1_v2 <- daylight_hours(temp_2022_pt1)
temp_2022_pt2_v2 <- daylight_hours(temp_2022_pt2)
temp_2022_v2 <- rbind(temp_2022_pt1_v2, temp_2022_pt2_v2)

days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

Sowing_date <- lubridate::ymd('2023-02-15') #Sowing date the same for NAM2 data

temp_2023_summary <- temp_2023_v2 %>%
  mutate(Date = date(Datetime)) %>%
  group_by(Date) %>%
  summarise(Mean_day_temp = mean(Temp)) %>%
  mutate(DAS = days.after.sowing(Date, Sowing_date))

Sowing_date <- lubridate::ymd('2022-09-21') #Sowing date the same for NAM2 data

temp_2022_summary <- temp_2022_v2 %>%
  mutate(Date = date(Datetime)) %>%
  group_by(Date) %>%
  summarise(Mean_day_temp = mean(Temp)) %>%
  mutate(DAS = days.after.sowing(Date, Sowing_date))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#PLOT

ggplot(temp_2023_summary, aes(x = Date, y = Mean_day_temp)) +
  geom_line(col = "black") +
  geom_line(data = temp_2022_summary, col = "blue") +
  my_theme

ggplot(temp_2023_summary, aes(x = DAS, y = Mean_day_temp)) +
  geom_line(col = "black") +
  geom_line(data = temp_2022_summary, col = "blue") +
  labs(x = "Days after sowing", y = "Mean daytime temperature (oC)") +
  my_theme

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

night_hours <- function(temp_table){
  temp_table <- temp_table %>%
    select(Index, Datetime, Temp) %>%
    mutate(Hour = hour(Datetime), Temp = as.numeric(Temp)) %>%
    filter(Hour %in% c(0:7, 20:23))
}

temp_2023_v2 <- night_hours(temp_2023)
temp_2022_pt1_v2 <- night_hours(temp_2022_pt1)
temp_2022_pt2_v2 <- night_hours(temp_2022_pt2)
temp_2022_v2 <- rbind(temp_2022_pt1_v2, temp_2022_pt2_v2)

days.after.sowing <- function(test_date, sowing_date){
  return(as.numeric(interval(start = sowing_date, end = test_date), "days"))
} 

Sowing_date <- lubridate::ymd('2023-02-15') #Sowing date the same for NAM2 data

temp_2023_night_summary <- temp_2023_v2 %>%
  mutate(Date = date(Datetime)) %>%
  group_by(Date) %>%
  summarise(Mean_day_temp = mean(Temp)) %>%
  mutate(DAS = days.after.sowing(Date, Sowing_date))

Sowing_date <- lubridate::ymd('2022-09-21') #Sowing date the same for NAM2 data

temp_2022_night_summary <- temp_2022_v2 %>%
  mutate(Date = date(Datetime)) %>%
  group_by(Date) %>%
  summarise(Mean_day_temp = mean(Temp)) %>%
  mutate(DAS = days.after.sowing(Date, Sowing_date))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#PLOT

ggplot(temp_2023_night_summary, aes(x = Date, y = Mean_day_temp)) +
  geom_line(col = "blue", lty = 2) +
  geom_line(data = temp_2022_night_summary, col = "blue", lty = 2) +
  geom_line(data = temp_2022_summary, col = "red") +
  geom_line(data = temp_2023_summary, col = "red") +
  labs(x = "Date", y = "Mean temperature (oC)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
  my_theme
