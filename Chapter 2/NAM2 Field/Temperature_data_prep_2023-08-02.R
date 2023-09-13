# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#Temp data
#2022-08-02
#Get some daily mean temperatures from Church Farm weather station data
#'
#'Capital letters for column names
#'# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'Metadata
#'Daily
#'date - date or datetime
# AirTC_Avg / mean_temp - average air temperature of the datalogger in °C
#'
#'
#'Will take average of hourly AirTC_Avg to calculate daily temperature
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#LOAD
library(tidyverse)
library(lubridate) #sort out the dates

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN
Average_air_temp_2023_05_to_2023_07 <- read_excel("Weather_data/Average air temp 2023-05 to 2023-07.xlsx", 
                                                       col_types = c("date", "guess", "numeric"))

#AVERAGE
#already a daily average
#sort date formatting
temperature_by_day <- Average_air_temp_2023_05_to_2023_07 %>%
  mutate(Date = as_datetime(ifelse(is.na(Date), dmy_hm(Date_time), Date))) %>%
  select(Date, Mean_temp)


#Calculate thermal time
#NB this is problematic in January with negative numbers but thankfully I don't have to worry
#about this in July
temperature_by_day <- temperature_by_day %>%
  arrange(Date) %>%
  mutate(Cumulative_temp = cumsum(Mean_temp))

cumulative_temp_1may <- temperature_by_day %>%
  filter(Date == ymd("2023-05-01")) %>%
  select(Cumulative_temp) %>%
  unlist()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#THEME

theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 24), plot.title = element_text(size = 24), legend.title = element_text(size=24),
                  axis.text = element_text(size = 20), legend.text = element_text(size = 20), strip.text = element_text(size = 20),
                  panel.grid = element_blank(), legend.position = "right")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#PLOT

#Choose interval

# 3 months
focus_interval <- interval(start = ymd("2023-05-01"), end = ymd("2023-07-31"), tzone = "GMT")

#July
# focus_interval <- interval(start = ymd("2023-07-01"), end = ymd("2023-07-31"), tzone = "GMT")


temp_plot <- temperature_by_day %>%
  filter(Date %within% focus_interval) %>%
  ggplot(aes(x=as_date(Date), y=Mean_temp)) +
  geom_line()

#THIS PLOT
temp_plot + labs(x = "Date", y = "Mean daily air temperature (oC)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
  my_theme

temperature_by_day %>%
  filter(Date %within% focus_interval) %>%
  ggplot(aes(x=Date, y=Cumulative_temp - cumulative_temp_1may)) +
  geom_line() +
  geom_point(size = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Calculate

# Before June 10
focus_interval <- interval(start = ymd("2023-05-01"), end = ymd("2023-06-09"), tzone = "GMT")
# After June 10
focus_interval <- interval(start = ymd("2023-06-10"), end = ymd("2023-07-31"), tzone = "GMT")

temperature_by_day %>%
  filter(Date %within% focus_interval) %>%
  summarise(mean(Mean_temp))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#WRITE
write_csv(temperature_by_day, file = "daily_mean_temperature_2023.csv")
