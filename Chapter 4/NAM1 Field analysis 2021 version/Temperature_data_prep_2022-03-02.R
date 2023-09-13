# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#Temp data
#2022-03-02
#Get some daily mean temperatures for 2021 and 2022
#Use WOW data.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(lubridate)
require(readxl)
library(tidyverse)
library(lubridate) #sort out the dates

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN
setwd("U:/Field senescence re-analysis/weather_data")
Air_temperature_Andrewsfield_MayJuly2021 <- read_excel("Air_temperature_Andrewsfield_MayJuly2021.xlsx", 
                                                       col_types = c("date", "numeric"))
Air_temperature_Andrewsfield_MayJuly2022 <- read_excel("Air_temperature_Andrewsfield_MayJuly2022.xlsx", 
                                                       col_types = c("date", "numeric"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#AVERAGE
temperature_by_day <- Air_temperature_Andrewsfield_MayJuly2022 %>%
  mutate(date = date(Obs_time)) %>%
  group_by(date) %>%
  summarise(mean_temp=mean(Air_temp, na.rm = TRUE))

#Calculate thermal time
#NB this is problematic in January with negative numbers but thankfully I don't have to worry
#about this in July
temperature_by_day <- temperature_by_day %>%
  arrange(date) %>%
  mutate(cumulative_temp = cumsum(mean_temp))

cumulative_temp_1may <- temperature_by_day %>%
  filter(date == "2022-05-01") %>%
  select(cumulative_temp) %>%
  unlist()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#THEME

theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 24), plot.title = element_text(size = 24), legend.title = element_text(size=24),
                  axis.text = element_text(size = 20), legend.text = element_text(size = 20), strip.text = element_text(size = 20),
                  panel.grid = element_blank(), legend.position = "right")
#Choose interval

# 3 months
focus_interval <- interval(start = ymd("2022-05-20"), end = ymd("2022-07-31"), tzone = "GMT")

#July
focus_interval <- interval(start = ymd("2022-07-01"), end = ymd("2022-07-31"), tzone = "GMT")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#PLOT
temp_plot <- temperature_by_day %>%
  filter(date %within% focus_interval) %>%
  ggplot(aes(x=date, y=mean_temp)) +
  geom_line()

temp_plot + labs(x = "Date in 2022", y = "Mean daily air temperature (oC)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %d") +
  my_theme

temperature_by_day %>%
  filter(date %within% focus_interval) %>%
  ggplot(aes(x=date, y=cumulative_temp - cumulative_temp_1may)) +
  geom_line() +
  geom_point(size = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#WRITE
write_csv(temperature_by_day, file = "daily_mean_temperature_2022.csv")


