# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Calculate senescence metrics
#'16/03/2022
#'Linear interpolation method
#'The true progress of senescence may follow a curve
#'However there are probably too few data points to fit a convincing curve
#'Therefore we will use linear interpolation
#'A tried and tested method since the time of the Greeks
#'Using the r function stats::approx
#'
#'19/07/2022
#'Run with 2022 NAM-1 NILs data
#'USE DAYS AFTER 1 MAY
#'Remove stats demo
#'
#'23/08/2023
#'USE Degree days after heading
#'
#'Run from RAGT internal
#'Built in R 4.0.4 and tidyverse (mainly readr, dplyr and ggplot2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(tidyverse)
require(lubridate)
require(MESS)
require(ggpubr)
theme_set(theme_bw(base_size = 8))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN
setwd("U:/Field senescence re-analysis/Processed_data")

Long_data <- read_csv("Field_results_NAM_2022_long_2023-08-23.csv")
Wide_data <- read_csv("Field_results_NAM_2022_wide_2023-08-23.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
#LABELS
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))

table_name <- "Field_results_NAM_"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY
Long_data <- Long_data %>%
  mutate(Year = 2022)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#AUC NORMALISATION TEST

#Mean senescence values are dataset-specific as they depend on how many data are collected.
#Therefore if I want to compare between sowing dates I'll need to normalise to a range in days after sowing or heading

#normalise to a range in days after sowing or heading
Long_data %>% group_by(Date) %>%
  summarise(max(Degree_days_after_heading),
            min(Degree_days_after_heading),
            mean(Degree_days_after_heading)
  )
#End AUC at minimum scoring endpoint - ensures interpolation only, no extrapolation.
#Start point is less important as additional score points are 0 so do not add to total auc
auc_limit <- Long_data %>% filter(Date == max(Date)) %>%
  summarise(
    min_DDAH = floor(min(Degree_days_after_heading))
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#RUN

Long_data <- Long_data %>%
  mutate(Control = ifelse(Line_name %in% c("SKYFALL", "RGT BOWEN"), TRUE, FALSE),
         Genotype = interaction(NAMA1, NAMB1))

Summary_Degree_days_after_heading <- Long_data %>%
  group_by(DATA_ID) %>%
  summarise(
    Leaf_TT10 = approx(Leaf_senescence, Degree_days_after_heading, c(10), ties = "ordered")$y,
    Leaf_TT30 = approx(Leaf_senescence, Degree_days_after_heading, c(30), ties = "ordered")$y,
    Leaf_TT50 = approx(Leaf_senescence, Degree_days_after_heading, c(50), ties = "ordered")$y,
    Leaf_TT70 = approx(Leaf_senescence, Degree_days_after_heading, c(70), ties = "ordered")$y,
    Leaf_TT90 = approx(Leaf_senescence, Degree_days_after_heading, c(90), ties = "ordered")$y,
    Dur_Leaf_Sen = ifelse(is.na(Leaf_TT10), Leaf_TT90 - min(Degree_days_after_heading),
                          Leaf_TT90 - Leaf_TT10),
    Rate_Leaf_Sen = 80/Dur_Leaf_Sen,
    AUC_Leaf_Sen = auc(Degree_days_after_heading, Leaf_senescence, type = "linear", to = auc_limit$min_DDAH),
    Ear_TT10 = approx(Ear_senescence, Degree_days_after_heading, c(10), ties = "ordered")$y,
    Ear_TT30 = approx(Ear_senescence, Degree_days_after_heading, c(30), ties = "ordered")$y,
    Ear_TT50 = approx(Ear_senescence, Degree_days_after_heading, c(50), ties = "ordered")$y,
    Ear_TT70 = approx(Ear_senescence, Degree_days_after_heading, c(70), ties = "ordered")$y,
    Ear_TT90 = approx(Ear_senescence, Degree_days_after_heading, c(90), ties = "ordered")$y,
    Dur_Ear_Sen = ifelse(is.na(Ear_TT10), Ear_TT90 - min(Degree_days_after_heading),
                          Ear_TT90 - Ear_TT10),
    Rate_Ear_Sen = 80/Dur_Ear_Sen,
    AUC_Ear_Sen = auc(Degree_days_after_heading, Ear_senescence, type = "linear", to = auc_limit$min_DDAH)
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#LOOK
#What do data look like
#Days
summary(Summary_Degree_days_after_heading)

key_plots <- list(
  ggplot(Summary_Degree_days_after_heading, aes(x = Dur_Leaf_Sen)) +
    geom_histogram(binwidth = 1),
  ggplot(Summary_Degree_days_after_heading, aes(x = Rate_Leaf_Sen)) +
    geom_histogram(binwidth = 1),
  ggplot(Summary_Degree_days_after_heading, aes(x = AUC_Leaf_Sen)) +
    geom_histogram(binwidth = 50),
  ggplot(Summary_Degree_days_after_heading, aes(x = AUC_Leaf_Sen, y = Dur_Leaf_Sen)) +
    geom_point()
)

key_plots[[1]]
key_plots[[2]]
key_plots[[3]]
key_plots[[4]]

arrange_2 <- ggarrange(key_plots[[1]], key_plots[[2]], key_plots[[3]], key_plots[[4]])

Wide_data_Degree_days_after_heading <- inner_join(Wide_data, Summary_Degree_days_after_heading, by = c("DATA_ID"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#SAVE
#Save the data tables with calculated metrics
# setwd("../NAM_Figures")
# ggsave(filename = paste("NAM_statistics_1may_Leaf_",today(),".png",sep=""), arrange_2, width = 240, height = 180, units = "mm", device = "png")
# ggsave(filename = paste("NAM_statistics_1may_Leaf_",today(),".pdf",sep=""), arrange_2, width = 240, height = 180, units = "mm", device = "pdf")


setwd("../Processed_data")
# write_csv(Wide_data_Days_after_heading, file = paste(table_name, "wide_DAH_", lubridate::today(), ".csv", sep = ""))
write_csv(Wide_data_Degree_days_after_heading, file = paste(table_name, "2022_", "wide_DDAH_", lubridate::today(), ".csv", sep = ""))
