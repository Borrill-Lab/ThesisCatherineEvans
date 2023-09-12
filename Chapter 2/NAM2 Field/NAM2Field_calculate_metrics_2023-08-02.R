# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Calculate senescence metrics
#'RUN AFTER NAM2Field_data_prep_2023-08-02
#'
#'16/03/2022
#'Linear interpolation method
#'The true progress of senescence may follow a curve
#'However there are probably too few data points to fit a convincing curve
#'Therefore we will use linear interpolation
#'A tried and tested method since the time of the Greeks
#'Using the r function stats::approx
#'
#'17/03/2022
#'Use Days_after_heading and Degree_days_after_heading instead of days_after_1may
#'Make new file and cut the clutter
#'Get rid of Mean_leaf as AUC is better
#'
#'02/08/2023
#'Adapt for NAM-2 phenotyping data 2023
#'
#'03/08/2023
#'With Degree_days_after_heading
#'
#'04/08/2023
#'With Dosage and A.B.D
#'
#'Run from U: drive
#'Built in R 4.0.4 and tidyverse (mainly readr, dplyr and ggplot2)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(tidyverse)
require(lubridate)
require(MESS) #for auc()
require(ggpubr) #for ggarrange()
theme_set(theme_bw(base_size = 8))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN
get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/2023-08-02_clean_data", sep = ""))

Long_data <- read_csv("Field_results_NAM2_long_2023-08-04.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

Wide_data <- read_csv("Field_results_NAM2_wide_2023-08-04.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#LABELS
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
table_name <- "Field_results_NAM2_"

#Mean senescence values are dataset-specific as they depend on how many data are collected.
#Therefore if I want to compare between sowing dates I'll need to normalise to a range in days after sowing or heading

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#RUN

Summary_Days_after_heading <- Long_data %>%
  filter(!is.na(Leaf_senescence) & Days_after_heading > 0) %>%
  arrange(Days_after_heading) %>%
  group_by(Index) %>%
  summarise(
    Leaf_TT10 = approx(Leaf_senescence, Days_after_heading, c(10), ties = "ordered")$y,
    Leaf_TT30 = approx(Leaf_senescence, Days_after_heading, c(30), ties = "ordered")$y,
    Leaf_TT50 = approx(Leaf_senescence, Days_after_heading, c(50), ties = "ordered")$y,
    Leaf_TT70 = approx(Leaf_senescence, Days_after_heading, c(70), ties = "ordered")$y,
    Leaf_TT90 = approx(Leaf_senescence, Days_after_heading, c(90), ties = "ordered")$y,
    Dur_Leaf_Sen = ifelse(is.na(Leaf_TT10), Leaf_TT90 - min(Days_after_heading),
                          Leaf_TT90 - Leaf_TT10),
    Rate_Leaf_Sen = 80/Dur_Leaf_Sen,
    AUC_Leaf_Sen = auc(Days_after_heading, Leaf_senescence, type = "linear")
  )

Summary_Degree_days_after_heading <- Long_data %>%
  filter(!is.na(Leaf_senescence) & Degree_days_after_heading > 0) %>%
  arrange(Degree_days_after_heading) %>%
  group_by(Index) %>%
  summarise(
    Leaf_TT10 = approx(Leaf_senescence, Degree_days_after_heading, c(10), ties = "ordered")$y,
    Leaf_TT30 = approx(Leaf_senescence, Degree_days_after_heading, c(30), ties = "ordered")$y,
    Leaf_TT50 = approx(Leaf_senescence, Degree_days_after_heading, c(50), ties = "ordered")$y,
    Leaf_TT70 = approx(Leaf_senescence, Degree_days_after_heading, c(70), ties = "ordered")$y,
    Leaf_TT90 = approx(Leaf_senescence, Degree_days_after_heading, c(90), ties = "ordered")$y,
    Dur_Leaf_Sen = ifelse(is.na(Leaf_TT10), Leaf_TT90 - min(Degree_days_after_heading),
                          Leaf_TT90 - Leaf_TT10),
    Rate_Leaf_Sen = 80/Dur_Leaf_Sen,
    AUC_Leaf_Sen = auc(Days_after_heading, Leaf_senescence, type = "linear")
  )

#LOOK
#What do data look like
summary(Summary_Days_after_heading)

key_plots <- list(
  ggplot(Summary_Days_after_heading, aes(x = Dur_Leaf_Sen)) +
    geom_histogram(binwidth = 5),
  ggplot(Summary_Days_after_heading, aes(x = Rate_Leaf_Sen)) +
    geom_histogram(binwidth = 2.5),
  ggplot(Summary_Days_after_heading, aes(x = AUC_Leaf_Sen)) +
    geom_histogram(binwidth = 250),
  ggplot(Summary_Days_after_heading, aes(x = AUC_Leaf_Sen, y = Rate_Leaf_Sen)) +
    geom_point()
)

arrange_2 <- ggarrange(key_plots[[1]], key_plots[[2]], key_plots[[3]], key_plots[[4]])
arrange_2

summary(Summary_Degree_days_after_heading)

key_plots <- list(
  ggplot(Summary_Degree_days_after_heading, aes(x = Dur_Leaf_Sen)) +
    geom_histogram(binwidth = 50),
  ggplot(Summary_Degree_days_after_heading, aes(x = Rate_Leaf_Sen)) +
    geom_histogram(binwidth = 0.025),
  ggplot(Summary_Degree_days_after_heading, aes(x = AUC_Leaf_Sen)) +
    geom_histogram(binwidth = 100),
  ggplot(Summary_Degree_days_after_heading, aes(x = AUC_Leaf_Sen, y = Rate_Leaf_Sen)) +
    geom_point()
)

arrange_3 <- ggarrange(key_plots[[1]], key_plots[[2]], key_plots[[3]], key_plots[[4]])

arrange_3


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#COMBINE
#Everything in 1 table :-)
Wide_data_Days_after_heading <- inner_join(Wide_data, Summary_Days_after_heading, by = "Index")
Wide_data_Degree_days_after_heading <- inner_join(Wide_data, Summary_Degree_days_after_heading, by = "Index")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#SAVE
#Save those demo figures
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/Figures", sep = ""))
ggsave(filename = paste("NAM_statistics_DAH_Leaf_",today(),".pdf",sep=""), arrange_2, width = 240, height = 180, units = "mm", device = "pdf")
ggsave(filename = paste("NAM_statistics_DDAH_Leaf_",today(),".pdf",sep=""), arrange_3, width = 240, height = 180, units = "mm", device = "pdf")

#Save the data tables with calculated metrics
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/2023-08-02_clean_data", sep = ""))
write_csv(Wide_data_Days_after_heading, file = paste(table_name, "wide_DAH_", lubridate::today(), ".csv", sep = ""))
write_csv(Wide_data_Degree_days_after_heading, file = paste(table_name, "wide_DDAH_", lubridate::today(), ".csv", sep = ""))
