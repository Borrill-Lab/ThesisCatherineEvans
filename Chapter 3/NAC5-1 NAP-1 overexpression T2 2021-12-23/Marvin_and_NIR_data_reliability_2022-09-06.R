#'Checking reliability of grain mass and NIR data
#'
#'Also see NL_exploratory_2022-07-01.R
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(dplyr)
require(readr)
require(ggplot2)
require(lubridate)
require(tidyr)
theme_set(theme_bw())
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER")

NAC5_NAP_Overexpression_Marvin_clean <- read_csv("2022-08-12_clean_data/NAC5_NAP_Overexpression_Marvin_clean.csv")

Kronos_Marvin_clean <- read_csv("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data/Kronos_Marvin_clean.csv")

Kronos_NIR_Extensive_Report_Wheat_clean <- read_csv("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data/Kronos_NIR_Extensive Report_Wheat_clean.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#ANALYSE MARVIN DATA

Marvin_all <- rbind(Kronos_Marvin_clean, NAC5_NAP_Overexpression_Marvin_clean) %>%
  mutate(Date_Marvin = dmy(Date_Marvin))

#Consistent names for these as calibrators
Marvin_all[16,1] <- "CAL_49"
Marvin_all[17,1] <- "CAL_50"
Marvin_all[19,1] <- "CAL_85"

calibrators <- c("CAL_49", "CAL_50", "CAL_85")

#Which seeds were measured more than once?
Counts <- Marvin_all %>%
  group_by(ID) %>%
  count() %>%
  filter(n>1)

Marvin_calibrators <- Marvin_all %>%
  filter(ID %in% calibrators & Seed_num != 347)
#Remove one outlier taken when some seeds had been removed from the CAL_50 packet

Marvin_pairs <- Marvin_all %>%
  filter(ID %in% Counts$ID & !(ID %in% calibrators))

Kronos_Marvin_pared_back <- Kronos_Marvin_clean %>%
  filter(!(ID %in% Counts$ID))

#Calibrator samples; not really sure about these.
Marvin_calibrators %>%
  ggplot(aes(x = Date_Marvin, y = Seed_num, col = ID)) +
  geom_point()

Marvin_calibrators %>%
  ggplot(aes(x = Date_Marvin, y = Seed_num, col = Machine)) +
  geom_point() +
  facet_wrap(vars(ID), scales = "free_y")

Marvin_calibrators %>%
  ggplot(aes(x = ID, y = Seed_num)) +
  geom_point()

variable_list <- c("Seed_num",    "Weight",      "TGW",         "O_area",      "O_width",     "Min_width",  
"Max_width",   "O_length",    "Min_length",  "Max_length")

for(variable in variable_list){
  plot <- Marvin_calibrators %>%
    ggplot(aes_string(x = "Date_Marvin", y = variable, col = "Machine")) +
    geom_point() +
    facet_wrap(vars(ID), scales = "free_y")
  print(plot)
}

#Within-machine comparison
Marvin_calibrators %>%
  group_by(ID, Machine) %>%
  summarise(across(Seed_num:Max_length, mean))

Marvin_calibrators %>%
  group_by(ID, Machine) %>%
  summarise(across(Seed_num:Max_length, sd))

#ALL
Marvin_calibrators %>%
  group_by(ID) %>%
  summarise(across(Seed_num:Max_length, mean))

Marvin_calibrators %>%
  group_by(ID) %>%
  summarise(across(Seed_num:Max_length, sd))

#PAIRS
plots <- list(0,0,0,0,0,0,0,0,0,0)
plots_machine <- list(0,0,0,0,0,0,0,0,0,0)

for(i in 1:length(variable_list)){
  Pairs_Date <- Marvin_pairs %>%
    filter(Machine == "Marvin_2") %>%
    pivot_wider(id_cols = ID, names_from = Date_Marvin, values_from = all_of(variable_list[i])) %>%
    filter(!is.na(`2022-04-27`))
  
  plot <- ggplot(Pairs_Date, aes(x = `2022-04-27`, y = `2022-08-16`)) +
    geom_point(alpha = 0.2, size = 4) +
      ggpubr::stat_cor() +
    ggtitle(variable_list[i]) +
    geom_abline(slope = 1, intercept = 0)
  
  plots[[i]] <- plot
  
  Pairs_Machine <- Marvin_pairs %>%
    pivot_wider(id_cols = ID, names_from = Machine, values_from = all_of(variable_list[i]), values_fn = mean) %>%
    filter(!is.na(`Marvin_1`))
  
  plot <- ggplot(Pairs_Machine, aes(x = Marvin_1, y = Marvin_2)) +
    geom_point(alpha = 0.2, size = 4) +
    ggpubr::stat_cor() +
    ggtitle(variable_list[i]) +
    geom_abline(slope = 1, intercept = 0)
  
  plots_machine[[i]] <- plot
  
  print(variable_list[i])
  print(mean(Pairs_Machine$Marvin_1 - Pairs_Machine$Marvin_2))
  print(sd(Pairs_Machine$Marvin_1 - Pairs_Machine$Marvin_2))
}

ggpubr::ggarrange(plotlist = plots)

ggpubr::ggarrange(plotlist = plots_machine)

#Oh dear, Marvin_1 and Marvin_2 are calibrated differently; Marvin_1 consistently claims seeds are larger

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#ADJUST
#Make a new dataset appropriate for analysis
#Adjust area, length & width values by average/systematic machine difference to match Marvin 2
#Choose Marvin 2 data where there are two
#Choose the most recent Marvin value where there are two

#KRONOS
#Easy as already all Marvin 2
#Choose the most recent Marvin value where there are two
Kronos_Marvin_sorted <- Kronos_Marvin_clean %>%
  mutate(Date_Marvin = lubridate::dmy(Date_Marvin)) %>%
  arrange(desc(Date_Marvin))

Kronos_pairs <- Kronos_Marvin_sorted[duplicated(Kronos_Marvin_sorted$ID),]
Kronos_Marvin_unique <- Kronos_Marvin_sorted[!duplicated(Kronos_Marvin_sorted$ID),]

#OVEREXPRESSION
Overexpression_Marvin_sorted <- NAC5_NAP_Overexpression_Marvin_clean %>%
  mutate(Date_Marvin = lubridate::dmy(Date_Marvin)) %>%
  arrange(desc(Machine), desc(Date_Marvin))

Overexpression_pairs <- Overexpression_Marvin_sorted[duplicated(Overexpression_Marvin_sorted$ID),]
Overexpression_Marvin_unique <- Overexpression_Marvin_sorted[!duplicated(Overexpression_Marvin_sorted$ID),] 

#ADJUST
#numbers from "pairs" above
area_1vs2 = 2.90072
width_1vs2 = 0.3235283
length_1vs2 = 0.3409747

Overexpression_unique_adjusted <- Overexpression_Marvin_unique %>%
  mutate(O_area_adj = ifelse(Machine == "Marvin_1", O_area - area_1vs2, O_area),
         O_width_adj = ifelse(Machine == "Marvin_1", O_width - width_1vs2, O_width),
         O_length_adj = ifelse(Machine == "Marvin_1", O_length - length_1vs2, O_length)
  )

setwd("./2022-08-12_clean_data/")
write_csv(Kronos_Marvin_unique, file = paste('Kronos_Marvin_unique_', today(), '.csv', sep = ""))
write_csv(Overexpression_unique_adjusted, file = paste('Overexpression_Marvin_unique_', today(), '.csv', sep = ""))
setwd("..")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#ANALYSE NIR DATA

NIR_data <- Kronos_NIR_Extensive_Report_Wheat_clean

#Which seeds were measured more than once?
Counts <- NIR_data %>%
  group_by(`Sample ID`) %>%
  count() %>%
  filter(n>1)

NIR_calibrator <- NIR_data %>%
  filter(`Sample ID` == "CAL_50")

variable_list <- c("`Predicted Moisture %`",          "`Predicted Protein Dry basis %`",
                   "`Predicted Starch Dry basis %`",  "`Predicted NDF Dry basis %`",     "`Predicted Dry gluten As is %`",
                   "`Predicted Hardness -`")

for(variable in variable_list){
  plot <- NIR_calibrator %>%
    ggplot(aes_string(x = "Date_NIR", y = variable, col = "`Sample Tray`")) +
    geom_jitter(width = 0.1)
  print(plot)
}

NIR_pared_back <- NIR_data %>%
  filter(!(`Sample ID` %in% Counts$`Sample ID`))

#Join datas

Kronos_combined <- inner_join(Kronos_Marvin_pared_back, NIR_pared_back, by = c("ID" = "Sample ID"))

#What if moisture content affects weight measurements?
#Calculate dry weights
Kronos_combined <- Kronos_combined %>%
  mutate(Weight_dry = Weight - (Weight * `Predicted Moisture %` /100),
         TGW_dry = TGW - (TGW * `Predicted Moisture %` /100))

ggplot(Kronos_combined, aes(x = `Predicted Moisture %`, y = Weight)) +
  geom_point() +
  ggpubr::stat_cor()
ggplot(Kronos_combined, aes(x = `Predicted Moisture %`, y = TGW)) +
  geom_point() +
  ggpubr::stat_cor()

#Does it change the values much? No.
ggplot(Kronos_combined, aes(x = Weight, y = Weight_dry)) +
  geom_point() +
  ggpubr::stat_cor()

ggplot(Kronos_combined, aes(x = TGW, y = TGW_dry)) +
  geom_point() +
  ggpubr::stat_cor()

#Does it improve the correlation with seed number? Yes, although only very subtly
ggplot(Kronos_combined, aes(x = Seed_num, y = TGW)) +
  geom_point() +
  ggpubr::stat_cor()

ggplot(Kronos_combined, aes(x = Seed_num, y = TGW_dry)) +
  geom_point() +
  ggpubr::stat_cor()

ggplot(Kronos_combined, aes(x = Seed_num, y = Weight)) +
  geom_point() +
  ggpubr::stat_cor()

ggplot(Kronos_combined, aes(x = Seed_num, y = Weight_dry)) +
  geom_point() +
  ggpubr::stat_cor()

# Moisture content varies little in this experiment so has negligible effect on TGW

variable_list <- c("Seed_num", "Weight", "TGW",
  "`Predicted Moisture %`",          "`Predicted Protein Dry basis %`",
                   "`Predicted Starch Dry basis %`",  "`Predicted NDF Dry basis %`",     "`Predicted Dry gluten As is %`",
                   "`Predicted Hardness -`")
#Poor coverage
Kronos_combined <- Kronos_combined %>%
mutate(Poor_coverage = grepl("p", Comment))

for(variable in variable_list){
  plot <- Kronos_combined %>%
    ggplot(aes_string(x = "Poor_coverage", y = variable, col = "`Sample Tray`")) +
    geom_jitter(width = 0.1)
  print(plot)
}

ggplot(Kronos_combined, aes(x = Weight, y = `Predicted Protein Dry basis %`, col = `Poor_coverage`)) +
  geom_point() +
  geom_vline(xintercept = 16)
