#'Kronos poster graphs SOURCE
#'
#'26/09/2022
#'Build graphs for Kronos data for Crop Gen Seminar on 03/10/2022
#'Based on Kronos_exploratory_2022-09-07
#'
#'06/10/2022
#'Build graphs for JIC ASM poster
#'New themes.
#'
#'09/03/2023
#'Added local_file_path
#'
#'This SOURCE script loads packages and data.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(car) #Anova for non-orthogonal data
require(emmeans) #marginal means estimates for when p-values are not enough
require(arm)
require(lubridate)
require(tidyverse) #Load tidyverse last to avoid masking
require(ggpubr) #for ggarrange
theme_set(theme_bw())

local_file_path <- stringr::str_extract(getwd(), ".*OneDrive?-?Norwich?BioScience?Institutes/")
setwd(paste(local_file_path, "NAC transgenics/2021-03 NAC phenotyping", sep = ""))
source("Functions_for_phenotype_data.R")

setwd(paste(local_file_path, "NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/2022-03-07_clean_data", sep = ""))

Kronos_SPAD_yday <- read_csv(paste("Kronos_SPAD_yday_", gene, "_", date, ".csv", sep= ""))
Kronos_phenotyping_yday <- read_csv(paste("Kronos_phenotyping_yday_", gene, "_", date, ".csv", sep= ""))

#THEME
theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 24), plot.title = element_text(size = 24), legend.title = element_text(size=24),
                  axis.text = element_text(size = 20), legend.text = element_text(size = 20), strip.text = element_text(size = 20),
                  panel.grid = element_blank(), legend.position = "right")

#LABELS
genotype_lab <- paste(gene, "Genotype")
leaf_senescence_lab <- "Days to 25% leaf senescence"
peduncle_senescence_lab <- "Days to 100% peduncle senescence"
SPAD_lab <- "Flag leaf chlorophyll content (SPAD)"
genotype_limits = c("Y:Y.Y:Y",
                    "X:X.Y:Y",
                    "Y:Y.X:X",
                    "X:X.X:X")
genotype_names <- c(
  `X:X.X:X` = "aabb",
  `Y:Y.X:X` = "AAbb",
  `X:X.Y:Y` = "aaBB",
  `Y:Y.Y:Y` = "AABB"
)

#Scales for ggplot2
genotype_scale_fill <- scale_fill_manual(values = c("Y:Y.Y:Y" = "#ffffffff",
                                         "X:X.Y:Y" = "#deebf7ff",
                                         "Y:Y.X:X" = "#9ecae1ff",
                                         "X:X.X:X" = "#3182bdff"), 
                                         labels = genotype_names)
genotype_scale_col <- scale_colour_manual(values = c("Y:Y.Y:Y" = "#000000ff",
                                                   "X:X.Y:Y" = "#deebf7ff",
                                                   "Y:Y.X:X" = "#9ecae1ff",
                                                   "X:X.X:X" = "#3182bdff"), 
                                        labels = genotype_names)
genotype_scale_col_2 <- scale_colour_manual(values = c("Y:Y.Y:Y" = "#000000ff",
                                                     "X:X.X:X" = "#3182bdff"), 
                                          labels = genotype_names)
genotype_scale_x <- scale_x_discrete(limits = c("Y:Y.Y:Y",
                                                "X:X.Y:Y",
                                                "Y:Y.X:X",
                                                "X:X.X:X"), labels = genotype_names)
genotype_scale_shape <- scale_shape(limits = c("Y:Y.Y:Y",
                                      "X:X.X:X"),
                                    labels = genotype_names)
date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
week_scale <- scale_x_continuous(breaks = seq(0,9,1))
subgenome_scale <- scale_colour_manual(values=c("A"="#3182bdff","B"="#62c2acff","D"="#e5f58dff"))
spad_scale <- scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60))

variable_list_full <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen",
                        "Peduncle_senescence_DAH", "Heading_date", "Height", "Tiller_1",
                        "Seed_num",    "Weight",      "TGW", "O_length",
                        "`Predicted Moisture %`",          "`Predicted Protein Dry basis %`",
                        "`Predicted Starch Dry basis %`",  "`Predicted NDF Dry basis %`",
                        "`Predicted Dry gluten As is %`", "`Predicted Hardness -`")
variable_labels <- c(
  Leaf_senescence_DAH = "Days to 25% leaf senescence",
  TT_SPAD30 = "Days to onset of leaf senescence (SPAD=30)",
  Dur_Leaf_Sen = "Duration of leaf senescence",
  AUC_Leaf_Sen = "AUC of leaf senescence",
  Peduncle_senescence_DAH = "Days to 100% peduncle senescence",
  Heading_date = "Heading date",
  Height = "Height (cm)",
  Tiller_1 = "Main tiller number",
  Seed_num = "Seed number",
  Weight = "Seed mass (g)",
  TGW = "Thousand grain weight (g)",
  O_length = "Average grain length (mm)",
  `Predicted Moisture %` = "Predicted Moisture (%)",
  `Predicted Protein Dry basis %` = "Predicted Protein Dry basis (%)",
  `Predicted Starch Dry basis %` = "Predicted Starch Dry basis (%)",
  `Predicted NDF Dry basis %` = "Predicted NDF Dry basis (%)",
  `Predicted Dry gluten As is %` = "Predicted Dry gluten (%)",
  `Predicted Hardness -` = "Predicted Hardness",
  Light = "Distance from Light",
  SPAD = "Flag leaf chlorophyll content (SPAD)"
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OUTLIERS
print(paste("Gene:", gene, "Date:", date))

# summary(Kronos_wide_yday)

Kronos_wide_yday <- Kronos_phenotyping_yday %>%
  mutate(A = factor(sub("\\?", NA, A)), B = factor(sub("\\?", NA, B)), Genotype = interaction(A, B)) %>%
  mutate(Column = as.numeric(Block)*2 + (as.numeric(Position)-1) %/% 4 -2,
         Block = as.factor(Block)) %>%
  tidyr::extract(Line_name, "Cross", regex = "^(X[0-9]{3})", remove = FALSE)

Kronos_SPAD_yday <- Kronos_SPAD_yday %>%
  mutate(A = factor(sub("\\?", NA, A)), B = factor(sub("\\?", NA, B)), Genotype = interaction(A, B)) %>%
  mutate(Column = as.numeric(Block)*2 + (as.numeric(Position)-1) %/% 4 -2,
         Block = as.factor(Block)) %>%
  tidyr::extract(Line_name, "Cross", regex = "^(X[0-9]{3})", remove = FALSE)

Kronos_outliers <- Kronos_wide_yday %>%
  filter(is.na(Plant_num)|grepl("snapped", Comments)|grepl("damage", Comments)|is.na(Genotype))

# View(Kronos_outliers)

#Remove outliers and NIR readings I don't trust very much as not enough seed
Kronos_analysis <- Kronos_wide_yday %>%
  anti_join(Kronos_outliers, by = "Index") %>%
  mutate(across(
    `Predicted Moisture %`:`Predicted Hardness -`,
    .fns =  ~ ifelse((Weight > 16 & Sowing_date == "2021-12-01"), .x, NA)
  ))

Kronos_SPAD_analysis <- Kronos_SPAD_yday %>%
  anti_join(Kronos_outliers, by = "Index")