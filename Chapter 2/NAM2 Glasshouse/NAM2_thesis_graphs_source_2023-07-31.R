# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'NAM2 SEB graphs SOURCE
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
#'19/06/2023
#'Rewrite for NAM2 glasshouse Exp2 data
#'
#'30/06/2023
#'Simplify genotype labels to single / double / triple mutant
#'
#'31/07/2023
#'Pretty graphs for thesis. Themes only, doesn't load data files.
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

get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()

setwd(paste(local_file_path, "NAC transgenics/2021-03 NAC phenotyping", sep = ""))
source("Functions_for_phenotype_data.R")

#THEME

#THESIS THEME
my_theme <- theme_linedraw(base_size = 20) + theme(panel.grid = element_blank(),
                                                   line = element_line(colour = "#000000"),
                                                   axis.text = element_text(colour = "#000000"),
                                                   strip.background = element_rect(fill = "#cccccc"),
                                                   strip.text = element_text(colour = "#000000"),
                                                   legend.position = "right")
theme_set(my_theme)

#LABELS
genotype_lab <- paste(gene, "Genotype")
leaf_senescence_lab <- "Days to 25% leaf senescence"
peduncle_senescence_lab <- "Days to 100% peduncle senescence"
SPAD_lab <- "Flag leaf chlorophyll content (SPAD)"

#Genotype
genotype_limits = c("G1",
                    "G5",
                    "G3",
                    "G2",
                    "G7",
                    "G6",
                    "G4",
                    "G8")
genotype_names <- c(
  G1 = "WT",
  G2 = "dd",
  G3 = "bb",
  G4 = "bbdd",
  G5 = "aa",
  G6 = "aadd",
  G7 = "aabb",
  G8 = "aabbdd"
)

#Scales for ggplot2
#Haven't fully thought through colour scheme here
genotype_scale_fill <- scale_fill_brewer(type = "seq", palette = "Blues",
                                         limits = genotype_limits,
                                         labels = genotype_names)
genotype_scale_col <- scale_colour_brewer(type = "seq", palette = "Blues",
                                          limits = genotype_limits,
                                        labels = genotype_names)
genotype_scale_col_2 <- scale_colour_manual(values = c("G1" = "#000000",
                                                     "G8" = "#2171b5"),
                                          labels = genotype_names)
genotype_scale_x <- scale_x_discrete(limits = genotype_limits, labels = genotype_names)
genotype_scale_x_2 <- scale_x_discrete(limits = c("G1",
                                                  "G8"),
                                       labels = genotype_names)
genotype_scale_shape_2 <- scale_shape(limits = c("G1",
                                      "G8"),
                                    labels = genotype_names)

#'Simplify genotype labels to single / double / triple mutant
dosage_limits = c("Wild type",
                  "single",
                  "double",
                  "triple"
                  )
dosage_names <- c(
  `Wild type` = "Wild type",
  single = "NAM-2 single mutant",
  double = "NAM-2 double mutant",
  triple = "NAM-2 triple mutant"
)

#Scales for ggplot2
#Colour scheme based on ColourBrewer "Blues" 4-colour palette, replacing WT with white/black
dosage_scale_fill <- scale_fill_manual(values = c(
  "Wild type" = "#ffffff",
  "single" = "#bdd7e7",
  "double" = "#6baed6",
  "triple" = "#2171b5"
),
limits = dosage_limits,
labels = dosage_names)
dosage_scale_col <- scale_colour_manual(values = c(
  "Wild type" = "#000000",
  "single" = "#bdd7e7",
  "double" = "#6baed6",
  "triple" = "#2171b5"
),
limits = dosage_limits,
labels = dosage_names)
dosage_scale_col_2 <- scale_colour_manual(values = c("Wild type" = "#000000",
                                                     "triple" = "#2171b5"),
                                          labels = dosage_names)
dosage_scale_x <- scale_x_discrete(limits = dosage_limits, labels = dosage_names)
dosage_scale_shape_2 <- scale_shape(limits = c("Wild type",
                                               "triple"),
                                    labels = dosage_names)

#Other variables
#Careful with limits on continuous scales as they CUT OFF THE DATA VALUES. Use coord_cartesian instead
date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
week_scale <- scale_x_continuous(breaks = seq(0,9,1))
subgenome_scale <- scale_colour_manual(values=c("A"="#3182bdff","B"="#62c2acff","D"="#e5f58dff"))
spad_scale <- scale_y_continuous(breaks = c(0, 20, 40, 60))

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
  SPAD = "Flag leaf chlorophyll content (SPAD)",
  Flag_leaf_length = "Flag leaf length (mm)",
  Flag_leaf_width = "Flag leaf width (mm)"
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#For the field data
DDAH_lab = "Thermal Time after Heading (°C Days)"
variable_list_field = c("AUC_Leaf_Sen", "Dur_Leaf_Sen", "Leaf_TT10", "Leaf_TT30", "Leaf_TT50", "Leaf_TT70", "Leaf_TT90",
                       "`Peduncle_senescence.2023-07-17`",
                       "`Peduncle_senescence.2023-07-20`",
                       "`Cover.2023-03-06`",
                       "`Height.2023-05-05`",
                       "EE_date")
variable_list_unquote_field = c("AUC_Leaf_Sen", "Dur_Leaf_Sen", "Leaf_TT10", "Leaf_TT30", "Leaf_TT50", "Leaf_TT70", "Leaf_TT90",
                        "Peduncle_senescence.2023-07-17",
                        "Peduncle_senescence.2023-07-20",
                        "Cover.2023-03-06",
                        "Height.2023-05-05",
                        "EE_date")
#Modified for similar peduncle senescence stats
variable_list_truncated_field <- paste(str_trunc(str_split(variable_list_unquote_field, "\\.", simplify = TRUE)[,1], width = 12, side = "right", ellipsis = "-"),
                                       str_trunc(str_split(variable_list_unquote_field, "\\.", simplify = TRUE)[,2], width = 3, side = "left", ellipsis = "-"),
                                       sep = "")
variable_labels_field = c(AUC_Leaf_Sen = "AUC of leaf senescence",
                          Dur_Leaf_Sen = "Duration of leaf senescence (°C days)",
                          Leaf_TT10 = "°C days to 10% Leaf senescence",
                          Leaf_TT30 = "°C days to 30% Leaf senescence",
                          Leaf_TT50 = "°C days to 50% Leaf senescence",
                          Leaf_TT70 = "°C days to 70% Leaf senescence",
                          Leaf_TT90 = "°C days to 90% Leaf senescence",
                          `Peduncle_senescence.2023-07-17` = "Senesced peduncles (out of 10) 17/07",
                          `Peduncle_senescence.2023-07-20` = "Senesced peduncles (out of 10) 20/07",
                          `Cover.2023-03-06` = "Cover 2023-03-06",
                          `Height.2023-05-05`= "Height 2023-05-05",
                          EE_date = "Date of ear emergence"
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Automate standard plot options
save_small_svg <- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 160,
    height = 90,
    units = "mm",
    device = "svg"
  )
}
save_big_svg <- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 320,
    height = 180,
    units = "mm",
    device = "svg"
  )
}

save_pdf_1 <- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_", lubridate::today(), ".pdf", sep = ""),
    plot,
    width = 190,
    height = 110,
    units = "mm",
    device = "pdf"
  )
}

save_poster_pdf_2 <- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_poster_", lubridate::today(), ".pdf", sep = ""),
    plot,
    width = 90,
    height = 110,
    units = "mm",
    device = "pdf"
  )
}

save_A4_svg<- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_A5_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 210,
    height = 297,
    units = "mm",
    device = "svg"
  )
}

save_A5_svg <- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_A4_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 210,
    height = 148,
    units = "mm",
    device = "svg"
  )
}

save_panel_svg <- function(plot, plot_name, gene, n_panel_cols, n_panel_rows, ratio = 0.75){
  ggsave(
    filename = paste(gene, "_", plot_name, "_", lubridate::today(), ".svg", sep = ""),
    plot,
    width = 50 + 160,
    height = 50 + (160/n_panel_cols) * ratio * n_panel_rows,
    units = "mm",
    device = "svg"
  )
}