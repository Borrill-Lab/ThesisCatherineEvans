# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#Source script for
#NAM_figures_part1
#NAM_figures_part2

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#LOAD
require(tidyverse)
require(car) #Anova for non-orthogonal data
# require(emmeans) #marginal means estimates for when p-values are not enough
# require(arm)
require(GGally) #ggpairs()
require(ggpubr) #stat_cor()
theme_set(theme_bw())

# setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/RAGT Field Senescence 2021/Data")
# setwd("U:/RAGT Senescence, protein and yield/RAGT Field Senescence 2021/Data")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/RAGT Field Senescence 2021/Data")

setwd("U:/Field senescence re-analysis/Processed_data")

Long_data <- read_csv("Field_results_NAM_2021_long_2023-08-23.csv",
                      col_types = cols("NAMA1"= col_factor(),
                                       "NAMB1"= col_factor())
)

Wide_data <- read_csv("Field_results_NAM_2021_wide_DDAH_2023-08-23.csv",
                      col_types = cols("NAMA1"= col_factor(),
                                       "NAMB1"= col_factor())
)
setwd("..")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#TIDY
#EDIT
#Add some extra columns for plotting
Wide_data <- Wide_data %>%
  mutate(Control = ifelse(Line_name %in% c("SKYFALL", "RGT BOWEN"), TRUE, FALSE),
         NAMA1 = fct_relevel(NAMA1, "Non functional mutant"),
         Genotype = interaction(NAMA1, NAMB1),
         Protein_yield = Yield * Protein/100)

Wide_data_control <- Wide_data %>%
  filter(Control == TRUE)

Wide_data_no_control <- Wide_data %>%
  filter(Control == FALSE)

Long_data <- Long_data %>%
  mutate(Control = ifelse(Line_name %in% c("SKYFALL", "RGT BOWEN"), TRUE, FALSE),
         NAMA1 = fct_relevel(NAMA1, "Non functional mutant"),
         Genotype = interaction(NAMA1, NAMB1))

Long_data_control <- Long_data %>%
  filter(Control == TRUE)

Long_data_no_control <- Long_data %>%
  filter(Control == FALSE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#THESIS THEME
my_theme <- theme_linedraw(base_size = 20) + theme(panel.grid = element_blank(),
                                                   line = element_line(colour = "#000000"),
                                                   axis.text = element_text(colour = "#000000"),
                                                   strip.background = element_rect(fill = "#cccccc"),
                                                   strip.text = element_text(colour = "#000000"),
                                                   legend.position = "right")
theme_set(my_theme)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#LABELS
genotype_lab <- "NAM-1 Genotype"
leaf_senescence_lab <- "Leaf senescence (%)"
ear_senescence_lab <- "Ear senescence (%)"
genotype_limits <- c(
  "Non functional mutant.Wild type",
  "Wild type.Wild type",
  "Non functional mutant.Functional mutant",
  "Wild type.Functional mutant"
)
genotype_names <- c(
  `Non functional mutant.Wild type` = "aabbDD",
  `Wild type.Wild type` = "AAbbDD (WT)",
  `Non functional mutant.Functional mutant` = "aaBBDD",
  `Wild type.Functional mutant` = "AABBDD"
)

#Scales for ggplot2
# genotype_scale_fill <- scale_fill_brewer(type = "div", palette = "PRGn",
#                                          limits = genotype_limits, labels = genotype_names)
# genotype_scale_col <- scale_colour_brewer(type = "div", palette = "PRGn",
#                                           limits = genotype_limits, labels = genotype_names)
genotype_scale_col <- scale_colour_manual(values= c(`Non functional mutant.Wild type` = "#c2a5cf",
                                                    `Wild type.Wild type` = "#7b3294",
                                                    `Non functional mutant.Functional mutant` = "#a6dba0",
                                                    `Wild type.Functional mutant` = "#008837"),
                                          limits = genotype_limits, labels = genotype_names)
genotype_scale_fill <- scale_fill_manual(values= c(`Non functional mutant.Wild type` = "#c2a5cf",
                                                   `Wild type.Wild type` = "#7b3294",
                                                   `Non functional mutant.Functional mutant` = "#a6dba0",
                                                   `Wild type.Functional mutant` = "#008837"),
                                         limits = genotype_limits, labels = genotype_names)
genotype_scale_x <- scale_x_discrete(limits = genotype_limits, labels = str_wrap(genotype_names, width = 6))
date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))

da1m_lab = "Days after 1st May 2022"
DDAH_lab = "Thermal Time after Heading (oC Days)"
variable_list_leaf = c("AUC_Leaf_Sen", "Dur_Leaf_Sen", "Leaf_TT10", "Leaf_TT30", "Leaf_TT50", "Leaf_TT70", "Leaf_TT90",
                       "Leaf_curling", "Peduncle_Sen")
variable_labels_leaf = c("AUC Leaf senescence",
                    "Duration Leaf senescence (°C days)",
                    "Days to 10% Leaf senescence (°C days)",
                    "Days to 30% Leaf senescence (°C days)",
                    "Days to 50% Leaf senescence (°C days)",
                    "Days to 70% Leaf senescence (°C days)",
                    "Days to 90% Leaf senescence (°C days)",
                    "Leaf curling on 22/06/22 (%)",
                    "Senesced peduncles on 08/07/22 (out of 10)"
)

variable_list_ear = c("AUC_Ear_Sen", "Dur_Ear_Sen", "Ear_TT10", "Ear_TT30", "Ear_TT50", "Ear_TT70", "Ear_TT90")
variable_labels_ear = c("AUC Ear senescence",
                    "Duration Ear senescence (°C days)",
                    "Days to 10% Ear senescence (°C days)",
                    "Days to 30% Ear senescence (°C days)",
                    "Days to 50% Ear senescence (°C days)",
                    "Days to 70% Ear senescence (°C days)",
                    "Days to 90% Ear senescence (°C days)"
)

variable_list <- c("EE_date",
                   "Leaf_TT70",
                   "AUC_Leaf_Sen",
                   "Dur_Leaf_Sen", 
                   "Ear_TT70",
                   "AUC_Ear_Sen",
                   "Dur_Ear_Sen",
                   "`Peduncle_senescence.2022-07-12`")
variable_labels = c("Ear emergence date",
                    "Days to 70% Leaf senescence (°C days)",
                    "AUC Leaf senescence",
                         "Duration Leaf senescence (°C days)",
                    "Days to 70% Ear senescence (°C days)",
                         "AUC Ear senescence",
                         "Duration Ear senescence (°C days)",
                    "Senesced peduncles on 08/07/22 (out of 10)"
)

variable_list_unquote <- variable_list #NEEDS ADJUSTING WHEN NIR values added to escape weird characters
variable_list_truncated <- c(str_trunc(variable_list_unquote, width = 12, side = "right", ellipsis = "-"))


#For the field data
DDAH_lab = "Thermal Time after Heading (°C Days)"
variable_list_field <- c("EE",
                   "Leaf_TT70",
                   "AUC_Leaf_Sen",
                   "Dur_Leaf_Sen", 
                   "Ear_TT70",
                   "AUC_Ear_Sen",
                   "Dur_Ear_Sen",
                   "Yield",
                   "Protein",
                   "Protein_yield")
variable_list_unquote_field<- variable_list_field
variable_list_truncated_field <- c(str_trunc(variable_list_unquote_field, width = 12, side = "right", ellipsis = "-"))

variable_labels_field = c(AUC_Leaf_Sen = "AUC of leaf senescence",
                          Dur_Leaf_Sen = "Duration of leaf senescence (°C days)",
                          Leaf_TT10 = "°C days to 10% leaf senescence",
                          Leaf_TT30 = "°C days to 30% leaf senescence",
                          Leaf_TT50 = "°C days to 50% leaf senescence",
                          Leaf_TT70 = "°C days to 70% leaf senescence",
                          Leaf_TT90 = "°C days to 90% leaf senescence",
                          AUC_Ear_Sen = "AUC of ear senescence",
                          Dur_Ear_Sen = "Duration of ear senescence (°C days)",
                          Ear_TT70 = "°C days to 70% ear senescence",
                          `Peduncle_senescence.2022-07-12` = "Senesced peduncles (out of 10) 12/07",
                          EE_date = "Date of ear emergence",
                          EE_REL = "Ear emergence (days after 1st May)",
                          EE = "Ear emergence R1 (days after 1st May)",
                          Yield = "Yield",
                          Protein = "GPC R1",
                          Protein_yield = "Protein yield R1"
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
