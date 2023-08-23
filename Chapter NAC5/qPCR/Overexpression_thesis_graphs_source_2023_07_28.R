# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'Overexpression thesis graphs SOURCE
#'
#'28/09/2022
#'Build graphs for Overexpression data for Crop Gen Seminar on 03/10/2022
#'Based on Overexpression_exploratory_2022-09-06
#'
#'06/10/2022
#'Update for ASM poster
#'
#'20/06/2023
#'Adapt for graphs to go in thesis for Overexpression data
#'
#'28/07/2023
#'Themes only for qPCR data
#'
#'This SOURCE script loads packages and data. Last run on NBI desktop.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(car) #Anova for non-orthogonal data
#require(emmeans) #marginal means estimates for when p-values are not enough
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#THEME
my_theme <- theme_linedraw(base_size = 20) + theme(panel.grid = element_blank(),
        line = element_line(colour = "#000000"),
        axis.text = element_text(colour = "#000000"),
        strip.background = element_rect(fill = "#cccccc"),
        strip.text = element_text(colour = "#000000"),
        legend.position = "right")
theme_set(my_theme)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#AXIS SCALES
genotype_lab <- c(control = "Non-transformed control", single = "Single-copy", multi = "Multi-copy")

if(gene == "NAP"){
line_name_limits <- c("7.con1", "10.con1", "10.5", "10.9", "10.1", "7.9")
} else {
  line_name_limits <- c("8.con1", "5.4", "8.9", "8.26", "8.2", "8.3")
}
genotype_scale_fill <- scale_fill_manual(values=c("hom"="#fecc5cff", "null"="#ffffffff"), limits = c("null", "hom"), labels = genotype_lab)
genotype_scale_col <- scale_colour_manual(values=c("hom"="#fecc5cff", "null"="#000000ff"), limits = c("null", "hom"), labels = genotype_lab)
genotype_scale_shape <- scale_shape(limits = c("null", "hom"), labels = genotype_lab)

genotype_scale_fill_multi <- scale_fill_manual(values=c("control"="#ffffffff", "single"="#fecc5cff", "multi" = "#e31a1cff"), labels = genotype_lab)
genotype_scale_col_multi <- scale_colour_manual(values=c("single"="#fecc5cff", "control"="#000000ff", "multi" = "#e31a1cff"), limits = c("single", "control", "multi"))
genotype_scale_shape_multi <- scale_shape(limits = c("single", "control", "multi"))

date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
week_scale <- scale_x_continuous(breaks = seq(0,9,1))
spad_scale <- scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60), expand = expansion(0, 0))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#AXIS LABELS

Line_name_lab <- paste(gene, "Line name")
leaf_senescence_lab <- "Days to 25% leaf senescence"
peduncle_senescence_lab <- "Days to 100% peduncle senescence"
SPAD_lab <- "Flag leaf chlorophyll content (SPAD)"

variable_list_full <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen",
                        "Peduncle_senescence_DAH", "Heading_date", "Height", "Tiller_1",
                        "Seed_num",    "Weight",      "TGW", "O_length")
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
  Light = "Distance from Light",
  SPAD = "Flag leaf chlorophyll content (SPAD)",
  Week = "Weeks after heading",
  Line_name = paste(gene, "Line name")
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
  
