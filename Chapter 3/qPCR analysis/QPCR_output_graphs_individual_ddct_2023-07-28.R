# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'31.05.21
#'Draw graphs for qPCR data
#'Here, NAC-5A qPCR experiment with individual samples
#'
#'01.06.21
#'Add theme and graphs for Plant Theme seminar
#'
#'03/06/21
#'Use for NAP-A pooled samples run 25/05/21. Simplified plot section.
#'
#'21/09/21
#'Use for NAP-A individual samples run 15/09/21 and 16/09/21.
#'
#'28/10/21
#'Use data from Pfaffl 2001 method
#'
#'26/09/22
#'Specifically plot individual plants that were used for the next experiment
#'
#'29/09/22
#'Omit CTA8.9
#'
#'28/07/23
#'Make pretty graphs for thesis
#'
#'05/12/23
#'Corrections: use new data from QPCR_run_me_NAPA_2023-12-05.R and QPCR_run_me_NAC5A_2023-12-05.R
#'
#'Built under R 4.0.5
#'
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(tidyverse)
library(reshape2)
library(lubridate)
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
source("QPCR_functions.R")

#'EDIT section
#'Choose file names and parameters
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")

#Date is the date raw data were run through QPCR_run_me.R
date = "2023-12-06"

genotype_key <- read_excel("NAPA_genotype_key.xlsx")
genotype_key2 <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/NAC-5A_genotype_key.xlsx")
genotype_key <- rbind(genotype_key, genotype_key2)

NAP_analysis<- read_csv(paste("qPCR_NAPA_individual_analysis_pfaffl_", date, ".csv", sep = "")) #Pfaffl method fold_change

NAC5_analysis<- read_csv(paste("qPCR_NAC5A_individual_analysis_pfaffl_", date, ".csv", sep = "")) #Pfaffl method fold_change

#RUN section
#Organise theme and font sizes
# theme_set(theme_bw())
# my_theme <- theme(axis.title = element_text(size = 18), plot.title = element_text(size = 32), axis.text = element_text(size = 16), 
#                   panel.grid = element_blank(), legend.position = "right", legend.text = element_text(size = 16),
#                   legend.title = element_text(size=18)) #Axis text reduced from 24 to 18 for fitting
# genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
# #Line_name_labs <- c(CTA5.con1="5.con", CTA8.con1="8.con", CTA5.4="5.4", CTA8.9="8.9", CTA8.26="8.26", CTA8.2="8.2", CTA8.3="8.3", CTA8.1 = "8.1", CTA8.10 = "8.10", CTA8.23 = "8.23", CTA8.25 = "8.25", CTA5.5 = "5.5", CTA5.1 = "5.1", CTA5.2="5.2")

#Filter and add Line_name and Parental Line_name data 
analysis_combined <- full_join(NAP_analysis, NAC5_analysis)

analysis_filtered <- analysis_combined %>%
  filter(Count.x >=2 & Count.y >=2) %>%
  mutate(Sample_2 = str_remove(Sample, "CTA"), Sample_3 = str_remove(Sample, "CTA")) %>% #Don't want to lose sample column
  tidyr::separate(Sample_2, into = c("Line_name", "PlantNum"), sep = "-") %>%
  left_join(genotype_key, by = "Line_name")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Run for NAP
#Sort out the theme labels
gene = "NAP"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/Overexpression_thesis_graphs_source_2023_07_28.R', echo=TRUE)

# #Individual parent plants of T2 experiment
# analysis_filtered <- analysis_filtered %>%
#   filter(Sample %in% c("CTA5.4-9", "CTA5.4-2", "CTA10.5-6", "CTA10.5-7", "CTA10.9-5", "CTA10.9-2"))

#Boxplot
boxplot_CENAP_7F_2R <- analysis_filtered %>%
  filter(Primer.y == "CENAP-7F-2R") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line", y = str_wrap("Relative transcript level of NAP-A1", width = 25))
boxplot_CENAP_7F_2R

boxplot_CENAP_8F_13R <- analysis_filtered %>%
  filter(Primer.y == "CENAP-8F-13R") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  #tweak breaks
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line", y = str_wrap("Relative transcript level of NAP-1, all homoeologs", width = 25))
boxplot_CENAP_8F_13R

boxplot_CENAP_5F_5R <- analysis_filtered %>%
  filter(Primer.y == "CENAP-5F-5R") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line", y = str_wrap("Relative transcript level of 3*FLAG-NAP-1 construct", width = 25))
boxplot_CENAP_5F_5R


#Save all my plots
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")

#ggsave(filename = paste("plot_CENAP_7F_2R_",today(),".svg",sep=""), plot_CENAP_7F_2R, width = 90, height = 90, units = "mm", device = "svg")

ggarrange(boxplot_CENAP_5F_5R, boxplot_CENAP_7F_2R, boxplot_CENAP_8F_13R)
save_panel_svg(ggarrange(boxplot_CENAP_5F_5R, boxplot_CENAP_7F_2R, boxplot_CENAP_8F_13R), plot_name = "Expression_individual_pfaffl", gene = "NAP",
              n_panel_cols = 2, n_panel_rows = 2, ratio = 1.25)
#ta da

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Repeat for NAC5
#Sort out the theme labels
gene = "NAC-5"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/Overexpression_thesis_graphs_source_2023_07_28.R', echo=TRUE)

#remove the extreme outlier
analysis_filtered <- analysis_filtered %>%
  filter(Sample != "CTA8.con1-1")

#Boxplot
boxplot_CENAP_7F_2R <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-7") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  #tweak breaks
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line", y = str_wrap("Relative transcript level of NAC5-A1", width = 25))
boxplot_CENAP_7F_2R

boxplot_CENAP_8F_13R <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-3") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  #tweak breaks
  scale_y_continuous(breaks = seq(0, 6, 1)) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line", y = str_wrap("Relative transcript level of NAC5-1, all homoeologs", width = 25))
boxplot_CENAP_8F_13R

boxplot_CENAP_5F_5R <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-5") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  #tweak breaks
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "none") +
  labs(x = "Line", y = str_wrap("Relative transcript level of 3*FLAG-NAC5-1 construct", width = 26))
boxplot_CENAP_5F_5R


#Save all my plots
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")

ggarrange(boxplot_CENAP_5F_5R, boxplot_CENAP_7F_2R, boxplot_CENAP_8F_13R)
save_panel_svg(ggarrange(boxplot_CENAP_5F_5R, boxplot_CENAP_7F_2R, boxplot_CENAP_8F_13R), plot_name = "Expression_individual_pfaffl", gene = "NAC5",
               n_panel_cols = 2, n_panel_rows = 2, ratio = 1.25)

#Get legend

boxplot_CENAP_7F_2R <- analysis_filtered %>%
  filter(Primer.y == "CENAC5-7") %>%
  ggplot(aes(x = Line_name, y = fold_change, fill = Parental_genotype, group = Line_name)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_x_discrete(limits = line_name_limits) +
  my_theme +
  genotype_scale_fill_multi +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Line", y = str_wrap("Relative transcript level of NAC5-A1", width = 25))
boxplot_CENAP_7F_2R
  
ggsave(filename = paste("boxplot_CENAP_7F_2R_",today(),".svg",sep=""), boxplot_CENAP_7F_2R, width = 90, height = 90, units = "mm", device = "svg")

