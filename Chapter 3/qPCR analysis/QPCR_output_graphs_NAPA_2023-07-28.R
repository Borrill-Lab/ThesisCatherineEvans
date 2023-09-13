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
#'28/07/23
#'Make pretty graphs for thesis
#'
#'Built under R 4.0.5
#'
#'
#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(tidyverse)
library(reshape2)
setwd("C:/Users/evansc/OneDrive - University of Warwick/NAC transgenics/NAC expression R")
source("QPCR_functions.R")

#that crucial part
gene = "NAP"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/Overexpression_thesis_graphs_source_2023_07_28.R', echo=TRUE)

#'EDIT section
#'Choose file names and parameters
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")


analysis_CENAP_2<- read_csv("qPCR_NAPA_CENAP-2_analysis_20210603.csv")
analysis_CENAP_5 <- read_csv("qPCR_NAPA_CENAP-5_analysis_20210603.csv")
genotype_key <- read_excel("NAPA_genotype_key.xlsx")

#RUN section
#Organise theme and font sizes
# theme_set(theme_bw())
# my_theme <- theme(axis.title = element_text(size = 18), plot.title = element_text(size = 32), axis.text = element_text(size = 16), 
#                   panel.grid = element_blank(), legend.position = "right", legend.text = element_text(size = 16),
#                   legend.title = element_text(size=18)) #Axis text reduced from 24 to 18 for fitting
# genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
#Line_name_labs <- c(CTA5.con1="5.con", CTA8.con1="8.con", CTA5.4="5.4", CTA8.9="8.9", CTA8.26="8.26", CTA8.2="8.2", CTA8.3="8.3", CTA8.1 = "8.1", CTA8.10 = "8.10", CTA8.23 = "8.23", CTA8.25 = "8.25", CTA5.5 = "5.5", CTA5.1 = "5.1", CTA5.2="5.2")


#Filter and add Line_name and Parental Line_name data 
#Separate step irrelevant for pooled samples
analysis_CENAP_2_filtered <- analysis_CENAP_2%>%
  filter(Count.x >=2 & Count.y >=2 & Primer.y == "CENAP-2") %>%
  mutate(Sample_2 = str_remove(Sample, "CTA")) %>% #Don't want to lose sample column
  #separate(Sample_2, into = c("Line_name", "PlantNum"), sep = "-") %>%
  mutate(Line_name = Sample) %>%
  left_join(genotype_key, by = c("Sample_2" = "Line_name"))

analysis_CENAP_5_filtered <- analysis_CENAP_5 %>%
  filter(Count.x >=2 & Count.y >=2 & Primer.y == "CENAP-5") %>%
  mutate(Sample_2 = str_remove(Sample, "CTA")) %>% #Don't want to lose sample column
  #separate(Sample_2, into = c("Line_name", "PlantNum"), sep = "-") %>%
  #complete(Line_name, PlantNum, Primer.x, Primer.y) %>%
  mutate(Line_name = Sample) %>%
  left_join(genotype_key, by = c("Sample_2" = "Line_name"))

#Bar chart please. Here's the pared-down simple version.

line_name_limits <- c("10.con1", "7.19", "10.3", "10.5", "10.9", "7.7", "7.8", "7.9","10.1", "10.8")

plot1 <- ggplot(data = analysis_CENAP_2_filtered, aes(x = Sample_2, y = two_power_DDCt, fill = Parental_genotype)) +
  geom_col(col = "black") +
  my_theme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_errorbar(aes(ymin = error_min, ymax = error_max)) +
  labs(y = str_wrap("Relative transcript level of NAP-A1", width = 25), x = "Line") +
  genotype_scale_fill_multi +
  scale_x_discrete(limits = line_name_limits) +
  theme(legend.position = "none")

plot2 <- ggplot(data = analysis_CENAP_5_filtered, aes(x = Sample_2, y = two_power_DDCt, fill = Parental_genotype)) +
  geom_col(col = "black") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_errorbar(aes(ymin = error_min, ymax = error_max)) +
  labs(y = str_wrap("Relative transcript level of 3*FLAG-NAP-1 construct", width = 25), x = "Line") +
  genotype_scale_fill_multi +
  scale_x_discrete(limits = line_name_limits) +
  theme(legend.position = "none")

ggarrange(plot2, plot1, nrow = 2)

setwd("C:/Users/evansc/OneDrive - Norwich Bioscience Institutes/NAC transgenics/NAC expression R")
save_panel_svg(ggarrange(plot2, plot1, nrow = 2), plot_name = "Expression_pooled_DDCt", gene = gene,
               n_panel_cols = 1, n_panel_rows = 2, ratio = 0.75)
