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
#'Repeat for NAC5A
#'Also see qPCR_graphs_230321.R
#'
#'Built under R 4.0.5
#'
#'
#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(tidyverse)
library(reshape2)
setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics/NAC expression R")
source("QPCR_functions.R")

#that crucial part
gene = "NAC-5"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R/Overexpression_thesis_graphs_source_2023_07_28.R', echo=TRUE)

#'EDIT section
#'Choose file names and parameters

setwd("C:/Users/evansc/OneDrive - Norwich Bioscience Institutes/NAC transgenics/NAC expression R")
genotype_key <- read_excel("NAC-5A_genotype_key.xlsx")

#RUN section
#Organise theme and font sizes
# theme_set(theme_bw())
# my_theme <- theme(axis.title = element_text(size = 18), plot.title = element_text(size = 32), axis.text = element_text(size = 16), 
#                   panel.grid = element_blank(), legend.position = "right", legend.text = element_text(size = 16),
#                   legend.title = element_text(size=18)) #Axis text reduced from 24 to 18 for fitting
# genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
#Line_name_labs <- c(CTA5.con1="5.con", CTA8.con1="8.con", CTA5.4="5.4", CTA8.9="8.9", CTA8.26="8.26", CTA8.2="8.2", CTA8.3="8.3", CTA8.1 = "8.1", CTA8.10 = "8.10", CTA8.23 = "8.23", CTA8.25 = "8.25", CTA5.5 = "5.5", CTA5.1 = "5.1", CTA5.2="5.2")

#READ
setwd("C:/Users/evansc/OneDrive - University of Warwick/NAC transgenics")

Primer_5 <- read_csv("qPCR_CENAC55.csv")

Primer_7 <- read_csv("qPCR_CENAC57.csv")




#Filter and add Line_name and Parental Line_name data 
#Separate step irrelevant for pooled samples
analysis_filtered_5 <- Primer_5 %>%
mutate(Line_name = str_remove(Line_name, "CTA")) %>%
  left_join(genotype_key, by = c("Line_name"))
analysis_filtered_7 <- Primer_7 %>%
  mutate(Line_name = str_remove(Line_name, "CTA")) %>%
  left_join(genotype_key, by = c("Line_name"))

#Bar chart please. Here's the pared-down simple version.

line_name_limits <- c("5.con1", "8.con1", "5.4", "8.9", "8.10", "8.23", "8.25", "8.26",
                      "5.1", "5.2", "5.5", "8.1", "8.2", "8.3")

#serious fudging
analysis_filtered_7 <- analysis_filtered_7[c(1:12, 14:16),]
analysis_filtered_5 <- analysis_filtered_5[c(1:12, 14:16),]

plot1 <- ggplot(data = analysis_filtered_7, aes(x = Line_name, y = DDCt, fill = Parental_genotype)) +
  geom_col(col = "black") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_errorbar(aes(ymin = DDCt - Min_y, ymax = DDCt + Max_y)) +
  labs(y = str_wrap("Relative transcript level of NAC5-A1", width = 25), x = "Line") +
  genotype_scale_fill_multi +
  scale_x_discrete(limits = line_name_limits) +
  scale_y_continuous(labels = function(x) paste(x, ".0", sep="")) +
  theme(legend.position = "none")

plot2 <- ggplot(data = analysis_filtered_5, aes(x = Line_name, y = DDCt, fill = Parental_genotype)) +
  geom_col(col = "black") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_errorbar(aes(ymin = DDCt - Min_y, ymax = DDCt + Max_y)) +
  labs(y = str_wrap("Relative transcript level of 3*FLAG-NAC5-A1 construct", width = 25), x = "Line") +
  genotype_scale_fill_multi +
  scale_x_discrete(limits = line_name_limits) +
  theme(legend.position = "none")

ggarrange(plot2, plot1, nrow = 2)

setwd("C:/Users/evansc/OneDrive - Norwich Bioscience Institutes/NAC transgenics/NAC expression R")
save_panel_svg(ggarrange(plot2, plot1, nrow = 2), plot_name = "Expression_pooled_DDCt", gene = gene,
               n_panel_cols = 1, n_panel_rows = 2, ratio = 0.75)



plot_5 <- ggplot(data = Primer_5, aes(x = Line_name, y = DDCt, fill = Parental_genotype)) +
  geom_col() +
  #scale_x_discrete(limits = Primer_5$Line_name, labels = genotype_labs) +
  my_theme +
  genotype_scale +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45)) +
  labs(x = "Line", y = "Relative transcript level of construct") +
  geom_errorbar(aes(ymin = DDCt - Min_y, ymax = DDCt + Max_y), width = 0.5)

plot_5 + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.y = element_text())



plot_7 <- ggplot(data = Primer_7, aes(x = Line_name, y = DDCt, fill = Parental_genotype)) +
  geom_col() +
  #scale_x_discrete(limits = Primer_7$Line_name, labels = genotype_labs) +
  my_theme +
  genotype_scale +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Line", y = "Relative transcript level of NAC-5A") +
  geom_errorbar(aes(ymin = DDCt - Min_y, ymax = DDCt + Max_y), width = 0.5)

plot_7 + scale_y_continuous(labels = function(x) paste(x, ".0", sep="")) #So that the digits are the same as the other plot!
