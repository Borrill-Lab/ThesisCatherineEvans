library(readr); library(ggplot2); library(reshape2); library(dplyr)
theme_set(theme_bw())
setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-03 NAC phenotyping")
SPAD_table <- read_csv("NAPA_SPAD_2021-10-29.csv")

#Organise theme and font sizes
theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 24), plot.title = element_text(size = 32), axis.text = element_text(size = 16), 
                  panel.grid = element_blank(), legend.position = "right", legend.text = element_text(size = 16),
                  legend.title = element_text(size=24))

#Labels.
genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
genotype_scale_colour <- scale_colour_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
line_scale <- scale_color_manual(values=c("CTA7.con1"="#ca0020","CTA10.con1" = "#f4a582", "CTA10.5" ="#d1e5f0","CTA10.9"="#92c5de", "CTA7.9"="#4393c3", "CTA10.1"="#2166ac"))
SPAD_labs <- c(`0`="0", `1`="18-20", `2`="25-27", `3`="32-34", `4`="39-41", `5`="46-48")


SPAD_plot <- ggplot(SPAD_table, aes(x=Week, y=SPAD, group = Line_name, shape=Parental_genotype, color=Parental_genotype))
SPAD_plot <- SPAD_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  labs(x="Days after Heading", y="SPAD", shape = "Parental_genotype", color="Parental_genotype") +
  scale_x_continuous(breaks = seq(0,5,1), labels=SPAD_labs) +
  my_theme
SPAD_plot + genotype_scale_colour+
  theme(legend.position = "none")

#Use NAC-5 data for below from NAC5A_poster_graphs_230321
SPAD_plot <- ggplot(filter(Table_SPAD, variable != "SPAD_4648DAH"), aes(x=variable, y=SPAD, group = Line_name, shape=Parental_genotype, color=Parental_genotype))
SPAD_plot <- SPAD_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  labs(x="Days after Heading", y="SPAD", shape = "Line", color="Line") +
  scale_x_discrete(labels=SPAD_labs) +
  my_theme
SPAD_plot + genotype_scale_colour +
#   scale_shape(limits = c("CTA5.con1", "CTA8.con1", "CTA8.2", "CTA8.3"), labels = genotype_labs) +
  theme(legend.position = "none")

