#NAP-A graphs and statistical analysis
#29/10/2021
#'Initial graph sections
#'Run ANOVA and draw boxplots
#'Draw a line graph for SPAD values
#'
#'31/01/2022
#'Facet correlation plots by line
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)

#'SOURCE section
#'Packages and source files
library(readr)
library(ggplot2)
library(dplyr)
library(arm)
library(stringr)
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics needs tidying/NAC phenotyping/Functions_for_phenotype_data.R')
setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics/NAC phenotyping")
theme_set(theme_bw())

everything_table <- read_csv("NAPA_EVERYTHING_2021-10-29.csv")
SPAD_table <- read_csv("NAPA_SPAD_2021-10-29.csv")
priority_table <- everything_table %>%
  filter(Priority == TRUE) %>%
  mutate(Block = as.factor(Block))

#'EDIT section
#'Key parameters


#Labels
genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
genotype_scale_colour <- scale_colour_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))


#'PLOT section
#'Graphs galore
#'Not saved because I'll save them manually for now.
#'All the plants individually
ggplot(everything_table, aes(y=Heading_date, x=Plant_num, col=Line_name)) +
  geom_point(shape =2) +
  geom_point(mapping = aes(y=Leaf_senescence)) +
  geom_point(mapping = aes(y=Peduncle_senescence), shape = 3) +
  facet_wrap(facets = vars(Line_name), scales="free_x") +
  scale_x_continuous(breaks = seq(1,12,1))

#Boxplots for priority lines
basic_plots <- list(
  CENAP_5F_5R = ggplot(priority_table, aes(x=Line_name, y=`CENAP-5F-5R`, fill = Parental_genotype)),
  CENAP_7F_2R = ggplot(priority_table, aes(x=Line_name, y=`CENAP-7F-2R`, fill = Parental_genotype)),
  CENAP_8F_13R = ggplot(priority_table, aes(x=Line_name, y=`CENAP-8F-13R`, fill = Parental_genotype)),
  Heading_date = ggplot(priority_table, aes(x=Line_name, y=`Heading_date`, fill = Parental_genotype)),
  Leaf_senescence_DAH = ggplot(priority_table, aes(x=Line_name, y=`Leaf_senescence_DAH`, fill = Parental_genotype)),
  Peduncle_senescence_DAH = ggplot(priority_table, aes(x=Line_name, y=`Peduncle_senescence_DAH`, fill = Parental_genotype)),
  Leaf_senescence = ggplot(priority_table, aes(x=Line_name, y=`Leaf_senescence`, fill = Parental_genotype)),
  Peduncle_senescence = ggplot(priority_table, aes(x=Line_name, y=`Peduncle_senescence`, fill = Parental_genotype)),
  Grain_mass = ggplot(priority_table, aes(x=Line_name, y=`Grain_mass`, fill = Parental_genotype))
)
linear_models <- list(
  CENAP_5F_5R = lm(log(`CENAP-5F-5R`, 2) ~ Line_name + Block, data = priority_table),
  CENAP_7F_2R = lm(log(`CENAP-7F-2R`, 2) ~ Line_name + Block, data = priority_table),
  CENAP_8F_13R = lm(log(`CENAP-8F-13R`,2) ~ Line_name + Block, data = priority_table),
  Heading_date = lm(`Heading_date` ~ Line_name + Block, data = priority_table),
  Leaf_senescence_DAH = lm(`Leaf_senescence_DAH` ~ Line_name + Block, data = priority_table),
  Peduncle_senescence_DAH = lm(`Peduncle_senescence_DAH` ~ Line_name + Block, data = priority_table),
  Leaf_senescence = lm(`Leaf_senescence` ~ Line_name + Block, data = priority_table),
  Peduncle_senescence = lm(`Peduncle_senescence` ~ Line_name + Block, data = priority_table),
  Grain_mass = lm(`Grain_mass` ~ Line_name + Block, data = priority_table)
)

#Where to put the TUKEY HSD letters
label_yvals <- c(1000, 420, 45, 57, 42, 48, 84, 91, 4)

#Plot everything at once
for(i in 1:length(basic_plots)){
  print(i)
  par(mfrow=c(2,2))
  plot(linear_models[[i]], which = c(1,2))
  LABELS <- plot_anova_and_return_tukey_labels(linear_models[[i]])
  par(mfrow=c(1,1))
  print(basic_plots[[i]] +
    geom_boxplot() +
    genotype_scale +
    geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = label_yvals[i], fill = "white", size = 6, label.size = 0))
}

#Let's double check the Peduncle senescence and Grain mass where 'Block' is significant
LABELS <- plot_anova_and_return_tukey_labels(linear_models[[6]], explanatory = "Block")
ggplot(priority_table, aes(x=Block, y=`Peduncle_senescence_DAH`)) +
  geom_boxplot() +
  genotype_scale +
  geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = label_yvals[6], fill = "white", size = 6, label.size = 0)


LABELS <- plot_anova_and_return_tukey_labels(linear_models[[7]], explanatory = "Block")
ggplot(priority_table, aes(x=Block, y=`Grain_mass`)) +
  geom_boxplot() +
  genotype_scale +
  geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = label_yvals[7], fill = "white", size = 6, label.size = 0)

#Scatter plots
basic_plots_scatter <- list(
  Heading_date = ggplot(priority_table, aes(x=log(`CENAP-7F-2R`), y=Heading_date, colour = Parental_genotype)),
  Leaf_senescence = ggplot(priority_table, aes(x=log(`CENAP-7F-2R`), y=Leaf_senescence, colour = Parental_genotype)),
  Peduncle_senescence = ggplot(priority_table, aes(x=log(`CENAP-7F-2R`), y=Peduncle_senescence, colour = Parental_genotype)),
  Leaf_senescence_DAH = ggplot(priority_table, aes(x=log(`CENAP-7F-2R`), y=Leaf_senescence_DAH, colour = Parental_genotype)),
  Peduncle_senescence_DAH = ggplot(priority_table, aes(x=log(`CENAP-7F-2R`), y=Peduncle_senescence_DAH, colour = Parental_genotype)),
  Grain_mass = ggplot(priority_table, aes(x=log(`CENAP-7F-2R`), y=Grain_mass, colour = Parental_genotype))
)

for(i in 1:length(basic_plots_scatter)){
  print(basic_plots_scatter[[i]] +
          geom_point() +
          genotype_scale_colour
          )
}


for(i in 1:length(basic_plots_scatter)){
  print(basic_plots_scatter[[i]] +
          geom_point() +
          facet_wrap(facets = vars(Line_name)) +
          genotype_scale_colour
  )
}

#Pearson's correlations to go with the scatter plots
cor.test(priority_table$Heading_date, log(priority_table$`CENAP-7F-2R`), alternative = 'two.sided', method = "pearson")
cor.test(priority_table$Leaf_senescence_DAH, log(priority_table$`CENAP-7F-2R`), alternative = 'two.sided', method = "pearson")
cor.test(priority_table$Peduncle_senescence_DAH, log(priority_table$`CENAP-7F-2R`), alternative = 'two.sided', method = "pearson")
cor.test(priority_table$Grain_mass, log(priority_table$`CENAP-7F-2R`), alternative = 'two.sided', method = "pearson")


#another plot
ggplot(priority_table, aes(x=Heading_date, y=Leaf_senescence_DAH, colour = log(`CENAP-7F-2R`))) +
  geom_point()

