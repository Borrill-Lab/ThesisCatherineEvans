#Demonstrate that shelves are different

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(car) #Anova for non-orthogonal data
require(emmeans) #marginal means estimates for when p-values are not enough
require(arm)
require(lubridate)
require(tidyverse) #Load tidyverse last to avoid masking
require(ggpubr) #for ggarrange
theme_set(theme_bw())

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-03 NAC phenotyping")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-03 NAC phenotyping")
source("Functions_for_phenotype_data.R")

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Overexpression CER/2022-08-12_clean_data")
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER/2022-08-12_clean_data")

#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date the scripts were run through Overexpression_data_preparation.R
gene = "NAC5"
date = "2022-09-06"

Overexpression_SPAD_yday <- read_csv(paste("Overexpression_SPAD_yday_", gene, "_", date, ".csv", sep= ""))
Overexpression_phenotyping_yday <- read_csv(paste("Overexpression_phenotyping_yday_", gene, "_", date, ".csv", sep= ""))

#LABELS
Line_name_lab <- paste(gene, "Line_name")
leaf_senescence_lab <- "Days to 25% leaf senescence"
peduncle_senescence_lab <- "Days to 100% peduncle senescence"
SPAD_lab <- "Flag leaf chlorophyll content (SPAD)"
Line_name_limits <- c(
  "Non functional mutant.Wild type",
  "Wild type.Wild type",
  "Non functional mutant.Functional mutant",
  "Wild type.Functional mutant"
)
Line_name_names <- c(
  `Non functional mutant.Wild type` = "aabbDD",
  `Wild type.Wild type` = "AAbbDD (WT)",
  `Non functional mutant.Functional mutant` = "aaBBDD",
  `Wild type.Functional mutant` = "AABBDD"
)

#Scales for ggplot2
Line_name_scale_fill <- scale_fill_brewer(type = "div", palette = "PRGn",
                                          limits = Line_name_limits, labels = Line_name_names)
Line_name_scale_col <- scale_colour_brewer(type = "div", palette = "PRGn",
                                           limits = Line_name_limits, labels = Line_name_names)
Line_name_scale_x <- scale_x_discrete(limits = Line_name_limits, labels = Line_name_names)
date_scale <- scale_x_date(date_breaks = "2 days")
percent_scale <- scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 20))
week_scale <- scale_x_continuous(breaks = seq(0,9,1))


# Overexpression NAC-5 phenotyping

## Outliers

print(paste("Gene:", gene, ", Date:", date))

# summary(Overexpression_wide_yday)

Overexpression_wide_yday <- Overexpression_phenotyping_yday %>%
  mutate(Shelf = as.factor(paste("Shelf", ((Block -1) %/% 5) + 1)),
         Block = as.factor(Block),
         Mildew = grepl("mildew", Comment))

Overexpression_SPAD_yday <- Overexpression_SPAD_yday %>%
  mutate(Shelf = as.factor(paste("Shelf", ((Block -1) %/% 5) + 1)),
         Block = as.factor(Block),
         Mildew = grepl("mildew", Comment))

Overexpression_outliers <- Overexpression_wide_yday %>%
  filter(is.na(Plant_num)|grepl("damage", Comment)|grepl("tagged", Comment)|snapped_pre_ped_sen)
#grepl("mildew", Comment)
# View(Overexpression_outliers)

Overexpression_analysis <- Overexpression_wide_yday %>%
  anti_join(Overexpression_outliers, by = "Index")

Overexpression_SPAD_analysis <- Overexpression_SPAD_yday %>%
  anti_join(Overexpression_outliers, by = "Index")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Spatial

#VARIABLES to test
variable_list <- c("Heading_date", "Tiller_1", "Seed_num")
variable_label <- c("Heading date", "Tiller number", "Seed number")
label_shelf <- c("High-spec lights 300umol", "Standard-spec lights 300umol")

plots <- list(0,0,0,0)

for(i in 1:length(variable_list)){
  plot2 <-
    ggplot(Overexpression_wide_yday, aes_string(x="Shelf", y = variable_list[i], group = "Shelf", color="Shelf")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(height = 0, width = 0.3)) +
    #scale_x_discrete(labels = str_wrap(label_shelf, width = 15)) +
    scale_color_manual(values = c(`Shelf 1`="#08519c", `Shelf 2` = "#6baed6")) +
    labs(y = variable_label[i]) +
    theme_minimal() +
    theme(legend.position = "none")
  plots[[i]] <- plot2
}

Averages_plot <- ggplot(Overexpression_SPAD_yday, aes(x=Week, y=SPAD, group = Shelf, color=Shelf))
Averages_plot <- Averages_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  week_scale +
  scale_color_manual(values = c(`Shelf 1`="#08519c", `Shelf 2` = "#6baed6")) +
  labs(y = "Flag leaf chlorophyll content (SPAD)") #+
  #scale_color_discrete(labels = str_wrap(label_shelf, width = 15))
plots[[4]] <- Averages_plot

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]])

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]]) %>% annotate_figure(top = "Wheat grown in Extended height CER, Building 18")
