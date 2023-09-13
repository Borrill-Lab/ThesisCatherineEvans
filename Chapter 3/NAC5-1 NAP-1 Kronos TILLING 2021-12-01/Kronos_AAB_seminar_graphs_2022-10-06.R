#'Kronos seminar graphs
#'
#'26/09/2022
#'Build graphs for Kronos data for Crop Gen Seminar on 03/10/2022
#'Based on Kronos_exploratory_2022-09-07
#'
#'06/10/2022
#'Adapt graphs for ASM poster
#'new sizes
#'
#'26/10/2022
#'Adapt graphs for AAB seminar
#'Best of both Crop Gen seminar and ASM poster
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5 (wait no, this computer only has 4.0.4)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'Differences from analysis in Kronos_exploratory_2022-09-07
#'


#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date the scripts were run through Kronos_data_preparation.R
gene = "NAC5"
gene = "NAP"
date = "2022-09-07"
#
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/Kronos_poster_graphs_source_2022-10-06.R', echo=TRUE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Kronos_SPAD_analysis <- Kronos_SPAD_analysis %>%
  mutate(Day_approx = Week*7)

#PART 2 plot a SPAD plot
Averages_plot <- ggplot(filter(Kronos_SPAD_analysis, Genotype %in% c("X:X.X:X","Y:Y.Y:Y")), aes(x=Week, y=SPAD, group = Genotype, color=Genotype, shape = Genotype))
Averages_plot <- Averages_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  scale_x_continuous(breaks = seq(0, 9, 1), labels = ~(.x*7)) +
  genotype_scale_col_2 +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 30), x = "Days after Heading") +
  my_theme + theme()
Averages_plot + scale_x_continuous(breaks = seq(0, 9, 1), labels = ~(.x*7))

# Explore testing individual timepoints with a Kruskal-Wallis non-parametric test
compare <- ggpubr::compare_means(SPAD ~ Light,
                      data = filter(Kronos_SPAD_analysis, Genotype %in% c("X:X.X:X","Y:Y.Y:Y")),
                      group.by = "Week")
# Averages_plot + stat_compare_means(aes(group = Genotype), label =  "p.signif", label.y = 30)

Averages_plot_wilcox <- ggline(data = filter(Kronos_SPAD_analysis,  Genotype %in% c("X:X.X:X","Y:Y.Y:Y")), 
       x = "Day_approx", y = "SPAD", add = "mean_se",
       color = "Genotype", shape = "Genotype") +
  stat_compare_means(aes(group = Genotype), label = "p.signif",
                     label.y = 58, method = "wilcox.test", size = 20 / .pt) +
  genotype_scale_col_2 +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Days after Heading") +
  my_theme +
  theme(legend.position = "none")
Averages_plot_wilcox

Facet_leaf_sen <- Kronos_SPAD_yday %>%
  filter(!is.na(Plant_num) & !is.na(SPAD) & Genotype %in% c("X:X.X:X","Y:Y.Y:Y")) %>%
  mutate(Outlier = Index %in% Kronos_outliers$Index) %>%
  ggplot(aes(x=Date, y = SPAD, group = Plant_num, col = Genotype)) +
  geom_point(shape=1) +
  geom_line() +
  facet_wrap(facets = vars(Light)) +
  genotype_scale_col
Facet_leaf_sen

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#PART 3 some boxplots
Kronos_analysis <- Kronos_wide_yday %>%
  mutate(SPAD_initial = (SPAD_heading + SPAD_1 + SPAD_2 + SPAD_3)/4) %>%
  anti_join(Kronos_outliers, by = "Index") %>%
  filter(Genotype %in% c("X:X.X:X", "Y:Y.Y:Y", "X:X.Y:Y", "Y:Y.X:X"))

#VARIABLES to test
variable_list <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen", "Peduncle_senescence_DAH", "Heading_date", "Height", "Tiller_1",
                   "Seed_num",    "Weight",      "TGW", "O_length",
                   "`Predicted Moisture %`",          "`Predicted Protein Dry basis %`",
                   "`Predicted Starch Dry basis %`", "SPAD_initial")
variable_list_unquote <- c("Leaf_senescence_DAH", "TT_SPAD30","Dur_Leaf_Sen", "AUC_Leaf_Sen", "Peduncle_senescence_DAH", "Heading_date", "Height", "Tiller_1",
                           "Seed_num",    "Weight",      "TGW", "O_length",
                           "Predicted Moisture %",          "Predicted Protein Dry basis %",
                           "Predicted Starch Dry basis %", "SPAD_initial")
variable_list_truncated <- c(str_trunc(variable_list_unquote, width = 12, side = "right", ellipsis = "-"))

#MODELS
linear_models <- list(
  Leaf_senescence_DAH = lm(Leaf_senescence_DAH ~ Light + Block + Genotype, data = Kronos_analysis),
  TT_SPAD30 = lm(TT_SPAD30 ~ Light + Block + Genotype, data = Kronos_analysis),
  Dur_Leaf_Sen = lm(Dur_Leaf_Sen ~ Light + Block + Genotype, data = Kronos_analysis),
  AUC_Leaf_Sen = lm(AUC_Leaf_Sen ~ Light + Block + Genotype, data = Kronos_analysis),
  Peduncle_senescence_DAH = lm(Peduncle_senescence_DAH ~ Light + Block + Genotype, data = Kronos_analysis),
  Heading_date = lm(Heading_date ~ Light + Block + Genotype, data = Kronos_analysis),
  Height = lm(Height ~ Light + Block + Genotype, data = Kronos_analysis),
  Tiller_1 = lm(Tiller_1 ~ Light + Block + Genotype, data = Kronos_analysis),
  Seed_num = lm(Seed_num ~ Light + Block + Genotype, data = Kronos_analysis),
  Weight = lm(Weight ~ Light + Block + Genotype, data = Kronos_analysis),
  TGW = lm(TGW ~ Light + Block + Genotype, data = Kronos_analysis),
  O_length = lm(O_length ~ Light + Block + Genotype, data = Kronos_analysis),
  M = lm(`Predicted Moisture %` ~ Light + Block + Genotype, data = Kronos_analysis),
  P = lm(`Predicted Protein Dry basis %` ~ Light + Block + Genotype, data = Kronos_analysis),
  S = lm(`Predicted Starch Dry basis %` ~ Light + Block + Genotype, data = Kronos_analysis),
  SPAD_initial = lm(SPAD_initial ~ as.factor(Light) + Block + Genotype, data = Kronos_analysis)
)


# #Where to put the TUKEY HSD letters
# label_yvals <- c(1000, 420, 45, 57, 42, 48, 84, 91, 4)

#Plot everything at once
for(i in 1:length(variable_list)){
  print(names(linear_models)[i])
  par(mfrow=c(1,2))
  plot(linear_models[[i]], which = c(1,2))
  #Crucial ANOVA and "compact letter display" labels
  labels <- plot_anova_and_return_labels_emm(linear_models[[i]], explanatory = "Genotype", alpha = 0.05)
  #axes
  ymin <- min(Kronos_analysis[variable_list_unquote[i]], na.rm = TRUE)
  ymax <- max(Kronos_analysis[variable_list_unquote[i]], na.rm = TRUE)
  label_y <- ymax + 0.12*(ymax - ymin)
  axis_ymin <- ymin
  axis_ymax <- ymax + 0.18*(ymax - ymin)
  #final boxplot
  par(mfrow=c(1,1))
  plot <-
    ggplot(Kronos_analysis, 
           aes_string(x = "Genotype", y = variable_list[i], fill = "Genotype")
    ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(height = 0)) +
    # genotype_scale_x +
    # genotype_scale_fill +
    labs(y = str_wrap(variable_labels[variable_list_unquote[i]], width = 23), x = NULL, fill = genotype_lab) +
    geom_label(
      data = labels,
      x = labels[, 1],
      label = labels[, 2],
      y = label_y,
      fill = "white",
      size = 20 / .pt,
      label.size = 0
    ) +
    coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
    genotype_scale_fill + genotype_scale_x +
    my_theme +
    theme(legend.position = "none")
  print(plot)
  save_small_svg(plot, str_c("Boxplot", variable_list_truncated[i], sep = "_"), gene = gene)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#SAVE graphs
setwd("../Figures")
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

save_poster_pdf_1 <- function(plot, plot_name, gene){
  ggsave(
    filename = paste(gene, "_", plot_name, "_poster", lubridate::today(), ".pdf", sep = ""),
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

#SAVE SAVE SAVE
setwd("../Figures")
save_poster_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox", gene = gene)
save_small_svg(Averages_plot_wilcox, "Averages_plot_wilcox", gene = gene)
