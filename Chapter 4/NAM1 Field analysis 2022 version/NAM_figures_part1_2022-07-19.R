# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Analysis of NAM data from field 2022 PART 1
#'
#'19/07/2022
#'Based on NAM_figures_2022-03-11.Rmd
#'NAM 2022 data
#'Leaf senescence stats only
#'No models - script was already too long
#'
#'Run from RAGT internal
#'Made in R 4.2.0
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#SOURCE
#Everything we need
setwd("U:/Field senescence re-analysis")
source("scripts/NAM_figures_source_2022-07-20.R")

source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_scripts/NAM2_thesis_graphs_source_2023-07-31.R', encoding = 'UTF-8', echo=TRUE)

setwd("U:/Field senescence re-analysis/Figures")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#PLOT

### Tile plots
for(i in 1:length(variable_list)){
  plot2 <- 
    Wide_data %>%
    ggplot(aes_string(x="X", y="Y", fill = variable_list[i])) +
    geom_tile() +
    scale_fill_gradient(low = "green4", high = "gold") +
    # facet_wrap(vars(Trial), nrow = 2, ncol = 1) +
    theme_minimal() +
    labs(fill = str_wrap(variable_labels[i], width = 20))
  print(plot2)
  save_small_svg(plot2, "tile_plot_2022", variable_list[i])
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Senescence
Facet_leaf_sen <- Long_data_no_control %>%
  ggplot(aes(x=Degree_days_after_heading, y = Leaf_senescence, group = DATA_ID, col = Genotype)) +
  geom_point(shape=1) +
  geom_line() +
  genotype_scale_col +
  facet_grid(rows = vars(Genotype), cols = vars(Rep)) +
  labs(col = genotype_lab, y = leaf_senescence_lab) +
  ylim(c(0,100))

Facet_leaf_sen
Facet_leaf_sen_2 <- Facet_leaf_sen +
  geom_point(data = Long_data_control, col = "black", aes(shape = Line_name)) +
  geom_line(data = Long_data_control, col = "black")

Facet_ear_sen <- Long_data_no_control %>%
  ggplot(aes(x=Degree_days_after_heading, y = Ear_senescence, group = DATA_ID, col = Genotype)) +
  geom_point(shape=1) +
  geom_line() +
  genotype_scale_col +
  facet_grid(rows = vars(Genotype), cols = vars(Rep)) +
  labs(col = genotype_lab, y = ear_senescence_lab) +
  ylim(c(0,100))

Facet_ear_sen
Facet_ear_sen_2 <- Facet_ear_sen +
  geom_point(data = Long_data_control, col = "black", aes(shape = Line_name)) +
  geom_line(data = Long_data_control, col = "black")
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

Averages_plot_leaf <- ggplot(Long_data_no_control, aes(x=Date, y=Leaf_senescence, group = Genotype, color=Genotype))
Averages_plot_leaf <- Averages_plot_leaf +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  genotype_scale_col +
  labs(color = genotype_lab, y = leaf_senescence_lab) +
  ylim(c(0,100))
Averages_plot_leaf

Averages_plot_leaf_2 <- Averages_plot_leaf +
  stat_summary(data = Long_data_control, aes(group = Line_name, shape = Line_name), fun.data = "mean_se", geom="point", size=3, col = "black") +
  stat_summary(data = Long_data_control, aes(group = Line_name), fun.data = "mean_se", geom="line", col = "black") +
  stat_summary(data = Long_data_control, aes(group = Line_name), fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  scale_shape_manual(values = c("RGT BOWEN"=15, "SKYFALL"=17))

Averages_plot_ear <- ggplot(Long_data_no_control, aes(x=Date, y=Ear_senescence, group = Genotype, color=Genotype))
Averages_plot_ear <- Averages_plot_ear +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  genotype_scale_col +
  labs(color = genotype_lab, y = ear_senescence_lab) +
  ylim(c(0,100))
Averages_plot_ear

Averages_plot_ear_2 <- Averages_plot_ear +
  stat_summary(data = Long_data_control, aes(group = Line_name, shape = Line_name), fun.data = "mean_se", geom="point", size=3, col = "black") +
  stat_summary(data = Long_data_control, aes(group = Line_name), fun.data = "mean_se", geom="line", col = "black") +
  stat_summary(data = Long_data_control, aes(group = Line_name), fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  scale_shape_manual(values = c("RGT BOWEN"=15, "SKYFALL"=17))

save_panel_svg(ggarrange(Averages_plot_leaf, Averages_plot_ear, ncol = 1), "Averages_plot_2022", "NAM1",
               n_panel_cols = 1, n_panel_rows = 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Correlation plot
Wide_data_no_control %>%
  dplyr::select(
                Leaf_TT10, Leaf_TT50, Leaf_TT90,
                Dur_Leaf_Sen,
                AUC_Leaf_Sen,
                Ear_TT10, Ear_TT50, Ear_TT90,
                Dur_Ear_Sen,
                AUC_Ear_Sen,
                `Leaf_curling.2022-06-22`,
                `Peduncle_senescence.2022-07-08`,
                `Peduncle_senescence.2022-07-12`
                ) %>%
  ggcorr_with_theme()

#Comparisons
#WITH + WITHOUT FACETS
#Switch to calculated metrics
#No controls because they didn't come with EE data

basic_plots <- list(
  ggplot(Wide_data_no_control, aes_string(x = "AUC_Leaf_Sen", y="Dur_Leaf_Sen", col = "Genotype")),
  ggplot(Wide_data_no_control, aes_string(x = "AUC_Leaf_Sen", y="AUC_Ear_Sen", col = "Genotype")),
  ggplot(Wide_data_no_control, aes(x = AUC_Leaf_Sen, y=`Leaf_curling.2022-06-22`, col = Genotype)),
  ggplot(Wide_data_no_control, aes(x = AUC_Leaf_Sen, y=`Peduncle_senescence.2022-07-12`, col = Genotype)),
  ggplot(Wide_data_no_control, aes(x =AUC_Ear_Sen, y=`Peduncle_senescence.2022-07-12`, col = Genotype)),
  ggplot(Wide_data_no_control, aes(x = `Ear_senescence.2022-07-12`, y=`Peduncle_senescence.2022-07-12`, col = Genotype))
  #ggplot(Wide_data_no_control, aes_string(x = "AUC_Leaf_Sen", y="Yield", col = "Genotype")),
  #ggplot(Wide_data_no_control, aes_string(x = "Dur_Leaf_Sen", y="Yield", col = "Genotype")),
  #ggplot(Wide_data_no_control, aes_string(x = "AUC_Ear_Sen", y="Yield", col = "Genotype")),
  #ggplot(filter(Wide_data_no_control, Rep == 1), aes_string(x = "Yield", y="Protein", col = "Genotype")),
  #ggplot(filter(Wide_data_no_control, Rep == 1), aes_string(x = "AUC_Leaf_Sen", y="Protein", col = "Genotype")),
  #ggplot(filter(Wide_data_no_control, Rep == 1), aes_string(x = "Dur_Leaf_Sen", y="Protein", col = "Genotype")),
  #ggplot(filter(Wide_data_no_control, Rep == 1), aes_string(x = "AUC_Ear_Sen", y="Protein", col = "Genotype"))
)
for(i in 1:length(basic_plots)){
  print(
    basic_plots[[i]] +
      geom_point() +
      genotype_scale_col +
      labs(col = genotype_lab) +
      stat_cor(method = "pearson", label.y = 18) +
      facet_grid(rows = vars(NAMA1), cols = vars(NAMB1))
  )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#SAVE PLOTS

setwd("./NAM_Figures")
ggsave(
  paste("NAM_Averages_plot_leaf_", lubridate::today(), ".png", sep = ""),
  Averages_plot_leaf,
  width = 160,
  height = 90,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Averages_plot_ear_", lubridate::today(), ".png", sep = ""),
  Averages_plot_ear,
  width = 160,
  height = 90,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Facet_leaf_sen_", lubridate::today(), ".png", sep = ""),
  Facet_leaf_sen,
  width = 320,
  height = 180,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Facet_ear_sen_", lubridate::today(), ".png", sep = ""),
  Facet_ear_sen,
  width = 320,
  height = 180,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Averages_plot_leaf_2_", lubridate::today(), ".png", sep = ""),
  Averages_plot_leaf_2,
  width = 160,
  height = 90,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Averages_plot_ear_2_", lubridate::today(), ".png", sep = ""),
  Averages_plot_ear_2,
  width = 160,
  height = 90,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Facet_leaf_sen_2_", lubridate::today(), ".png", sep = ""),
  Facet_leaf_sen_2,
  width = 320,
  height = 180,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM_Facet_ear_sen_2_", lubridate::today(), ".png", sep = ""),
  Facet_ear_sen_2,
  width = 320,
  height = 180,
  units = "mm",
  device = "png"
)
