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
#'02/08/2023
#'Adapt for NAM2 data from field 2023
#'
#'03/08/2023
#'Use Degree_days_after_heading
#'
#'Run from U: drive
#'Made in R 4.2.0
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#SOURCE
#Themes n stuff

gene = "NAM2"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_scripts/NAM2_thesis_graphs_source_2023-07-31.R', echo=TRUE)

#OPEN
local_file_path <- get_local_file_path()
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/2023-08-02_clean_data", sep = ""))
Long_data <- read_csv("Field_results_NAM2_long_2023-08-03.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

Wide_data <- read_csv("Field_results_NAM2_wide_DDAH_2023-08-03.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

#Tweak

Long_data <- Long_data %>%
  mutate(Rep = X)
Wide_data <- Wide_data %>%
  mutate(Rep = X)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#PLOT

### Tile plots
for(i in 1:length(variable_list_field)){
  plot2 <- 
    Wide_data %>%
    ggplot(aes_string(x="X", y="Y", fill = variable_list_field[i])) +
    geom_tile() +
    scale_fill_gradient(low = "gold", high = "green4") +
    # facet_wrap(vars(Trial), nrow = 2, ncol = 1) +
    theme_minimal() +
    labs(fill = str_wrap(variable_labels_field[i], width = 20))
  print(plot2)
  save_A5_svg(plot2, paste("NAM_tileplot_", variable_list_field[i], "_", lubridate::today(), ".png", sep = ""), gene = gene)
}

#Senescence
Facet_leaf_sen <- Long_data %>%
  filter(Genotype %in% c("G1", "G8")) %>%
  ggplot(aes(x=Degree_days_after_heading, y = Leaf_senescence, group = Index, col = Dosage)) +
  geom_point(shape=1) +
  geom_line() +
  dosage_scale_col +
  #facet_grid(rows = vars(`Vigour.2022-11-10`), cols = vars(Cross)) +
  facet_wrap(vars(Cross)) +
  labs(col = genotype_lab, y = SPAD_lab, x = DDAH_lab) +
  ylim(c(0,100))

Facet_leaf_sen + theme(legend.position = "none")

save_A5_svg(Facet_leaf_sen + theme(legend.position = "none"), "Facet_leaf_sen", gene = gene)

#Peduncle
Facet_ear_sen <- Long_data %>%
  filter(Genotype %in% c("G1", "G8")) %>%
  ggplot(aes(x=Degree_days_after_heading, y = Peduncle_senescence, group = Index, col = Dosage)) +
  geom_point(shape=1) +
  geom_line() +
  dosage_scale_col +
  #facet_grid(rows = vars(Genotype), cols = vars(Rep)) +
  facet_wrap(vars(Cross)) +
  labs(col = genotype_lab, y = peduncle_senescence_lab) +
  ylim(c(0,10))

Facet_ear_sen + theme(legend.position = "none")

save_A5_svg(Facet_ear_sen + theme(legend.position = "none"), "Facet_ear_sen", gene = gene)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Have to use date here for sensible averaging
Averages_plot_leaf <- ggplot(Long_data, aes(x=Date, y=Leaf_senescence, group = Genotype, color=Genotype))
Averages_plot_leaf <- Averages_plot_leaf +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  genotype_scale_col +
  labs(color = genotype_lab, y = leaf_senescence_lab) +
  ylim(c(0,100))
Averages_plot_leaf

Averages_plot_ear <- ggplot(Long_data, aes(x=Date, y=Peduncle_senescence, group = Genotype, color=Genotype))
Averages_plot_ear <- Averages_plot_ear +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  genotype_scale_col +
  labs(color = genotype_lab, y = peduncle_senescence_lab) +
  ylim(c(0,10))
Averages_plot_ear

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Correlation plot
Wide_data %>%
  dplyr::select(
    Leaf_TT10, Leaf_TT50, Leaf_TT90,
    Dur_Leaf_Sen,
    AUC_Leaf_Sen
  ) %>%
  ggcorr_with_theme()

#Comparisons
#WITH + WITHOUT FACETS
#Switch to calculated metrics
#No controls because they didn't come with EE data

basic_plots <- list(
  ggplot(Wide_data, aes_string(x = "AUC_Leaf_Sen", y="Dur_Leaf_Sen", col = "Genotype")),
  ggplot(Wide_data, aes(x = AUC_Leaf_Sen, y=`Peduncle_senescence.2023-07-20`, col = Genotype))#,
  #ggplot(Wide_data, aes_string(x = "AUC_Leaf_Sen", y="Yield", col = "Genotype")),
  #ggplot(Wide_data, aes_string(x = "Dur_Leaf_Sen", y="Yield", col = "Genotype")),
  #ggplot(Wide_data, aes_string(x = "AUC_Ear_Sen", y="Yield", col = "Genotype")),
  #ggplot(filter(Wide_data, Rep == 1), aes_string(x = "Yield", y="Protein", col = "Genotype")),
  #ggplot(filter(Wide_data, Rep == 1), aes_string(x = "AUC_Leaf_Sen", y="Protein", col = "Genotype")),
  #ggplot(filter(Wide_data, Rep == 1), aes_string(x = "Dur_Leaf_Sen", y="Protein", col = "Genotype")),
  #ggplot(filter(Wide_data, Rep == 1), aes_string(x = "AUC_Ear_Sen", y="Protein", col = "Genotype"))
)
for(i in 1:length(basic_plots)){
  print(
    basic_plots[[i]] +
      geom_point() +
      genotype_scale_col +
      labs(col = genotype_lab) +
      stat_cor(method = "pearson", label.y = 18) +
      #facet_grid(rows = vars(NAMA1), cols = vars(NAMB1))
      facet_wrap(vars(Cross))
  )
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#SAVE PLOTS

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2022-09-21 NAM2 Field 23/Figures")
ggsave(
  paste("NAM2_Averages_plot_leaf_", lubridate::today(), ".png", sep = ""),
  Averages_plot_leaf,
  width = 160,
  height = 90,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM2_Averages_plot_ear_", lubridate::today(), ".png", sep = ""),
  Averages_plot_ear,
  width = 160,
  height = 90,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM2_Facet_leaf_sen_", lubridate::today(), ".png", sep = ""),
  Facet_leaf_sen,
  width = 320,
  height = 180,
  units = "mm",
  device = "png"
)

ggsave(
  paste("NAM2_Facet_ear_sen_", lubridate::today(), ".png", sep = ""),
  Facet_ear_sen,
  width = 320,
  height = 180,
  units = "mm",
  device = "png"
)
