#'Overexpression seminar graphs SOURCE
#'
#'28/09/2022
#'Build graphs for Overexpression data for Crop Gen Seminar on 03/10/2022
#'Based on Overexpression_exploratory_2022-09-06
#'
#'29/09/2022
#'Use 2022-09-29 data with tidy Marvin data
#'
#'http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
#'
#'22/06/2022
#'Adapt for graphs to go in thesis for Overexpression CER data
#'Include Flag leaf length
#'Use Marvin adjusted lengths (see Marvin_and_NIR_data_reliability)
#'New graph saving protocols
#'14pt tukey letters
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date the scripts were run through Overexpression_data_preparation.R
gene = "NAC5"
#gene = "NAP"
date = "2022-09-29"

source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER/Overexpression_thesis_graphs_source_2023-06-22.R', echo=TRUE)
source('~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/NAC transgenics/2021-11-10 NAC5 NAP overexpression CER/2022-08-12_scripts/Overexpression_thesis_graphs_source_2023_06_22.R', echo=TRUE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#PART 1 focus on the single-copy transgenic lines with best overexpression
#KEEP ALL THE DATA
Overexpression_SPAD_analysis <- Overexpression_SPAD_analysis %>%
  # filter(Line_name %in% c("5.4_null", "5.4_hom",
  #                         "10.5_null", "10.5_hom", "10.9_null", "10.9_hom")) %>%
  tidyr::extract(Line_name, "Background", "(.*)_", remove = FALSE)

Overexpression_analysis <- Overexpression_analysis %>%
  # filter(Line_name %in% c("5.4_null", "5.4_hom", 
  #                         "10.5_null", "10.5_hom", "10.9_null", "10.9_hom")) %>%
  tidyr::extract(Line_name, "Background", "(.*)_", remove = FALSE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#PART 2 plot a SPAD plot
#FIGURE 9 @ 22/06/2023

Averages_plot <- ggplot(filter(Overexpression_SPAD_analysis, Genotype %in% c("null", "hom")), aes(x=Week, y=SPAD, group = Line_name, color=Genotype, shape = Genotype))
Averages_plot <- Averages_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  week_scale +
  genotype_scale_col +
  genotype_scale_shape +
  spad_scale +
  labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = str_wrap(variable_labels["Week"], width = 25)) +
  facet_grid(vars(Background), vars(Shelf)) +
  theme(legend.position = "top")


# Explore testing individual timepoints with a Kruskal-Wallis non-parametric test
#This thing doesn't work very well.
for(background in unique(Overexpression_SPAD_analysis$Background)){
  for(shelf in c("Shelf 1", "Shelf 2")){
    data = filter(Overexpression_SPAD_analysis, Background == background & Shelf == shelf)
    print(nrow(data))
    compare <- ggpubr::compare_means(SPAD ~ Genotype,
                      data = data,
                      group.by = "Week", method = "wilcox.test")
    print(compare)
  }
}

#For NAC-5, only 5.4 week 4 is significant. Add an asterisk afterwards

#SAVE this plot
save_A4_svg(Averages_plot, "Fig9_Averages_plot_Overexpression", gene = gene)

save_panel_svg(Averages_plot, "Fig9_Averages_plot_Overexpression", gene = gene,
               n_panel_cols = 2, n_panel_rows = 3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#PART 3 some boxplots
#Fig 10 @22/6/2023 senescence traits
#Fig 11 @22/6/2023 agronomic traits
#Fig 12 @22/6/2023 grain traits (leave space for Perl Stain images)

#VARIABLES to test
variable_list <- c("Leaf_senescence_DAH", "TT_SPAD30", "Dur_Leaf_Sen", "AUC_Leaf_Sen", "Peduncle_senescence_DAH", "Heading_date", "Tiller_1",
                   "Seed_num",    "Weight",      "TGW", "O_length_adj", "Flag_leaf_length")
variable_list_truncated <- c(str_trunc(variable_list, width = 12, side = "right", ellipsis = "-"))

#MODELS
linear_models <- list(
  Leaf_senescence_DAH = lm(Leaf_senescence_DAH ~ Shelf * Line_name, data = Overexpression_analysis),
  TT_SPAD30 = lm(TT_SPAD30 ~ Shelf * Line_name, data = Overexpression_analysis),
  Dur_Leaf_Sen = lm(Dur_Leaf_Sen ~ Shelf * Line_name, data = Overexpression_analysis),
  AUC_Leaf_Sen = lm(AUC_Leaf_Sen ~ Shelf * Line_name, data = Overexpression_analysis),
  Peduncle_senescence_DAH = lm(Peduncle_senescence_DAH ~ Shelf * Line_name, data = Overexpression_analysis),
  Heading_date = lm(Heading_date ~ Shelf * Line_name, data = Overexpression_analysis),
  Tiller_1 = lm(Tiller_1 ~ Shelf * Line_name, data = Overexpression_analysis),
  Seed_num = lm(Seed_num ~ Shelf * Line_name, data = Overexpression_analysis),
  Weight = lm(Weight ~ Shelf * Line_name, data = Overexpression_analysis),
  TGW = lm(TGW ~ Shelf * Line_name, data = Overexpression_analysis),
  O_length = lm(O_length ~ Shelf * Line_name, data = Overexpression_analysis),
  Flag_leaf_length = lm(Flag_leaf_length ~ Shelf * Line_name, data = Overexpression_analysis)
)


# #Where to put the TUKEY HSD letters
# label_yvals <- c(1000, 420, 45, 57, 42, 48, 84, 91, 4)
plotlist <- list()

#Plot everything at once
for(i in 1:length(variable_list)){
  print(names(linear_models)[i])
  par(mfrow=c(1,2))
  plot(linear_models[[i]], which = c(1,2))
  #Crucial ANOVA and "compact letter display" labels
  labels <- plot_anova_and_return_labels_emm(linear_models[[i]], explanatory = "Shelf*Line_name", alpha = 0.05, variables = 2)
  #axes
  ymin <- min(Overexpression_analysis[variable_list[i]], na.rm = TRUE)
  ymax <- max(Overexpression_analysis[variable_list[i]], na.rm = TRUE)
  label_y <- ymax + 0.12*(ymax - ymin)
  axis_ymin <- ymin
  axis_ymax <- ymax + 0.18*(ymax - ymin)
  #final boxplot
  par(mfrow=c(1,1))
  plot <-
    ggplot(Overexpression_analysis, 
           aes_string(x = "Line_name", y = variable_list[i], fill = "Genotype")
    ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(height = 0, width= 0.3)) +
    scale_x_discrete(limits = line_name_limits) +
    genotype_scale_fill_multi +
    labs(x = NULL, fill = Line_name_lab, y = str_wrap(variable_labels[variable_list[i]], width = 25)) +
    coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
    facet_wrap(vars(Shelf)) +
    my_theme +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90))
  print(plot)
  plot2 <- plot +
    #labels[, 1] is Shelf
    geom_label(
      data = labels,
      x = labels[, 2], #labels[,2] is Line name
      label = labels[, 3], #labels[,3] is emmeans::cld() letters
      y = label_y,
      fill = "white",
      size = 14 / .pt,
      label.size = 0
    ) #label doesn't work with facets
  print(plot2)
  plotlist[[i]] <- plot2
  save_panel_svg(plot2, str_c("Boxplot_facet_letters_all", variable_list_truncated[i], sep = "_"), gene = gene,
                 n_panel_cols = 2, n_panel_rows = 1)
}

setwd("./Figures/")
save_poster_pdf_1(Averages_plot_wilcox1, "Averages_plot_wilcox1", gene = gene)
save_poster_pdf_1(Averages_plot_wilcox2, "Averages_plot_wilcox2", gene = gene)

save_poster_pdf_1(Averages_plot_wilcox2 + theme(legend.position = "right"), "Averages_plot_wilcox2_legend", gene = gene)
