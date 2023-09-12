# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Analysis of NAM data from field 2022 PART 2
#'
#'19/07/2022
#'Based on NAM_figures_2022-03-11.Rmd and Poster figures 2022-03-20.R
#'NAM 2022 data
#'Leaf senescence stats only
#'Fig. 1 Methods for calculating AUC and Dur - added direct to poster
#'Fig. 3 Boxplots - added direct to poster
#'
#'03/08/2023
#'Adapt for NAM2 Field 2023
#'
#'Run from U:
#'Made in R 4.2.0
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#SOURCE
#Themes n stuff

gene = "NAM2"
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_scripts/NAM2_thesis_graphs_source_2023-07-31.R', echo=TRUE)

#OPEN
local_file_path <- get_local_file_path()
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/2023-08-02_clean_data", sep = ""))
Long_data <- read_csv("Field_results_NAM2_long_2023-08-04.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

Wide_data <- read_csv("Field_results_NAM2_wide_DDAH_2023-08-04.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

#Tweak

Long_data <- Long_data %>%
  mutate(Rep = X)
Wide_data <- Wide_data %>%
  mutate(Rep = X)

setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-21 NAM2 Field 23/Figures", sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#PLOT

#FIG 1

#Not needed, can use figure I already have

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#OUTLIERS - there are none
outlier_plots <- unlist(Wide_data[which(is.na(Wide_data$Leaf_TT10)),"Index"])
Facet_leaf_sen +
  geom_point(data = filter(Long_data, Index %in% outlier_plots), col = "black")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#MODELS and BOXPLOTS

#Factorial model (optional)
linear_models_factorial <- list(
  AUC = lm(AUC_Leaf_Sen ~ Genotype, data = Wide_data),
  Dur = lm(Dur_Leaf_Sen ~ Genotype, data = Wide_data),
  TT10 = lm(Leaf_TT10 ~ Genotype, data = Wide_data),
  TT30 = lm(Leaf_TT30 ~ Genotype, data = Wide_data),
  TT50 = lm(Leaf_TT50 ~ Genotype, data = Wide_data),
  TT70 = lm(Leaf_TT70 ~ Genotype, data = Wide_data),
  TT90 = lm(Leaf_TT90 ~ Genotype, data = Wide_data),
  Peduncle = lm(`Peduncle_senescence.2023-07-20` ~ Genotype, data = Wide_data)
)

for(i in 1:length(linear_models_factorial)){
  print(names(linear_models_factorial)[i])
  par(mfrow=c(1,2))
  plot(linear_models_factorial[[i]], which = c(1,2))
  print(car::Anova(linear_models_factorial[[i]],type=3))
  display(linear_models_factorial[[i]])
}

print(paste(target, "means"))
estimated_means <- emmeans(linear_models_factorial$AUC, ~ NAMA1 + NAMB1, nesting=NULL)
print(estimated_means)
print(confint(pairs(estimated_means), level = 0.95))

print(paste(target,"Dur", "means"))
estimated_means <- emmeans(linear_models_factorial$Dur, ~ Genotype, nesting=NULL)
print(estimated_means)
print(confint(pairs(estimated_means), level = 0.95))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#FIG 3
#BOXPLOTS not BOXPLOTS 'cause there aren't enough datapoints
#SAVE SAVE SAVE
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/Figures", sep = ""))

#VARIABLES to test in thesis graphs source

#MODELS

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters
model_formula <- " ~ Cross*Genotype"
#Adding both Rep (block) and Cross would break it (no df left)
n_variables <- 2
explanatory <- "Cross*Genotype"

#Iterate through variables
  
#Plot everything at once
for(i in 1:length(variable_list_field)){
  print(variable_list_field[i])
  linear_model <- lm(formula(paste(variable_list_field[i], model_formula, sep = "")), data = filter(Wide_data, Cross %in% c("X50", "X61")))
  
  par(mfrow=c(1,2))
  # plot(linear_model, which = c(1,2))
  #Crucial ANOVA and "compact letter display" labels
  labels <- plot_anova_and_return_labels_emm(linear_model, explanatory = explanatory, alpha = 0.05, variables = n_variables)
  #axes
  ymin <- min(Wide_data[variable_list_unquote_field[i]], na.rm = TRUE)
  ymax <- max(Wide_data[variable_list_unquote_field[i]], na.rm = TRUE)
  label_y <- ymax + 0.12*(ymax - ymin)
  axis_ymin <- ymin
  axis_ymax <- ymax + 0.18*(ymax - ymin)
  #final boxplot
  par(mfrow=c(1,1))
  
  plot <- ggplot(filter(Wide_data, Cross %in% c("X50", "X61")),
                 aes_string(x = "Genotype", y = variable_list_field[i], fill = "Dosage")
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(height = 0)) +
    #stat_summary(fun.data = "mean_se", col = "black") +
    genotype_scale_x +
    dosage_scale_fill +
    labs(y = str_wrap(variable_labels_field[variable_list_unquote_field[i]], width = 23), x = NULL) +
    facet_wrap(vars(Cross), ncol = 2) +
    geom_label(
      data = labels,
      x = labels[, n_variables],
      label = labels[, (n_variables + 1)],
      y = label_y,
      fill = "white",
      col = "black",
      size = 16 / .pt, #size 12 for combined plots, size 16 - 20 for individual plots
      label.size = 0
    ) +
    coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
    my_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  print(plot)
  plot_list[[list_index]] <- plot
  list_index <- list_index + 1
  #save_A5_svg(plot, str_c("Boxplot", "linename", "dosage", variable_list_truncated_field[i], sep = "_"), gene = gene)
}


#COMBINE
j <- 1
combined <- ggarrange(plot_list[[j]],
                      plot_list[[j+1]],
                      plot_list[[j+2]],
                      plot_list[[j+3]],
                      ncol = 1,
                      nrow = 4
)
combined

save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated_field[j], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 4, ratio = 1)

j <- 5
combined <- ggarrange(plot_list[[j]],
                      plot_list[[j+1]],
                      plot_list[[j+2]],
                      plot_list[[j+4]],
                      ncol = 1,
                      nrow = 4
)
combined

save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated_field[j], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 4, ratio = 1)

j <- 10
combined <- ggarrange(plot_list[[j]],
                      plot_list[[j+1]],
                      ncol = 1,
                      nrow = 2
)
combined

save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated_field[j], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 2, ratio = 0.75)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#FIG 4
#BOXPLOTS for other discrete variables
#Set model parameters to the desired variable

#SAVE SAVE SAVE
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/Figures", sep = ""))

#VARIABLES to test in thesis graphs source

#MODELS

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters
model_formula <- " ~ `Vigour.2022-11-10`"
#Adding both Rep (block) and Cross would break it (no df left)
n_variables <- 1
explanatory <- "`Vigour.2022-11-10`"

#Iterate through variables

#Plot everything at once
for(i in 1:length(variable_list_field)){
  print(variable_list_field[i])
  linear_model <- lm(formula(paste(variable_list_field[i], model_formula, sep = "")), data = Wide_data)
  
  par(mfrow=c(1,2))
  # plot(linear_model, which = c(1,2))
  #Crucial ANOVA and "compact letter display" labels
  labels <- plot_anova_and_return_labels_emm(linear_model, explanatory = explanatory, alpha = 0.05, variables = n_variables)
  #axes
  ymin <- min(Wide_data[variable_list_unquote_field[i]], na.rm = TRUE)
  ymax <- max(Wide_data[variable_list_unquote_field[i]], na.rm = TRUE)
  label_y <- ymax + 0.12*(ymax - ymin)
  axis_ymin <- ymin
  axis_ymax <- ymax + 0.18*(ymax - ymin)
  #final boxplot
  par(mfrow=c(1,1))
  
  plot <- ggplot(Wide_data,
                 aes_string(x = explanatory, y = variable_list_field[i])
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(height = 0)) +
    labs(y = str_wrap(variable_labels_field[variable_list_unquote_field[i]], width = 23), x = NULL) +
    #facet_wrap(vars(Cross), ncol = 2) +
    geom_label(
      data = labels,
      x = labels[, n_variables],
      label = labels[, (n_variables + 1)],
      y = label_y,
      fill = "white",
      size = 16 / .pt, #size 12 for combined plots, size 16 - 20 for individual plots
      label.size = 0
    ) +
    coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
    my_theme + theme(legend.position = "none")
  
  print(plot)
  plot_list[[list_index]] <- plot
  list_index <- list_index + 1
  save_A5_svg(plot, str_c("Boxplot", explanatory, variable_list_truncated_field[i], sep = "_"), gene = gene)
}


#FACET
combined <- ggarrange(plot_list[[10]],
                      plot_list[[11]],
                      plot_list[[12]],
                      plot_list[[1]],
                      plot_list[[2]],
                      plot_list[[4]],
                      plot_list[[7]],
                      plot_list[[8]],
                      ncol = 2,
                      nrow = 4
)
combined

save_panel_svg(combined, paste("Combined_boxplots", explanatory, sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 4, ratio = 0.75)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#FIG 5
#Scatter plots
for(i in 1:length(variable_list_field)){
  for(j in 1:length(variable_list_field)){
    scatter <- ggplot(Wide_data,
           aes_string(x = variable_list_field[i], y = variable_list_field[j])) +
      geom_point() +
      stat_cor(method = "pearson", label.y = 18) +
      facet_wrap(vars(Cross), ncol = 2) 
    print(scatter)
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#SAVE
setwd("./NAM_Figures")

#FIG 3
