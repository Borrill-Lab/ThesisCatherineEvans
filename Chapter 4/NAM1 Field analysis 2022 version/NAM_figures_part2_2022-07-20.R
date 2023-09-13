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
#'23/08/2023
#'Degree days, thesis theme
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

year = 2022
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#PLOT

#FIG 1
#Method for calculating AUC and Dur
#TEST
#Work with an example set of 1 plot (picking a nice looking one)
Example <- Long_data %>%
  filter(DATA_ID == 2755505) %>%
  dplyr::select(DATA_ID, Degree_days_after_heading, Leaf_senescence, Ear_senescence)

# #Here's what it looks like
# qplot(Example$Degree_days_after_heading, Example$Leaf_senescence, geom = c("point","line"))

#Use "approx" function to do linear interpolation
#Data MUST be in Degree_days_after_heading order

estimate <- approx(Example$Leaf_senescence, Example$Degree_days_after_heading, c(10, 30, 50, 70, 90), ties = "ordered")
estimate <- as_tibble(estimate) %>%
  rename(Leaf_senescence = x,
         Degree_days_after_heading = y) #Here this uses Date, NOT POSIXct


#Demo how I've made the stats

TT10_x <- filter(estimate, Leaf_senescence == 10)$Degree_days_after_heading
TT90_x <- filter(estimate, Leaf_senescence == 90)$Degree_days_after_heading

methods_figure <- ggplot(Example, aes(x = Degree_days_after_heading, y = Leaf_senescence)) +
  #AUC_Leaf_Sen
  geom_area(fill = "lightblue1") +
  geom_text(aes(x=700, y=30, label = "AUC_Leaf_Sen"), 
            col = "blue", size = 8) +
  #TT10 and TT90
  geom_point(data = filter(estimate, Leaf_senescence %in% c(10, 90)), 
             col = "red", size = 4) +
  geom_text(data = filter(estimate, Leaf_senescence %in% c(10, 90)), 
            aes(label = paste(Leaf_senescence, "%")),
            col = "red", nudge_x = 2, nudge_y = -4, size = 8) +
  #Dur_Leaf_Sen
  geom_line(data = filter(estimate, Leaf_senescence %in% c(10, 90)), 
            aes(x = Degree_days_after_heading, y = c(100, 100)),
            col = "red", size = 2) +
  geom_line(data = data.frame(x = c(TT10_x, TT10_x), y = c(10, 100)),
            aes(x = x, y = y),
            col = "red", linetype = 2, size = 1.2) +
  geom_line(data = data.frame(x = c(TT90_x, TT90_x), y = c(90, 100)),
            aes(x = x, y = y),
            col = "red", linetype = 2, size = 1.2) +
  geom_text(aes(x = 600, y = 90, label = c("Dur_Leaf_Sen")), 
            col = c("red"), size = 8) +
  #Actual plot
  geom_point(size = 2) +
  geom_line(size = 1.2) +
  #Use str_wrap to deal with long x-axis label
  percent_scale +
  #xlim(50, 76) +
  labs(y = leaf_senescence_lab, x = str_wrap(da1m_lab, width = 20)) +
  theme(plot.margin = unit(c(5,5,5,5), units = "mm"))

methods_figure

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#OUTLIERS - there are none
outlier_plots <- unlist(Wide_data[which(is.na(Wide_data$Leaf_TT10)),"DATA_ID"])
Facet_leaf_sen +
  geom_point(data = filter(Long_data_no_control, DATA_ID %in% outlier_plots), col = "black")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#MODELS and BOXPLOTS
target = "Leaf senescence"
print(target)

basic_plots <- list(
  AUC_Leaf_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=AUC_Leaf_Sen, fill = Genotype)),
  Dur_Leaf_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=Dur_Leaf_Sen, fill = Genotype)),
  Leaf_TT10 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Leaf_TT10, fill = Genotype)),
  Leaf_TT30 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Leaf_TT30, fill = Genotype)),
  Leaf_TT50 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Leaf_TT50, fill = Genotype)),
  Leaf_TT70 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Leaf_TT70, fill = Genotype)),
  Leaf_TT90 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Leaf_TT90, fill = Genotype)),
  Leaf_curling = ggplot(Wide_data_no_control, aes(x=Genotype, y=`Leaf_curling.2022-06-22`, fill = Genotype)),
  Peduncle = ggplot(Wide_data_no_control, aes(x=Genotype, y=`Peduncle_senescence.2022-07-08`, fill = Genotype))
)

linear_models_factorial <- list(
  AUC = lm(AUC_Leaf_Sen ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  Dur = lm(Dur_Leaf_Sen ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT10 = lm(Leaf_TT10 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT30 = lm(Leaf_TT30 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT50 = lm(Leaf_TT50 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT70 = lm(Leaf_TT70 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT90 = lm(Leaf_TT90 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  Leaf_curling = lm(`Leaf_curling.2022-06-22` ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  Peduncle = lm(`Peduncle_senescence.2022-07-08` ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  Height = lm(Height ~ NAMA1 * NAMB1, data = Wide_data_no_control)
)

linear_models <- list(
  AUC = lm(AUC_Leaf_Sen ~ Genotype, data = Wide_data_no_control),
  Dur = lm(Dur_Leaf_Sen ~ Genotype, data = Wide_data_no_control),
  TT10 = lm(Leaf_TT10 ~ Genotype, data = Wide_data_no_control),
  TT30 = lm(Leaf_TT30 ~ Genotype, data = Wide_data_no_control),
  TT50 = lm(Leaf_TT50 ~ Genotype, data = Wide_data_no_control),
  TT70 = lm(Leaf_TT70 ~ Genotype, data = Wide_data_no_control),
  TT90 = lm(Leaf_TT90 ~ Genotype, data = Wide_data_no_control),
  Leaf_curling = lm(`Leaf_curling.2022-06-22` ~ Genotype, data = Wide_data_no_control),
  Peduncle = lm(`Peduncle_senescence.2022-07-08` ~ Genotype, data = Wide_data_no_control)
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
estimated_means <- emmeans(linear_models_factorial$Dur, ~ NAMA1 * NAMB1, nesting=NULL)
print(estimated_means)
print(confint(pairs(estimated_means), level = 0.95))

label_yvals <- c(1700,20,65,70,70,71,74, 55, 6.5) #Height for Tukey Letters


for(i in 1:length(basic_plots)){
  par(mfrow=c(1,1))
  LABELS <- plot_anova_and_return_tukey_labels(linear_models[[i]], explanatory = "Genotype")
  formatted_plot <- basic_plots[[i]] +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.3) +
    genotype_scale_x +
    genotype_scale_fill +
    scale_y_continuous(expand = expansion(mult = 0.2, add = 0)) +
    labs(x = "", fill = genotype_lab) +
    geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = label_yvals[i], fill = "white", size = 6, label.size = 0) +
    #update y axis labels
    labs(y = str_wrap(variable_labels_leaf[i], width = 30))+
    theme(legend.position = "none")
  print(formatted_plot)
  save_small_svg(formatted_plot, variable_list[i], "NAM1_Boxplot")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Ear senescence
target = "Ear senescence"
print(target)
basic_plots <- list(
  AUC_Ear_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=AUC_Ear_Sen, fill = Genotype)),
  Dur_Ear_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=Dur_Ear_Sen, fill = Genotype)),
  Ear_TT10 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Ear_TT10, fill = Genotype)),
  Ear_TT30 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Ear_TT30, fill = Genotype)),
  Ear_TT50 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Ear_TT50, fill = Genotype)),
  Ear_TT70 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Ear_TT70, fill = Genotype)),
  Ear_TT90 = ggplot(Wide_data_no_control, aes(x=Genotype, y=Ear_TT90, fill = Genotype))
)
linear_models_factorial <- list(
  AUC = lm(AUC_Ear_Sen ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  Dur = lm(Dur_Ear_Sen ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT10 = lm(Ear_TT10 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT30 = lm(Ear_TT30 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT50 = lm(Ear_TT50 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT70 = lm(Ear_TT70 ~ NAMA1 * NAMB1, data = Wide_data_no_control),
  TT90 = lm(Ear_TT90 ~ NAMA1 * NAMB1, data = Wide_data_no_control)
)
linear_models <- list(
  AUC = lm(AUC_Ear_Sen ~ Genotype, data = Wide_data_no_control),
  Dur = lm(Dur_Ear_Sen ~ Genotype, data = Wide_data_no_control),
  TT10 = lm(Ear_TT10 ~ Genotype, data = Wide_data_no_control),
  TT30 = lm(Ear_TT30 ~ Genotype, data = Wide_data_no_control),
  TT50 = lm(Ear_TT50 ~ Genotype, data = Wide_data_no_control),
  TT70 = lm(Ear_TT70 ~ Genotype, data = Wide_data_no_control),
  TT90 = lm(Ear_TT90 ~ Genotype, data = Wide_data_no_control)
)

for(i in 1:length(linear_models_factorial)){
  print(names(linear_models_factorial)[i])
  par(mfrow=c(1,2))
  plot(linear_models_factorial[[i]], which = c(1,2))
  print(car::Anova(linear_models_factorial[[i]],type=3))
  display(linear_models_factorial[[i]])
}

print(paste(target,"AUC", "means"))
estimated_means <- emmeans(linear_models_factorial$AUC, ~ NAMA1 * NAMB1, nesting=NULL)
print(estimated_means)
print(confint(pairs(estimated_means), level = 0.95))

print(paste(target,"Dur", "means"))
estimated_means <- emmeans(linear_models_factorial$Dur, ~ NAMA1 * NAMB1, nesting=NULL)
print(estimated_means)
print(confint(pairs(estimated_means), level = 0.95))

label_yvals <- c(1300,18,68,70,71.5,73,75) #Height for Tukey Letters

for(i in 1:length(basic_plots)){
  par(mfrow=c(1,1))
  LABELS <- plot_anova_and_return_tukey_labels(linear_models[[i]], explanatory = "Genotype")
  formatted_plot <- basic_plots[[i]] +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.3) +
    genotype_scale_x +
    genotype_scale_fill +
    scale_y_continuous(expand = expansion(mult = 0.2, add = 0)) +
    labs(x = "", fill = genotype_lab) +
    geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = label_yvals[i], fill = "white", size = 6, label.size = 0) +
    #update y axis labels
    labs(y = str_wrap(variable_labels_ear[i], width = 30))+
    theme(legend.position = "none")
  print(formatted_plot)
  save_small_svg(formatted_plot, variable_list[i], "NAM1_Boxplot")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#FIG 3
#BOXPLOTS
basic_plots <- list(
  AUC_Leaf_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=AUC_Leaf_Sen, fill = Genotype)),
  Dur_Leaf_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=Dur_Leaf_Sen, fill = Genotype)),
  AUC_Ear_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=AUC_Ear_Sen, fill = Genotype)),
  Dur_Ear_Sen = ggplot(Wide_data_no_control, aes(x=Genotype, y=Dur_Ear_Sen, fill = Genotype))#,
  #Yield = ggplot(Wide_data_no_control, aes(x=Genotype, y=Yield, fill = Genotype)),
  #Protein = ggplot(filter(Wide_data_no_control, Rep == 1), aes(x=Genotype, y=Protein, fill = Genotype))
)

#Non-factorial version generates same pairwise comparisons by Tukey method
#As long as I keep the interaction term in the model
#I HAVE checked this, can get pairwise comps from factorial anova by using pairs(emmeans(linear_models[[i]])) 
#Use this as a work-around because my function plot_anova_and_return_tukey_labels() only takes 1 variable
linear_models <- list(
  AUC_Leaf_Sen = lm(AUC_Leaf_Sen ~ Rep + Genotype, data = Wide_data_no_control),
  Dur_Leaf_Sen = lm(Dur_Leaf_Sen ~ Rep + Genotype, data = Wide_data_no_control),
  AUC_Ear_Sen = lm(AUC_Ear_Sen ~ Rep + Genotype, data = Wide_data_no_control),
  Dur_Ear_Sen = lm(Dur_Ear_Sen ~ Rep + Genotype, data = Wide_data_no_control)#,
  #Yield = lm(Yield ~ Rep + Genotype, data = Wide_data_no_control),
  #Protein = lm(Protein ~ Genotype, data = filter(Wide_data_no_control, Rep == 1))
)


label_yvals <- c(1600,20,1300,16,9.8,13) #Height for Tukey Letters

formatted_plots <- list()

for(i in 1:length(basic_plots)){
  par(mfrow=c(1,1))
  LABELS <- plot_anova_and_return_tukey_labels(linear_models[[i]], explanatory = "Genotype")
  formatted_plots[[i]] <- basic_plots[[i]] +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.3) +
    genotype_scale_x +
    genotype_scale_fill +
    scale_y_continuous(expand = expansion(mult = 0.2, add = 0)) +
    labs(x = "", fill = genotype_lab) +
    geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = label_yvals[i], fill = "white", size = 8, label.size = 0) +
    #update y axis labels
    labs(y = str_wrap(variable_labels[i], width = 20))+
    theme(legend.position = "none")
}

#FACET
arranged_plots <- ggarrange(plotlist = formatted_plots, ncol = 2, nrow = 3)
arranged_plots

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#MODELS

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters
model_formula <- " ~ Genotype"
n_variables <- 1
explanatory <- "Genotype"

#Iterate through variables
# for(year in c(2021, 2022)){
  
  #Plot everything at once
  for(i in 1:length(variable_list_field)){
    
    print(variable_list_field[i])
    linear_model <- lm(formula(paste(variable_list_field[i], model_formula, sep = "")), data = filter(Wide_data_no_control, TRUE))
    
    par(mfrow=c(1,2))
    # plot(linear_model, which = c(1,2))
    #Crucial ANOVA and "compact letter display" labels
    labels <- plot_anova_and_return_labels_emm(linear_model, explanatory = explanatory, alpha = 0.05, variables = n_variables)
    #axes
    ymin <- min(filter(Wide_data_no_control, TRUE)[variable_list_unquote_field[i]], na.rm = TRUE)
    ymax <- max(filter(Wide_data_no_control, TRUE)[variable_list_unquote_field[i]], na.rm = TRUE)
    label_y <- ymax + 0.12*(ymax - ymin)
    axis_ymin <- ymin
    axis_ymax <- ymax + 0.18*(ymax - ymin)
    #final boxplot
    par(mfrow=c(1,1))
    
    plot <- ggplot(filter(Wide_data_no_control, TRUE),
                   aes_string(x = "Genotype", y = variable_list_field[i], fill = "Genotype")
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(height = 0)) +
      genotype_scale_x +
      genotype_scale_fill +
      labs(y = str_wrap(variable_labels_field[variable_list_unquote_field[i]], width = 23), x = NULL, fill = genotype_lab) +
      geom_label(
        data = labels,
        x = labels[, 1],
        label = labels[, 2],
        y = label_y,
        fill = "white",
        size = 16 / .pt, #size 12 for combined plots, size 16 - 20 for individual plots
        label.size = 0
      ) +
      coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
      my_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    print(plot)
    plot_list[[list_index]] <- plot
    list_index <- list_index + 1
    save_A5_svg(plot, str_c("Boxplot", "DDAH", variable_list_truncated_field[i], year, sep = "_"), gene = "NAM1")
  }
# }


#COMBINE PLOTS
length(variable_list_field) #10

#k = length(variable_list_field) - 1 #-1 because of flag leaf width only 2023

combined <- ggarrange(plot_list[[1]],
                      plot_list[[2]],
                      plot_list[[3]],
                      plot_list[[4]],
                      ncol = 2,
                      nrow = 2
)


save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated_field[1], sep = "_"), gene = "NAM1",
               n_panel_cols = 2, n_panel_rows = 2, ratio = 1)

combined <- ggarrange(plot_list[[5]],
                      plot_list[[6]],
                      plot_list[[7]],
                      plot_list[[8]],
                      ncol = 2,
                      nrow = 2
)


save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated_field[5], sep = "_"), gene = "NAM1",
               n_panel_cols = 2, n_panel_rows = 2, ratio = 1)

