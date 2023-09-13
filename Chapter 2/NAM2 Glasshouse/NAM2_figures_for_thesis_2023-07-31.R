# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'NAM2 figures for THESIS
#'
#'2023-06-19
#'Make figures as quickly as possible to draft for lab meeting
#'
#'2023-06-30
#'Simplify genotype labels to dosage (single / double / triple mutant) for SEB
#'
#'2023-07-31
#'Pretty graphs for thesis
#'
#'# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#EDIT GENE AND DATE TO SUIT OCCASION
#date is the date the scripts were run through NAM2_data_preparation.R
gene = "NAM2"
date = "2023-06-30"

get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()

#Need to source every time to load packages
source("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_scripts/NAM2_thesis_graphs_source_2023-07-31.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN

date = "2023-07-31"
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_clean_data", sep = ""))

NAM2_wide_yday_23 <- read_csv(paste("NAM2_phenotyping_yday", "_", date, ".csv", sep= ""))
NAM2_SPAD_yday_23 <- read_csv(paste("NAM2_SPAD_yday", "_", date, ".csv", sep= ""))

date = "2023-08-02"
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-14 NAM2 greenhouse/2023-01-12_clean_data", sep = ""))

NAM2_wide_yday_22 <- read_csv(paste("NAM2_phenotyping_yday", "_", date, ".csv", sep= ""))
NAM2_SPAD_yday_22 <- read_csv(paste("NAM2_SPAD_yday", "_", date, ".csv", sep= ""))

NAM2_wide_yday_23 <- NAM2_wide_yday_23 %>%
  mutate(Year = 2023)
NAM2_SPAD_yday_23 <- NAM2_SPAD_yday_23 %>%
  mutate(Year = 2023)
NAM2_wide_yday_22 <- NAM2_wide_yday_22 %>%
  mutate(Year = 2022)
NAM2_SPAD_yday_22 <- NAM2_SPAD_yday_22 %>%
  mutate(Year = 2022)

NAM2_wide_yday <- full_join(NAM2_wide_yday_23, NAM2_wide_yday_22)
NAM2_SPAD_yday <- full_join(NAM2_SPAD_yday_23, NAM2_SPAD_yday_22)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OUTLIERS
print(paste("Gene:", gene, "Date:", date))

NAM2_outliers <- NAM2_wide_yday %>%
  filter(is.na(Plant_num)|grepl("snapped", Comment)|grepl("damage", Comment)|grepl("chlorosis", Comment)|is.na(Genotype))

NAM2_analysis <- NAM2_wide_yday %>%
  anti_join(NAM2_outliers, by = "Index")

NAM2_SPAD_analysis <- NAM2_SPAD_yday %>%
  anti_join(NAM2_outliers, by = "Index")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#SPAD plot

#This plot isn't compatible with facet so I'll plot individually
make_wilcox_plot <- function(data){
  Averages_plot_wilcox <- ggline(data = data, 
                                 x = "Week_approx", y = "SPAD", add = "mean_se",
                                 color = "Dosage", shape = "Dosage") +
    stat_compare_means(aes(group = Genotype), label = "p.signif",
                       label.y = 60, method = "wilcox.test", size = 16 / .pt, hide.ns = TRUE) +
    dosage_scale_col_2 +
    dosage_scale_shape_2 +
    spad_scale +
    scale_x_discrete(limits = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    labs(y = str_wrap(variable_labels["SPAD"], width = 25), x = "Weeks after Heading") +
    my_theme +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 65))
  return(Averages_plot_wilcox)
}


#SAVE SAVE SAVE
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/Figures", sep = ""))
plots <- list()

Averages_plot_wilcox <- make_wilcox_plot(filter(NAM2_SPAD_analysis,  Genotype %in% c("G1","G8") & Cross == "X50" & Year == 2022)) +
  ggtitle("A (2022 X50)")
save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_X50_2022_", gene = gene)
save_small_svg(Averages_plot_wilcox, "Averages_plot_wilcox_X50_2022_", gene = gene)
plots[[1]] <- Averages_plot_wilcox

Averages_plot_wilcox <- make_wilcox_plot(filter(NAM2_SPAD_analysis,  Genotype %in% c("G1","G8") & Cross == "X61" & Year == 2022)) +
  ggtitle("C (2022 X61)")
save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_X61_2022_", gene = gene)
save_small_svg(Averages_plot_wilcox, "Averages_plot_wilcox_X61_2022_", gene = gene)
plots[[2]] <- Averages_plot_wilcox

Averages_plot_wilcox <- make_wilcox_plot(filter(NAM2_SPAD_analysis,  Genotype %in% c("G1","G8") & Cross == "X50" & Year == 2023)) +
  ggtitle("B (2023 X50)")
save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_X50_2023_", gene = gene)
save_small_svg(Averages_plot_wilcox, "Averages_plot_wilcox_X50_2023_", gene = gene)
plots[[3]] <- Averages_plot_wilcox

Averages_plot_wilcox <- make_wilcox_plot(filter(NAM2_SPAD_analysis,  Genotype %in% c("G1","G8") & Cross == "X61" & Year == 2023)) +
  ggtitle("D (2023 X61)")
save_pdf_1(Averages_plot_wilcox, "Averages_plot_wilcox_X61_2023_", gene = gene)
save_small_svg(Averages_plot_wilcox, "Averages_plot_wilcox_X61_2023_", gene = gene)
plots[[4]] <- Averages_plot_wilcox

ggarrange(plots[[1]], plots[[3]], plots[[2]], plots[[4]])

#SAVE
save_panel_svg(ggarrange(plots[[1]], plots[[3]], plots[[2]], plots[[4]]), plot_name = "Averages_plot_wilcox", gene = "NAM2",
               n_panel_cols = 2, n_panel_rows = 2, ratio = 0.75)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#BOXPLOTS

#SAVE SAVE SAVE
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/Figures", sep = ""))

#VARIABLES to test
variable_list <-
  c(
    "Heading_date",
    "Leaf_senescence_DAH",
    "TT_SPAD30",
    "Dur_Leaf_Sen",
    "AUC_Leaf_Sen",
    "Peduncle_senescence_DAH",
    "Height",
    "Flag_leaf_length",
    "Tiller_1",
    "Flag_leaf_width"
  )
variable_list_unquote <- variable_list #NEEDS ADJUSTING WHEN NIR values added to escape weird characters
variable_list_truncated <- c(str_trunc(variable_list_unquote, width = 12, side = "right", ellipsis = "-"))


#MODELS

#Collect plots
plot_list <- list()
list_index = 1

#Set model parameters
model_formula <- " ~ Row + Block + Cross*Genotype"
n_variables <- 2
explanatory <- "Cross*Genotype"

#Iterate through variables
for(year in c(2022, 2023)){
  
  #Plot everything at once
  for(i in 1:length(variable_list)){
    if(year == 2022 & i == 10){ #no flag leaf length in 2022
      next()
    }
    print(variable_list[i])
    linear_model <- lm(formula(paste(variable_list[i], model_formula, sep = "")), data = filter(NAM2_analysis, Year == year))
    
    par(mfrow=c(1,2))
    # plot(linear_model, which = c(1,2))
    #Crucial ANOVA and "compact letter display" labels
    labels <- plot_anova_and_return_labels_emm(linear_model, explanatory = explanatory, alpha = 0.05, variables = n_variables)
    #axes
    ymin <- min(filter(NAM2_analysis, Year == year)[variable_list_unquote[i]], na.rm = TRUE)
    ymax <- max(filter(NAM2_analysis, Year == year)[variable_list_unquote[i]], na.rm = TRUE)
    label_y <- ymax + 0.12*(ymax - ymin)
    axis_ymin <- ymin
    axis_ymax <- ymax + 0.18*(ymax - ymin)
    #final boxplot
    par(mfrow=c(1,1))

    plot <- ggplot(filter(NAM2_analysis, Year == year),
                   aes_string(x = "Genotype", y = variable_list[i], fill = "Dosage")
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(position = position_jitter(height = 0)) +
      genotype_scale_x +
      dosage_scale_fill +
      labs(y = str_wrap(variable_labels[variable_list_unquote[i]], width = 23), x = NULL, fill = genotype_lab) +
      facet_wrap(vars(Cross), ncol = 2) +
      geom_label(
        data = labels,
        x = labels[, 2],
        label = labels[, 3],
        y = label_y,
        fill = "white",
        size = 12 / .pt, #size 12 for combined plots, size 16 - 20 for individual plots
        label.size = 0
      ) +
      coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
      my_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

    # print(plot)
    plot_list[[list_index]] <- plot
    list_index <- list_index + 1
    #save_A5_svg(plot, str_c("Boxplot", "faceted", "dosage", variable_list_truncated[i], year, sep = "_"), gene = gene)
  }
}


#COMBINE PLOTS
length(variable_list) #10

k = length(variable_list) - 1 #-1 because of flag leaf width only 2023

combined <- ggarrange(plot_list[[1]],
          plot_list[[1 + k]],
          plot_list[[3]],
          plot_list[[3 + k]],
          ncol = 1,
          nrow = 4
          )


save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated[1], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 4, ratio = 1)

combined <- ggarrange(plot_list[[5]],
                      plot_list[[5 + k]],
                      plot_list[[4]],
                      plot_list[[4 + k]],
                      ncol = 1,
                      nrow = 4
)


save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated[5], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 4, ratio = 1)

combined <- ggarrange(plot_list[[2 + k]],
                      plot_list[[6 + k]],
                      ncol = 1,
                      nrow = 2
)


save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated[2], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 2, ratio = 1)

combined <- ggarrange(plot_list[[9]],
                      plot_list[[9 + k]],
                      plot_list[[7]],
                      plot_list[[7 + k]],
                      ncol = 1,
                      nrow = 4
)


save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated[9], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 4, ratio = 1)

combined <- ggarrange(plot_list[[8]],
                      plot_list[[8 + k]],
                      plot_list[[10 + k]],
                      ncol = 1,
                      nrow = 3
)

save_panel_svg(combined, paste("Combined_boxplots", variable_list_truncated[8], sep = "_"), gene = gene,
               n_panel_cols = 2, n_panel_rows = 3, ratio = 1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#BOXPLOTS for G1 and G8 only
#Broken down by cross for the variables where cross seems to matter

#SAVE SAVE SAVE
setwd(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/Figures", sep = ""))

#Plot NO ANALYSIS
for(i in 1:length(variable_list)){
  par(mfrow=c(1,2))

  ymin <- min(NAM2_analysis[variable_list_unquote[i]], na.rm = TRUE)
  ymax <- max(NAM2_analysis[variable_list_unquote[i]], na.rm = TRUE)
  label_y <- ymax + 0.12*(ymax - ymin)
  axis_ymin <- ymin
  axis_ymax <- ymax + 0.18*(ymax - ymin)
  #final boxplot
  par(mfrow=c(1,1))
  
  plot <- ggplot(filter(NAM2_analysis, Genotype %in% c("G1", "G8")),
                 aes_string(x = "Genotype", y = variable_list[i], fill = "Dosage")
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(height = 0)) +
    genotype_scale_x_2 +
    dosage_scale_fill +
    labs(y = str_wrap(variable_labels[variable_list_unquote[i]], width = 23), x = NULL, fill = genotype_lab) +
    facet_wrap(vars(Year, Cross), ncol = 4) +
    coord_cartesian(ylim = c(axis_ymin, axis_ymax)) +
    my_theme + theme(legend.position = "none")
  
  print(plot)
  save_A5_svg(plot, str_c("Boxplot", "G1G8", variable_list_truncated[i], sep = "_"), gene = gene)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

