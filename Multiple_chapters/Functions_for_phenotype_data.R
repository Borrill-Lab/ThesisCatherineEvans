#Labels and formatting
# genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
# genotype_scale_colour <- scale_colour_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))



#Functions
#'returns a, b, c letters for plots
generate_label_df <- function(TUKEY, variable){
  stopifnot(require(multcompView))
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  # print(Tukey.levels)
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
}

#prints out anova and tukey HSD post-hoc test and returns a, b, c letters for plots
plot_anova_and_return_tukey_labels <- function(linear_model, explanatory = "Line_name"){
  print(anova(linear_model))
  aov <- aov(linear_model)
  TUKEY <- TukeyHSD(x=aov, conf.level=0.95)
  plot(TUKEY, las=1 , col="brown")
  LABELS <- generate_label_df(TUKEY , explanatory)
  return(LABELS)
}

#runs ggcorr with some formatting tweaks
ggcorr_with_theme <- function(data){
  stopifnot(require(GGally))
  ggcorr(
    data,
    method = c("pairwise", "pearson"),
    label = TRUE,
    label_round = 2,
    layout.exp = 2,
    hjust = 1
  )
}

#'prints out anova and emmeans test and returns a, b, c letters for plots
#'equivalent to tukey hsd if using anova with 1 explanatory factor
#'but this function works for more types of linear models
#'@param linear_model an lm object
#'@param explanatory the variable to compare
#'@param alpha p-value threshold
#'@param type type 2 or 3 Anova
#'@param ... additional arguments to pass to cld() e.g. adjust="tukey"
#'@return compact letter display for boxplots
plot_anova_and_return_labels_emm <- function(linear_model, explanatory, alpha =  0.05, type = 2, variables = 1, ...){
  stopifnot(require(emmeans))
  print(car::Anova(linear_model,type=type))
  estimated_means <- emmeans(linear_model, formula(paste("~", explanatory)), nesting=NULL)
  print(plot(estimated_means, comparisons = TRUE))
  labels <- generate_label_df_emm(estimated_means, explanatory, alpha, variables, ...)
  return(labels)
}

#'returns a, b, c letters for plots
#'@param estimated_means an emmGrid object
#'@param explanatory the variable to compare
#'@param alpha p-value threshold
#'@param variables number of explanatory factors
#'@param ... additional arguments to pass to cld() e.g. adjust="tukey"
#'@return compact letter display for boxplots
generate_label_df_emm <- function(estimated_means, explanatory, alpha = 0.05, variables = 1, ...){
  stopifnot(require(multcomp)); stopifnot(require(emmeans))
  labels <- cld(estimated_means, alpha = alpha, Letters = letters, ...)
  print(labels) #full emmeans table with compact letter display alongside
  if(variables == 1){
  labels <- labels[,c(1,7)]
  labels <- labels[order(labels[,1]) , ]
  } else if(variables == 2){
    labels <- labels[,c(1,2,8)]
    labels <- labels[order(labels[,1], labels[,2]) , ]
  } else if(variables == 3){
    labels <- labels[,c(1,2,3,9)]
    labels <- labels[order(labels[,1], labels[,2], labels[,3]) , ]
  } else if(variables == 4){
    labels <- labels[,c(1,2,3,4,10)]
    labels <- labels[order(labels[,1], labels[,2], labels[,3]) , ]
  } else {
    print("Error: Too many variables")
    stop()
    }
  return(labels)
}

#'Calculates the difference between 2 dates in degree days
#'EXCLUDES the temperature on the start date
#'So that if you put in the same date for start and end, the answer is 0
#'Half days are currently rounded UP
#'@param start_date date(s) in lubridate format
#'@param end_date date(s) in lubridate format
#'@param temperature_by_day a dataframe with columns "date" (inc all dates between start_date and end_date) and "mean_temp" (mean temperature each day)
degree_day_interval <- function(start_date, end_date, temperature_by_day){
  stopifnot(require(lubridate)); stopifnot(require(tidyverse))
  if(length(start_date)==length(end_date)){
    degree_days = c()
    for(i in 1:length(start_date)){
      focus_interval <- interval(start = start_date[i]+days(1), end = end_date[i], tzone = "GMT")
      degree_days[i] <- temperature_by_day %>%
        filter(date %within% focus_interval) %>%
        summarise(sum(mean_temp)) %>%
        unlist()
    }
    return(degree_days)
  }else{
    print("Error: Start date and end date need to be vectors of the same length")
  }
  
}