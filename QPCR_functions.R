#'19/04/21
#'Functions to be called by the main script
#'I want to create a new pipeline in R which contains the same analysis steps, equations and output tables as the Excel method.
#'Some steps in Excel method are not repeatable, prone to human error or otherwise problematic.
#'I want to try an R-based method which is more automated, repeatable, consistent, and once made can be re-used quickly for future experiments.

#test test test

#'Turn sample plate into vertical layout table
#'@param layout_samples layout data frame or matrix
#'@param n number of wells
make_table_samples <- function(layout_samples, n = 384) {
  layout_samples_2 <- t(as.matrix(layout_samples))
  table_samples <-
    data.frame(
      Well = 1:n,
      Col = row(layout_samples_2)[1:n],
      Row = col(layout_samples_2)[1:n],
      Row_letter = LETTERS[col(layout_samples_2)[1:n]],
      Sample = layout_samples_2[1:n],
      stringsAsFactors = FALSE
    )
  return(table_samples)
}


#'Returns a table of wells to exclude.
#'Uses find_sd_outliers to choose outliers with high SD.
#'@param qPCR_table data.frame containing CT values
#'@param exclude_m manually excluded wells
#'@param n number of wells
exclude_wells <-
  function(qPCR_table,
           exclude_m,
           max_SD = 1,
           max_diff = 1) {
    exclude_a <- qPCR_table %>%
      filter(Amp.Status != "Amp") %>%
      select(Plate, `Well.Position`)
    exclude_sd <- c()
    high_sd_trio <- c()
    for (t in 1:max(qPCR_table$Trio)) {
      trio_info <- qPCR_table %>%
        filter(Trio == t)
      CT_values <- trio_info %>%
        select(CT) %>%
        unlist()
      outlier <- find_sd_outliers(CT_values, max_SD, max_diff)
      if (outlier != FALSE) {
        high_sd_trio <- c(high_sd_trio, t)
        outlier_well <- trio_info %>%
          filter(Rep == outlier) %>%
          select(Plate, `Well.Position`)
        exclude_sd <- rbind(exclude_sd, outlier_well)
      }
    }
    print("High SD Trios")
    print(qPCR_table %>%
            filter(Trio %in% high_sd_trio))
    #Fudge factors to get rid of duplicate rows.
    exclude_a <- data.frame(exclude_a, Reason = "Amp Status")
    exclude_sd <-
      data.frame(exclude_sd, Reason = "Differs from technical replicates")
    return(rbind(exclude_m, exclude_a, exclude_sd))
  }

#' Return 1, 2 or 3 for replicate to exclude, or FALSE if none should be excluded.
#' This is literally what I do by eye
#' @param CT_values three CT values to compare
#' @param max_SD maximum SD to accept
#' @param max_diff pairwise difference between CT to accept/reject
find_sd_outliers <- function(CT_values,
                             max_SD = 1,
                             max_diff = 1) {
  #Return undetermined/NA as FALSE
  for (x in 1:3) {
    if (is.na(CT_values[x]) | CT_values[x] == "Undetermined") {
      return(FALSE)
    }
  }
  CT_values <- as.double(CT_values)
  #Return false if SD is low
  if (sd(CT_values, na.rm = TRUE) <= 1) {
    return(FALSE)
  }
  #Can use dist() function but this is more clear for 3 values
  diff12 <- abs(CT_values[1] - CT_values[2])
  diff13 <- abs(CT_values[1] - CT_values[3])
  diff23 <- abs(CT_values[2] - CT_values[3])
  #Return a number if it is more than 1 away from the others, and the others are less than 1 apart
  if (diff12 <= max_diff & diff13 > max_diff & diff23 > max_diff) {
    return(3)
  } else if (diff12 > max_diff &
             diff13 <= max_diff & diff23 > max_diff) {
    return(2)
  } else if (diff12 > max_diff &
             diff13 > max_diff & diff23 <= max_diff) {
    return(1)
  }
  return(FALSE)
}

#' Calculate DDCT and standard deviation columns using means table
#' Containing column Means, SD, Primer, and Sample
#' @param qPCR_means data frame
#' @param target target Primer
#' @param reference reference Primer
#' @return table containing DDCt and other columns
analyse_ddct <- function(qPCR_means, target, reference, normalize) {
  target_means <- qPCR_means %>%
    filter(Primer == target)
  reference_means <- qPCR_means %>%
    filter(Primer == reference)
  analysis <-
    left_join(reference_means, target_means, by = "Sample")
  print("Target and Reference data")
  print(head(target_means))
  print(head(reference_means))
  analysis <- analysis %>%
    mutate(DCt = Mean.y - Mean.x)
  normalize_value <-
    as.double(analysis[which(analysis$Sample == normalize), "DCt"])
  print(normalize_value)
  #Here's the important bit.
  analysis2 <- analysis %>%
    mutate(
      DDCt = as.double(DCt) - normalize_value,
      two_power_DDCt = 2 ^ -(DDCt),
      SD_compound = sqrt(SD.y ^ 2 + SD.x ^ 2),
      error_min = 2 ^ -(DDCt + SD_compound),
      error_max = 2 ^ -(DDCt - SD_compound)
    )
  print(head(analysis2))
  return(analysis2)
}

#' Calculate Fold change in expression using Pfaffl method
#' Table containing column Means, SD, Primer, and Sample
#' @param qPCR_means data frame
#' @param target target Primer
#' @param reference reference Primer
#' @param primer_efficiencies table of primer efficiencies
#' @return table containing fold_change and other columns
analyse_pfaffl <-
  function(qPCR_means,
           target,
           reference,
           primer_efficiencies) {
    target_means <- qPCR_means %>%
      filter(Primer == target)
    reference_means <- qPCR_means %>%
      filter(Primer == reference)
    analysis <-
      left_join(reference_means, target_means, by = "Sample")
    target_efficiency <-
      unlist(primer_efficiencies %>% filter(X1 == target) %>% select(X2)) #choose primer efficiency values from table. These values should be around 1
    reference_efficiency <-
      unlist(primer_efficiencies %>% filter(X1 == reference) %>% select(X2))
    
    print("Target and Reference data")
    print(head(target_means))
    print(head(reference_means))
    
    normalize_sample <-
      analysis %>% #simply use the smallest expression as the control
      mutate(DCt = Mean.y - Mean.x) %>%
      filter(!is.na(DCt) & Count.x>1 & Count.y>1) %>%
      arrange(desc(DCt)) %>%
      head(1)
    
    print(normalize_sample) #show which sample has been used
    #Here's the important bit.
    analysis2 <- analysis %>% mutate(
      target_amount = (target_efficiency + 1) ^ (Mean.y - normalize_sample$Mean.y),
      reference_amount = (reference_efficiency + 1) ^ (Mean.x - normalize_sample$Mean.x),
      fold_change = reference_amount / target_amount,
      log_fold_change = log(fold_change, base = exp(1))
    )
    print(head(analysis2))
    return(analysis2)
  }

#' Generate labels for TUKEY test
#' Ready to put on a graph
#' Function credits to R GRAPH GALLERY
#' @param TUKEY TUKEY test object
#' @param variable Explanatory variable e.g. genotype
#' @return dataframe of a b c labels ready to put on a graph, with column treatment
generate_label_df <- function(TUKEY, variable) {
  require(multcompView)
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][, 4]
  print(Tukey.levels)
  Tukey.labels <-
    data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$treatment) ,]
  return(Tukey.labels)
}