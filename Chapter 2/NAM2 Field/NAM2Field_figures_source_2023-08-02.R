# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#Source script for
#'NAM2 Field thesis figures
#'
#'02/08/23
#'
#'# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#LOAD
require(tidyverse)
require(car) #Anova for non-orthogonal data
# require(emmeans) #marginal means estimates for when p-values are not enough
# require(arm)
require(GGally) #ggpairs()
require(ggpubr) #stat_cor()
theme_set(theme_bw())

#SOURCE
source('C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-03 NAC phenotyping/Functions_for_phenotype_data.R', echo=TRUE)

#OPEN
Long_data <- read_csv("Field_results_NAM2_long_2023-08-02.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

Wide_data <- read_csv("Field_results_NAM2_wide_DAH_2023-08-02.csv",
                      col_types = cols("Cross"= col_factor(),
                                       "Genotype"= col_factor(),
                                       "X"= col_factor()))

#'# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

