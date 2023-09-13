#KASP plate results
#14/12/2021

#Borrowing the plate layout scripts from the qPCR stuff

#'25/01/21
#'KASP run 20/01/21
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5


#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(dplyr); library(ggplot2)
library(tidyr)
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
source("QPCR_functions.R")


#'EDIT section
#'EDIT this to contain the files for current experiment
#'Choose RAW DATA FILES and PLATE LAYOUTS
#'#NB layout must be the size of the full plate, i.e. 384 wells
#'NAM2 Gen7 KASP run 14/12/2021
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC transgenics/2021-11-10 NAC5 NAP Kronos greenhouse/KASP genotyping 2022-01-20")
n_plates <- 2
layout_samples_names <- c("2022-01-20_Kronos_layout_Plate1.csv", "2022-01-20_Kronos_layout_Plate1.csv")
KASP_raw_names <-  c("2022-01-20_Kronos_NAC5_NAP_35c.csv","2022-01-20_Kronos_NAC5_NAP_40c.csv")

use_layout <- TRUE
#primers = c("CE1773", "CE0895-2", "CE0237") #use this if the whole plate has the same primer
layout_primers_names <- rep("2022-01-20_Kronos_primers_Plate1.csv", 3) #use this if there are different primers


#'RUN
#'Combine all n plates into a single file
#'Add primer and sample information
combined <- c()

for (plate in 1:n_plates){
  layout_samples_plate <- read_csv(layout_samples_names[plate], col_names = FALSE)
  
  KASP_plate <- read_csv(KASP_raw_names[plate], skip = 18)
  table_samples <- make_table_samples(layout_samples_plate, n=384)
  KASP_v2_plate <- KASP_plate %>%
    tidyr::extract(MasterWell, c("Row_letter", "Col"), "([A-Z])([0-9]{2})") %>%
    mutate(Col = as.integer(Col)) %>%
    inner_join(table_samples, by = c("Row_letter", "Col")) %>% #Sample names and row and column numbers
    mutate(Plate = plate) #Technical replicates (assuming they run acrossways)
  if(use_layout){
    layout_primers_plate <- read_csv(layout_primers_names[plate], col_names = FALSE)
    table_primers <- rename(make_table_samples(layout_primers_plate, n=384), Primer = Sample)
    KASP_v3_plate <- KASP_v2_plate %>%
      inner_join(table_primers, by = c("Well", "Col","Row", "Row_letter"))
  }else{
    KASP_v3_plate <- KASP_v2_plate %>%
      mutate(Primer = primers[plate])
  }
  combined <- rbind(combined, KASP_v3_plate)
}


#EXPLORE
want = c("X:X","Y:Y")
het = c("X:Y")

calls_table <- combined %>%
  select(MasterPlate, Well, Primer, Sample, Call) %>%
  pivot_wider(id_cols = c(Well, Primer, Sample), names_from = MasterPlate, values_from = Call) %>%
  arrange(Primer, Sample)

calls_table_NAC5 <- calls_table %>%
  filter(Primer %in% c("PB479", "PB480", "PB481"))

calls_table_NAP <- calls_table %>%
  filter(Primer %in% c("PB462", "PB463"))

# calls_table <- combined %>%
#   filter(Call != "NTC") %>%
#   select(Well, Sample, Primer, Call) %>%
#   pivot_wider(id_cols = c(Well, Sample), names_from = Primer, values_from = Call) %>%
#   mutate(Hom = (CE1773 %in% want & `CE0895-2` %in% want & CE0237 %in% want),
#          Het = (CE1773 %in% het & `CE0895-2` %in% het & CE0237 %in% het)) %>% #Ah bother I can't work out what
#   mutate(across(where(is.character), ~as.factor(.)))
# 
# #summary statistics
# summary(calls_table)
# 
# calls_table_2 <- calls_table[grep("Cadenza", calls_table$Sample, invert = TRUE),]
# 
# summary(calls_table_2)
# calls_table_counts <- calls_table_2 %>%  
#   filter(CE1773 != "?", `CE0895-2` != "?", CE0237 != "?") %>%
#   group_by(CE1773, `CE0895-2`, CE0237) %>%
#   count()

#WRITE
write_csv(combined, file = paste(lubridate::today(), "_KASP_Kronos_NAC5_NAP_Plate1.csv"))
write_csv(calls_table_NAC5, file = paste(lubridate::today(), "_KASP_Kronos_NAC5_Plate1_calls.csv"))
write_csv(calls_table_NAP, file = paste(lubridate::today(), "_KASP_Kronos_NAP_Plate1_calls.csv"))
