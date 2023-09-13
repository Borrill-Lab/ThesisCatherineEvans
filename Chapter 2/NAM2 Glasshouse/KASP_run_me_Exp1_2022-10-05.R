#KASP plate results
#05/10/22

#Borrowing the plate layout scripts from the qPCR stuff

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
#'NAM2 Greenhouse KASP run 05/10/2022
#'NEW primers!!
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2022-09-14 NAM2 greenhouse/2022-10-05 KASP NAM2 greenhouse")
n_plates <- 3
layout_samples_names <- c("2022-10-05_NAM2_greenhouse_plate1.csv","2022-10-05_NAM2_greenhouse_plate2.csv", "2022-10-05_NAM2_greenhouse_plate3.csv")
KASP_raw_names <-  c("2022-10-05_NAM2_greenhouse_plate1_30c.csv","2022-10-05_NAM2_greenhouse_plate2_35c.csv","2022-10-05_NAM2_greenhouse_plate3_30c_editedcalls.csv")

use_layout <- TRUE
# primers = c("CE1656-COM-2", "CE1436-COM-3", "CE0231-COM-2", "CE1773-2", "CE0895-COM-2", "CE0237") #use this if the whole plate has the same primer
layout_primers_names <- c("2022-10-05_NAM2_greenhouse_primers1.csv", "2022-10-05_NAM2_greenhouse_primers2.csv", "2022-10-05_NAM2_greenhouse_primers3.csv")

skip = 16 + 2

combined <- c()

for (plate in 1:n_plates){
  layout_samples_plate <- read_csv(layout_samples_names[plate], col_names = FALSE)
  
  KASP_plate <- read_csv(KASP_raw_names[plate], skip = skip) 
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

want = c("X:X","Y:Y")
het = c("X:Y")
  
calls_table <- combined %>%
  filter(Call != "NTC") %>%
  select(Well, Sample, Primer, Call) %>%
  pivot_wider(id_cols = c(Sample), names_from = Primer, values_from = Call) %>%
  mutate(Hom_set1 = (`CE1656-COM-2` %in% want & `CE1436-COM-3` %in% want & `CE0231-COM-2` %in% want),
         Hom_set2 = (`CE1773-2` %in% want & `CE0895-COM-2` %in% want & `CE0237` %in% want)) %>%
         # Het = (`CE1656-COM-2` %in% het & `CE1436-COM-3` %in% het & `CE0231-COM-2` %in% het)) %>% #Ah bother I can't work out what
  mutate(across(where(is.character), ~as.factor(.)))

#summary statistics
summary(calls_table)

# calls_table_2 <- calls_table[grep("Cadenza", calls_table$Sample, invert = TRUE),]
calls_table_2 <- calls_table %>%
  filter(Sample != "0") %>%
  arrange(Sample)

summary(calls_table_2)

#COUNT
calls_table_counts <- calls_table_2 %>%  
  filter(`CE1656-COM-2` != "?", `CE1436-COM-3` != "?", `CE0231-COM-2` != "?") %>%
  group_by(`CE1656-COM-2`, `CE1436-COM-3`, `CE0231-COM-2`) %>%
  count()

calls_table_counts

calls_table_counts <- calls_table_2 %>%  
  filter(`CE1773-2` != "?", `CE0895-COM-2` != "?", CE0237 != "?") %>%
  group_by(`CE1773-2`, `CE0895-COM-2`, CE0237) %>%
  count()

calls_table_counts

write.csv(calls_table_2, file = paste(lubridate::today(), "_KASP_NAM2_Greenhouse_calls.csv"))
