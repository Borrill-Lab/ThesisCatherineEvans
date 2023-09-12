#KASP plate results
#10/03/2023

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
#'NAM2 Exp2 KASP run 03/03/22 - 08/03/22

setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAM-2 TILLING/2023-02-15 NAM2 greenhouse/KASP NAM2 Exp2 2023-03-03")
n_plates <- 6
layout_samples_names <- c(rep("2023-03-07_layout_plates1-3.csv",3),rep("2023-03-07_layout_plates4-6.csv",3))
KASP_raw_names <-  c("2023-03-03_CE1656-COM-2_30c.csv","2023-03-03_CE1436-COM-3_30c.csv","2023-03-03_CE0231-COM-2_30c.csv",
                     "2023-03-07_CE1773-2_35c.csv","2023-03-07_CE0895-COM-2_35c.csv","2023-03-07_CE0237_30c.csv")
primers = c("CE1656-COM-2", "CE1436-COM-3", "CE0231-COM-2", "CE1773-2", "CE0895-COM-2", "CE0237") #use this if the whole plate has the same primer

combined <- c()

for (plate in 1:n_plates){
  layout_samples_plate <- read_csv(layout_samples_names[plate], col_names = FALSE)
  KASP_plate <- read_csv(KASP_raw_names[plate], skip = 18)
  table_samples <- make_table_samples(layout_samples_plate, n=384)
  KASP_v2_plate <- KASP_plate %>%
    tidyr::extract(MasterWell, c("Row_letter", "Col"), "([A-Z])([0-9]{2})") %>%
    mutate(Col = as.integer(Col)) %>%
    inner_join(table_samples, by = c("Row_letter", "Col")) %>% #Sample names and row and column numbers
    mutate(Plate = plate, Primer = primers[plate]) #Technical replicates (assuming they run acrossways)
  combined <- rbind(combined, KASP_v2_plate)
}

want = c("X:X","Y:Y")
het = c("X:Y")
  
calls_table <- combined %>%
  filter(Call != "NTC") %>%
  select(Well, Sample, Primer, Call) %>%
  pivot_wider(id_cols = c(Well, Sample), names_from = Primer, values_from = Call) %>%
  mutate(Hom_set1 = (`CE1656-COM-2` %in% want & `CE1436-COM-3` %in% want & `CE0231-COM-2` %in% want),
         Hom_set2 = (`CE1773-2` %in% want & `CE0895-COM-2` %in% want & `CE0237` %in% want)) %>%
         # Het = (`CE1656-COM-2` %in% het & `CE1436-COM-3` %in% het & `CE0231-COM-2` %in% het)) %>% #Ah bother I can't work out what
  mutate(across(where(is.character), ~as.factor(.)))

#summary statistics
summary(calls_table)

calls_table_2 <- calls_table[grep("Cadenza", calls_table$Sample, invert = TRUE),]
calls_table_2 <- calls_table_2 %>%
  filter(Sample != "0" & Sample != "POS") %>%
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

write.csv(calls_table_2, file = paste(lubridate::today(), "_KASP_NAM2_Exp2_calls.csv"))
