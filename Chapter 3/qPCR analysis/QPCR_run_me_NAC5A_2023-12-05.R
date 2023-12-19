#'19/04/21
#'I want to create a new pipeline in R which contains the same analysis steps, equations and output tables as the Excel method.
#'Some steps in Excel method are not repeatable, prone to human error or otherwise problematic.
#'I want to try an R-based method which is more automated, repeatable, consistent, and once made can be re-used quickly for future experiments.
#'
#'13/05/21
#'Edit to include multiple plates at once
#'
#'31/05/21
#'Iteration through all primers. Built to make files with 20210531.csv
#'
#'21/09/21
#'Run the script for NAP-A individual samples on 21/09/21. Updated to take primers from a table automatically as well as samples.
#'
#'28/10/21
#'Implement Pfaffl 2001 method in order to include primer efficiencies. The downside is that I don't know how to calculate the
#'error bars.
#'
#'17/12/21
#'NBI setwd()
#'
#'05/12/23
#'Corrections: normalize more appropriately
#'Repeat for NAC5A
#'
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5


#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(dplyr); library(ggplot2)
library(reshape2)
library(lubridate)
require(stringr)
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
source("QPCR_functions.R")


#'EDIT section
#'EDIT this to contain the files for current qPCR experiment
#'NAP-A experiment run on 15/9/21 and 16/9/21
#'Choose RAW DATA FILES and PLATE LAYOUTS
#'#NB layout must be the size of the full plate, i.e. 384 wells
# setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
n_plates <- 3
layout_samples_names <- c("Layout_Samples_Plate1_250321.csv", "Layout_Samples_Plate2_300321.csv", "Layout_Samples_Plate3_edit_300321.csv")
qPCR_raw_names <-  c("2021-03-25 qPCR NAC5 individual Plate 1 Catherine 45c.xls","2020-03-30 qPCR NAC5 individual Plate 2 Catherine 45c.xls", "2021-03-30 qPCR NAC5 individual Plate 3 Catherine 45c.xls")
primer_efficiencies <- read_csv("Primer_efficiencies_NAP_2021-07-27.csv", col_names = FALSE)

#Choose PRIMERS
#There should be one or more primers in each category. 
target_primers <- c("CENAC5-7", "CENAC5-5", "CENAC5-3")
reference_primers <- "PB441_Act2"
#NO
#normalize_sample <- "CTA10.con1-1"

#Choose WELL EXCLUSIONS
#This could even be a file with 3 columns, Plate, Well.Position and Reason.
exclude_m <- data.frame(Plate = c(1,1,1), Well.Position = c("O12","P21","H22"), Reason = "Melt curve") #Any to be excluded manually

#Choose TITLES
#Used for file names
today_date <- lubridate::today()
experiment_title <- "NAC5A_individual_"

#RUN section
#This section shouldn't need editing for a new qPCR experiment.

#Add columns for sample names, trios and reps
#Combine all n plates into a single file
qPCR_v2 <- c()
modifier <- 0 #Count trios across all plates, for ease of error checking
for (plate in 1:n_plates){
  layout_samples_plate <- read_csv(layout_samples_names[plate], col_names = FALSE)
  qPCR_plate <- read_excel(qPCR_raw_names[plate], skip = 46)
  table_samples <- make_table_samples(layout_samples_plate, n=384)
  qPCR_v2_plate <- qPCR_plate %>%
    inner_join(table_samples, by = "Well") %>% #Sample names and row and column numbers
    mutate(Plate = plate, Trio = modifier + (Well+2) %/% 3, Rep = (Well+2) %% 3 +1) #Technical replicates (assuming they run acrossways)
  modifier <- unlist(max(qPCR_v2_plate$Trio))
  qPCR_v2 <- rbind(qPCR_v2, qPCR_v2_plate)
}

#The Primer names were added in the raw files in the qPCR software
colnames(qPCR_v2)[3] <- "Primer"

colnames(qPCR_v2) <- gsub(" ", ".", colnames(qPCR_v2)) #Replace spaces with dots
            
#Denote unwanted wells for excluding.
exclude <- exclude_wells(qPCR_v2, exclude_m = exclude_m)

exclude_report <- qPCR_v2 %>%
  inner_join(exclude, by=c("Plate", "Well.Position")) %>%
  select(Plate, Well.Position, Primer, Sample,CT, Amp.Status, Reason) %>%
  arrange(Plate, Well.Position)

View(exclude_report) #Check it's OK before moving on

duplicate_outliers <- exclude %>%
  select(Plate, Well.Position) %>%
  duplicated() %>%
  which()

if(length(duplicate_outliers) != 0){
  exclude <- exclude[-duplicate_outliers,]
}
exclude <-  mutate(exclude, Include = FALSE)

#Unique Plate and Well Position values.
qPCR_v2.5 <- left_join(qPCR_v2, exclude, by = c("Plate", "Well.Position"))
qPCR_v2.5[which(is.na(qPCR_v2.5$Include)), "Include"] <- TRUE

#TIDY
#Tidy the file to add the necessary extra columns for the correction normalization
#Error "Expected 2 pieces. Missing pieces filled with `NA` in 48 rows" expected due to NTCs
genotype_key <- read_excel("C:/Users/evansc/OneDrive - Norwich Bioscience Institutes/NAC transgenics/NAC expression R/NAC-5A_genotype_key.xlsx")
qPCR_v2.6 <- qPCR_v2.5 %>%
  mutate(Sample_2 = str_remove(Sample, "CTA"), Sample_3 = str_remove(Sample, "CTA")) %>% #Don't want to lose sample column
  tidyr::separate(Sample_2, into = c("Line_name", "PlantNum"), sep = "-") %>%
  left_join(genotype_key, by = "Line_name")

#Calculate Mean and SD.
qPCR_v3 <- qPCR_v2.6 %>%
  group_by(Primer, Parental_genotype, Line_name, PlantNum, Trio, Sample) %>%
  filter(Include == TRUE) %>%
  summarise(Mean = mean(as.numeric(CT), na.rm=TRUE), SD = sd(as.numeric(CT), na.rm=TRUE), Count = sum(!is.na(CT)))

#Make it look more like the table in the book by comparing with the raw values. Where did this come from?
#NB means and SDs are from filtered values whereas these raw values are not.
qPCR_v4 <- qPCR_v2.6 %>%
  dcast(formula = Plate + Primer + Trio + Sample ~ Rep, value.var = "CT") %>%
  full_join(qPCR_v3)


#PFAFFL 2001 method
analysis <- data.frame()

table1 <- qPCR_v3 %>%
  group_by(Primer, Parental_genotype) %>%
  summarise(max(Mean), min(Mean), mean(Mean), median(Mean), n())
#Use mean value because if we've already averaged between replicates it's probably okay (and means and medians are v similar)
#Use single copy because controls don't amplify with some primers

for(target in target_primers){
  for(reference in reference_primers){
    #Normalize against mean of single-copy samples
    normalize_target <- table1 %>%
        filter(Primer == target & Parental_genotype == "single") %>%
        ungroup() %>%
        select(`mean(Mean)`) %>% unlist()
    normalize_reference <- table1 %>%
        filter(Primer == reference & Parental_genotype == "single") %>%
        ungroup() %>%
        select(`mean(Mean)`) %>% unlist()
    print(c("Target",target,"Reference",reference))
    print(c("Target normalization mean",target,"Reference normalization mean",reference))
    analysis <- rbind(analysis, analyse_pfaffl_v2(qPCR_v3,
                                                  target = target,
                                                  reference = reference,
                                                  primer_efficiencies = primer_efficiencies,
                                                  normalize_target = normalize_target,
                                                  normalize_reference = normalize_reference
                                                  )
                      )
  }
}

#Save some documents
write.csv(exclude_report, file = paste("qPCR_", experiment_title,"exclusions_pfaffl_", today_date, ".csv", sep=""))
write.csv(qPCR_v4, file = paste("qPCR_", experiment_title,"replicate_CT_pfaffl_",today_date, ".csv", sep=""))
write.csv(analysis, file = paste("qPCR_", experiment_title,"analysis_pfaffl_",today_date, ".csv", sep=""))
