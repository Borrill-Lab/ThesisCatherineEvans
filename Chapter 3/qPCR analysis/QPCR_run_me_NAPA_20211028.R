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
setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
source("QPCR_functions.R")


#'EDIT section
#'EDIT this to contain the files for current qPCR experiment
#'NAP-A experiment run on 15/9/21 and 16/9/21
#'Choose RAW DATA FILES and PLATE LAYOUTS
#'#NB layout must be the size of the full plate, i.e. 384 wells
setwd("/Users/u1984449/OneDrive - Norwich BioScience Institutes/NAC transgenics/NAC expression R")
n_plates <- 3
layout_samples_names <- c("Layout_Samples_NAPA_Plate1_2021-09-16.csv","Layout_Samples_NAPA_Plate2_2021-09-16.csv","Layout_Samples_NAPA_Plate3_2021-09-16.csv")
qPCR_raw_names <-  c("2021-09-15 qPCR NAP-A individual Plate 1 Catherine 45c.xls", "2021-09-16 qPCR NAP individual Plate 2 Catherine 45c.xls","2021-09-16 qPCR NAP individual Plate 3 Catherine 45c.xls")
layout_primers_names <- rep("Layout_Primers_NAPA_AllPlates_2021-09-16.csv", 3) #same for all plates, so same file three times
primer_efficiencies <- read_csv("Primer_efficiencies_2021-09-09.csv", col_names = FALSE)

#Choose PRIMERS
#There should be one or more primers in each category. 
target_primers <- c("CENAP-7F-2R", "CENAP-5F-5R", "CENAP-8F-13R")
reference_primers <- "PB441_Act2"
normalize_sample <- "CTA10.con1-1"

#Choose WELL EXCLUSIONS
#This could even be a file with 3 columns, Plate, Well.Position and Reason.
exclude_m <- data.frame(Plate = c(1,2,2,1,2,2,2,2,2,2,2,3), Well.Position = c("H17","H8","H13","D21","A24","B19","B20","J17","J18","J21","P18","N13"), Reason = c(rep("Passive reference signal changes near CT",3), rep("Multiple Tm peaks",9))) #Any to be excluded manually

#Choose TITLES
#Used for file names
today_date <- lubridate::today()
experiment_title <- "NAPA_individual_"

#RUN section
#This section shouldn't need editing for a new qPCR experiment.

#Add columns for sample names, trios and reps
#Combine all n plates into a single file
qPCR_v2 <- c()
modifier <- 0 #Count trios across all plates, for ease of error checking
for (plate in 1:n_plates){
  qPCR_plate <- read_excel(qPCR_raw_names[plate], skip = 46)
  layout_samples_plate <- read_csv(layout_samples_names[plate], col_names = FALSE)
  table_samples <- make_table_samples(layout_samples_plate, n=384)
  layout_primers_plate <- read_csv(layout_primers_names[plate], col_names = FALSE)
  table_primers <- rename(make_table_samples(layout_primers_plate, n=384), Primer = Sample)
  qPCR_v2_plate <- qPCR_plate %>%
    inner_join(table_samples, by = "Well") %>% #Sample names and row and column numbers
    inner_join(table_primers, by = c("Well", "Col","Row")) %>% #Primer names
    mutate(Plate = plate, Trio = modifier + (Well+2) %/% 3, Rep = (Well+2) %% 3 +1) #Technical replicates (assuming they run acrossways)
  modifier <- unlist(max(qPCR_v2_plate$Trio))
  qPCR_v2 <- rbind(qPCR_v2, qPCR_v2_plate)
}

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

#Calculate Mean and SD.
qPCR_v3 <- qPCR_v2.5 %>%
  group_by(Primer, Trio, Sample) %>%
  filter(Include == TRUE) %>%
  summarise(Mean = mean(as.numeric(CT), na.rm=TRUE), SD = sd(as.numeric(CT), na.rm=TRUE), Count = sum(!is.na(CT)))

#Make it look more like the table in the book by comparing with the raw values. Where did this come from?
#NB means and SDs are from filtered values whereas these raw values are not.
qPCR_v4 <- qPCR_v2.5 %>%
  dcast(formula = Plate + Primer + Trio + Sample ~ Rep, value.var = "CT") %>%
  full_join(qPCR_v3)

#PFAFFL 2001 method
analysis <- data.frame()
for(target in target_primers){
  for(reference in reference_primers){
    print(c("Target",target,"Reference",reference))
    analysis <- rbind(analysis, analyse_pfaffl(qPCR_v3, target = target, reference = reference, primer_efficiencies = primer_efficiencies))
  }
}

#Save some documents
write.csv(exclude_report, file = paste("qPCR_", experiment_title,"exclusions_pfaffl_", today_date, ".csv", sep=""))
write.csv(qPCR_v4, file = paste("qPCR_", experiment_title,"replicate_CT_pfaffl_",today_date, ".csv", sep=""))
write.csv(analysis, file = paste("qPCR_", experiment_title,"analysis_pfaffl_",today_date, ".csv", sep=""))
