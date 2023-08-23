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
#'This is the main script.
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.5


#'SOURCE section
#'Packages and source file QPCR_functions.R
library(readr); library(readxl)
library(dplyr); library(ggplot2)
library(reshape2)
setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics/NAC expression R")
source("QPCR_functions.R")


#'EDIT section
#'EDIT this to contain the files for current qPCR experiment
#'NAC-5A experiment run 25/03/21 and 30/03/21
#'Choose RAW DATA FILES and PLATE LAYOUTS
#'#NB layout must be the size of the full plate, i.e. 384 wells
setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics/NAC expression R")
n_plates <- 3
layout_samples_names <- c("Layout_Samples_Plate1_250321.csv", "Layout_Samples_Plate2_300321.csv", "Layout_Samples_Plate3_edit_300321.csv")
qPCR_raw_names <-  c("2021-03-25 qPCR NAC5 individual Plate 1 Catherine 45c.xls","2020-03-30 qPCR NAC5 individual Plate 2 Catherine 45c.xls", "2021-03-30 qPCR NAC5 individual Plate 3 Catherine 45c.xls")

#Choose PRIMERS
#There should be one or more primers in each category. 
target_primers <- c("CENAC5-7", "CENAC5-5", "CENAC5-3")
reference_primers <- "PB441_Act2"
normalize_sample <- "CTA5.4-1"

#Choose WELL EXCLUSIONS
#This could even be a file with 3 columns, Plate, Well.Position and Reason.
exclude_m <- data.frame(Plate = c(1,1,1), Well.Position = c("O12","P21","H22"), Reason = "Melt curve") #Any to be excluded manually


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
colnames(qPCR_v2)[3] <- "Primer"
colnames(qPCR_v2) <- gsub(" ", ".", colnames(qPCR_v2)) #Replace spaces with dots

#Denote unwanted wells for excluding.
exclude <- exclude_wells(qPCR_v2, exclude_m = exclude_m)

exclude_report <- qPCR_v2 %>%
  inner_join(exclude, by=c("Plate", "Well.Position")) %>%
  select(Plate, Well.Position, Primer, Sample,CT, Amp.Status, Reason) %>%
  arrange(Plate, Well.Position)
View(exclude_report)

View(exclude_report) #Check it's OK before moving on

duplicate_outliers <- exclude %>%
  select(Plate, Well.Position) %>%
  duplicated() %>%
  which()
exclude <- exclude[-duplicate_outliers,] %>%
  mutate(Include = FALSE)


#Unique Plate and Well Position values.
qPCR_v2 <- left_join(qPCR_v2, exclude[-duplicate_outliers,], by = c("Plate", "Well.Position"))
qPCR_v2[which(is.na(qPCR_v2$Include)), "Include"] <- TRUE

#Calculate Mean and SD.
qPCR_v3 <- qPCR_v2 %>%
  group_by(Primer, Trio, Sample) %>%
  filter(Include == TRUE) %>%
  summarise(Mean = mean(as.numeric(CT), na.rm=TRUE), SD = sd(as.numeric(CT), na.rm=TRUE), Count = sum(!is.na(CT)))

#Make it look more like the table in the book by comparing with the raw values. Where did this come from?
#NB means and SDs are from filtered values whereas these raw values are not.
qPCR_v4 <- qPCR_v2 %>%
  dcast(formula = Plate + Primer + Trio + Sample ~ Rep, value.var = "CT") %>%
  full_join(qPCR_v3)

#Got up to here 13/05/21
#31/05/21
analysis <- data.frame()
for(target in target_primers){
  for(reference in reference_primers){
    print(c("Target",target,"Reference",reference))
    analysis <- rbind(analysis, analyse_ddct(qPCR_v3, target = target, reference = reference, normalize = normalize_sample))
  }
}

#Save some documents
today_date <- "CTA5.4_20210531"
write.csv(exclude_report, file = paste("qPCR_NAC5A_","exclusions_", today_date, ".csv", sep=""))
write.csv(qPCR_v4, file = paste("qPCR_NAC5A_","replicate_CT_",today_date, ".csv", sep=""))
write.csv(analysis, file = paste("qPCR_NAC5A_","analysis_",today_date, ".csv", sep=""))
