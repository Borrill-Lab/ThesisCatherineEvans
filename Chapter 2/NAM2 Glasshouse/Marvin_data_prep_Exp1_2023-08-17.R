# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'Marvin data preparation
#'Combine rows for pooled Marvin stats
#'
#'13/08/2023
#'
#'17/08/2023
#'Use for NAM-2 2022 data
#'
#'Syntax rules: variables snake_case, Column Names in Caps, use pipe %>% where beneficial.
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'SOURCE
#'Packages and source files
require(lubridate)
require(readxl)
require(tidyverse) #Load tidyverse last to avoid masking

get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#OPEN 

setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-14 NAM2 greenhouse/2022-10-28_phenotyping_data", sep = ""))
Marvin_X50 <- read_excel("NAM2_X50 21_09_22.xls")

Marvin_X61 <- read_excel("NAM2_X61_Sown_21_09_22.xls")

date = "2023-08-02"
setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-14 NAM2 greenhouse/2023-01-12_clean_data", sep = ""))

NAM2_clean <- read_csv(paste("NAM2_phenotyping_clean", "_", date, ".csv", sep= ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY

Marvin <- full_join(Marvin_X50, Marvin_X61)

summary(Marvin)

length(unique(Marvin$ID)) #should equal num seed packets (currently 160 + 1 test - 1 missing X61)

Marvin_grouped <- Marvin %>%
  group_by(ID) %>%
  summarise(Seed_num = sum(`Main Seeds`), Weight = sum(`Weight(g)`))

summary(Marvin_grouped)

setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-14 NAM2 greenhouse/2023-01-12_clean_data", sep = ""))
write_csv(Marvin_grouped, file = paste('Marvin_clean_', today(), '.csv', sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#EXPLORE

#Figure out the best way to pool samples to get over 20g per sample
Marvin_grouped <- Marvin_grouped %>%
  filter(ID != "test") %>%
  mutate(Plant_num = as.numeric(ID))

Marvin_v2 <- NAM2_clean %>%
  select(Plant_num,
         Line_name,
         Cross,
         Genotype,
         Sample,
         Tiller_1,
         Rep) %>%
  left_join(Marvin_grouped, by = c("Plant_num")) %>%
  mutate(TKW = (Weight/Seed_num)*1000,
         Seed_per_tiller = (Seed_num/Tiller_1))

#Check out
qplot(Line_name, Seed_per_tiller, data = Marvin_v2)
qplot(TKW, Weight, data = Marvin_v2)
#Tiller number followed by TKW explain grain weight; seeds per tiller is fairly consistent
#X50_G3_46 high outlier - would be in range if it had 6 tillers rather than 5

Marvin_v2 %>%
  group_by(Line_name) %>%
  summarise(sum(Weight)/5)

#Rep numbers each set of lines from 1 to 10

Marvin_v3 <- Marvin_v2 %>%
  group_by(Cross, Genotype) %>%
  mutate(Pair = paste(Line_name, sample(rep(c("A", "B", "C", "D", "E"), 2), 10, replace = FALSE), sep = "_"))

Pair_summary <- Marvin_v3 %>%
  group_by(Line_name, Pair) %>%
  summarise(Weight = sum(Weight)) %>%
  mutate(NIR_valid = Weight > 20)

Pair_list <- Marvin_v3 %>%
  arrange(Pair) %>%
  mutate(Within = rep(c(1,2), 5)) %>%
  select(Cross, Genotype, Pair, Within, Sample, Plant_num, Weight) %>%
  pivot_wider(id_cols = c(Cross, Genotype, Pair), names_from = Within, values_from = c(Sample, Plant_num, Weight)) %>%
  mutate(Weight_total = Weight_1 + Weight_2)

setwd(paste(local_file_path, "NAM-2 TILLING/2022-09-14 NAM2 greenhouse/2023-01-12_clean_data", sep = ""))
write_csv(Pair_list, file = paste('Marvin_pairlist_', today(), '.csv', sep = ""))
