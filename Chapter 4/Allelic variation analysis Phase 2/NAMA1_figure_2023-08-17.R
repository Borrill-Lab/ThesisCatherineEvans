# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'NAMA1 figure
#'Make a stacked bar graph to represent NAM-A1 haplotype distributions
#'
#'17/08/2023
#'
#'
#'# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
require(readxl)
require(tidyverse)

get_local_file_path <- function(){
  local_file_path <- str_extract(getwd(), ".*OneDrive - Norwich BioScience Institutes/")
  if(is.na(local_file_path)){
    local_file_path <- str_extract(getwd(), ".*OneDrive-NorwichBioScienceInstitutes/")
  }
  print("The local file path is: "); print(local_file_path)
  return(local_file_path)
}
local_file_path <- get_local_file_path()

gene = "NAMA1"
#Need to source every time to load packages
source(paste(local_file_path, "NAM-2 TILLING/2023-02-15 NAM2 greenhouse/2023-06-15_scripts/NAM2_thesis_graphs_source_2023-07-31.R", sep = ""))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

NAMA1_haplotypes <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC genetic variation 2022/2022-10-06 RAGT KASP markers/2022-05-17 NAMA1 haplotypes.xlsx", 
                                                sheet = "R_summary")

NAMA1_haplotypes <- NAMA1_haplotypes %>%
  mutate(Haplotype = fct_rev(Haplotype),
         Source_category = fct_collapse(Source,
                               Published = c("Watkins",
                                             "INRA",
                                             "Elite_Europe",
                                             "Elite_Russia"),
                               Exome = c("RAGT_Exome"),
                               KASP = c("All",
                                        "France",
                                        "Germany", 
                                        "UK",
                                        "Czech",
                                        "Exotic",
                                        "Other")))

source_limits <- c("Watkins",
                   "INRA",
                   "Elite_Europe",
                   "Elite_Russia",
                   "RAGT_Exome",
                   "All",
                   "France",
                   "Germany", 
                   "UK",
                   "Czech",
                   "Exotic",
                   "Other")

stack_plot1 <- ggplot(NAMA1_haplotypes, aes(x=Source, y=Freq, fill=Haplotype))+
  geom_col(position = "stack") +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(type = "seq", palette = "YlGnBu") +
  scale_x_discrete(limits = source_limits) +
  theme(axis.text.x = element_text(angle = 90))

stack_plot1

source_limits <- c("Watkins",
                   "INRA",
                   "Elite_Europe",
                   "Elite_Russia",
                   "RAGT_Exome",
                   "France",
                   "Germany", 
                   "UK",
                   "Czech",
                   "Exotic",
                   "Other")

stack_plot2 <- ggplot(NAMA1_haplotypes, aes(x=Source, y=Count, fill=Haplotype))+
  geom_col(position = "stack") +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(type = "seq", palette = "YlGnBu") +
  scale_x_discrete(limits = source_limits) +
  theme(axis.text.x = element_text(angle = 90))
stack_plot2

save_A4_svg(ggarrange(stack_plot1, stack_plot2, ncol = 1), paste("haplotype_stackplot_", lubridate::today(), ".svg"), gene = "NAMA1")

