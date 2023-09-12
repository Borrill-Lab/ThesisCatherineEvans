# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'
#'19/10/2022
#'Try to make sense of KASP outputs for supervisory meeting
#'hmmmm
#'
#'18/08/2023
#'OH DEAR THE KASPS WERE SCRAMBLED this needs fixing, wasn't I suspicious all along?
#'
#'Built in R 4.2.0 and tidyverse 1.3.1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD
library(readr)
require(readxl)
require(ggpubr)
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

setwd(paste(local_file_path, "NAC genetic variation 2022/2022-10-06 RAGT KASP markers", sep = ""))

summary_SNP_Transcript_full_RAGT_Capture <- read_excel("../2022-07-01 RAGT Exome Capture/summary_SNP_Transcript_full_+_RAGT_Capture_edit_2023-08-18.xlsx")

KASP_results <- read_excel("Senescence associated_Nacgenes run on crossing parents to share with CE clean.xlsx", 
                                                                                            skip = 6)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# #18/08/23 UNSCRAMBLE
# 
# #The file I sent which had all the wrong names in "KASP name"
# scrambled_id_file <- read_excel("2022-07-22 NAC KASP markers Catherine Evans.xlsx")
# 
# #Remove all but names and primer sequence
# #Get name without A/G at end
# scrambled_id_subset <- scrambled_id_file %>%
#   select(`A + FAM`, COM, `KASP name`) %>%
#   extract(`KASP name`, "KASP_name_short", regex = "(.*)_[ACGT///]+$", remove = FALSE)
# 
# #The original Polymarker output which has the primers associated with the correct SNPs
# polymarker_output <- read_excel("../2022-07-01 RAGT Exome Capture/Polymarker_output_NAC.nosingletons_2022-07-06.xlsx")
# 
# #Join tables by the actual primer sequences
# key <- polymarker_output %>%
#   mutate(`A + FAM` = paste("GAAGGTGACCAAGTTCATGCT", A, sep = "")) %>%
#   left_join(scrambled_id_subset, by = c("A + FAM", "common" = "COM")) %>%
#   filter(!is.na(`KASP name`))
# 
# #annotate KASP results which are listed by "KASP name"
# KASP_results_primer_annot <- KASP_results %>%
#   mutate(KASP_name_short = ...1) %>%
#   left_join(key, by = c("KASP_name_short"))
# #ONLY USE Marker column from key
# KASP_results_annot <- KASP_results_primer_annot %>%
#   left_join(summary_SNP_Transcript_full_RAGT_Capture_2022_07_06,
#             by = c("Marker" = "Uploaded_variation"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#TIDY

#ONLY USE Marker column from key
KASP_results_annot <- KASP_results %>%
  left_join(summary_SNP_Transcript_full_RAGT_Capture,
                                by = c("Variant ID" = "Uploaded_variation", "gene.name"))

colnames(KASP_results_annot)

#OK so let's select the columns we might actually be interested in rn and drop misleading ones
  
KASP_results_annot_canonical <- KASP_results_annot %>%
  filter(CANONICAL == "YES") %>%
  dplyr::select(
    `Variant ID`,
    gene.name,
    #KASP output data
    All.Count.REF:Exotic.Freq.ALT,
    #SNP data
    `Chromosome/scaffold name`:SIFT, CANONICAL,
    #Public datasets frequency data
    He19_TOTAL:RAGT_ALT_FREQ
  )

KASP_frequencies <- KASP_results_annot_canonical %>%
  mutate(across(All.Freq.ALT:Exotic.Freq.ALT, .fns = ~.x/100)) %>%
  rename_with(.cols = He19_TOTAL:`He19_OTHER-SPECIES`, .fn = ~paste(.x, ".Freq.ALT", sep = "")) %>%
  rename(Pont19.Freq.ALT = Pont19_MAF,
         Watkins.Freq.ALT = Watkins_allele_freq,
         RAGT_Exome.Freq.ALT = RAGT_ALT_FREQ) %>%
  dplyr::select(`Variant ID`, REF.y, ALT.y, `Chromosome/scaffold name`:`Transcript stable ID`, gene.name, `Variant consequence`, `SIFT score`, ends_with("Freq.ALT"))

write_csv(KASP_frequencies, file = paste("KASP_frequencies_XP_v3_", lubridate::today(), ".csv", sep = ""))

KASP_frequencies_long <- KASP_frequencies %>%
  pivot_longer(cols = All.Freq.ALT:RAGT_Exome.Freq.ALT, names_to = c("Dataset", "spare1", "spare2"), values_to = "Frequency", names_sep = "\\.")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#PLOT

#Fresh start
KASP_frequencies <- read_csv("KASP_frequencies_XP_v3_2023-08-18.csv")

KASP_frequencies_long <- KASP_frequencies %>%
  pivot_longer(cols = All.Freq.ALT:RAGT_Exome.Freq.ALT, names_to = c("Dataset", "spare1", "spare2"), values_to = "Frequency", names_sep = "\\.")


dataset_limits <- c("All",        "France",     "Germany",    "UK",         "Czech",      "Exotic",
                    "RAGT_Exome", "He19_TOTAL", "He19_EuropeWest", "Pont19", "Watkins")

tile_plot <- ggplot(data=KASP_frequencies_long, aes(x=as.factor(`Variant ID`), y=Dataset, fill=Frequency)) +
  geom_tile()
tile_plot
tile_plot <- tile_plot +
  scale_fill_gradient(low = "white", high="navy") + #Colour scale
  scale_y_discrete(limits = dataset_limits) + #Scale of continents
  labs(title="Frequency of SNPs by dataset", x = "KASP `Variant ID` name", y="Dataset", fill="ALT Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) #Swivel the x-axis labels 90 degrees
tile_plot
tile_plot + facet_wrap(vars(gene.name), scale = "free_x", nrow = 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Try using ggarrange to improve facet plot?
# 
# plot_list <- c()
# 
# dataset_limits <- c("All",        "France",     "Germany",    "UK",         "Czech",      "Exotic",
#                     "RAGT_Exome")
# genes <- unique(KASP_frequencies_long_v2$gene.name)
# frequency_limits <- c(min(KASP_frequencies_long_v2$Frequency, na.rm = TRUE), max(KASP_frequencies_long_v2$Frequency, na.rm = TRUE))
# 
# for(i in 1:length(genes)){
#   tile_plot <-
#     KASP_frequencies_long_v2 %>%
#     filter(gene.name == genes[i]) %>%
#     ggplot(aes(x=as.factor(`Variant ID`), y=Dataset, fill=Frequency)) +
#     geom_tile() +
#     scale_fill_gradient(low = "white", high="navy", limits = frequency_limits) + #Colour scale
#     scale_y_discrete(limits = dataset_limits, expand = c(0,0)) + #Scale of continents
#     scale_x_discrete(expand = c(0,0)) +
#     labs(title=genes[i], x = "KASP `Variant ID` name", y="Dataset", fill="ALT Frequency") +
#     ggtitle(genes[i]) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle=90), legend.position = "none",
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.title.y = element_blank()
#           ) #Swivel the x-axis labels 90 degrees
#   #print(tile_plot)
#   plot_list[[i]] <- tile_plot
# }
# 
# ggarrange(plotlist = plot_list, ncol = 9, nrow= 1, widths = c(2,1,3,1,1,1,1,5,3))
# #NOPE THE FORMATTING OF THIS IS AWFUL
# #The egg or ggh4x package might solve this but I can't get packages to install rn

frequency_limits <- c(0,1)

tile_plot <-
  KASP_frequencies_long %>%
  ggplot(aes(x=fct_reorder(as.factor(`Variant ID`), gene.name, .desc = TRUE), y=Dataset, fill=Frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high="navy", limits = frequency_limits) + #Colour scale
  scale_y_discrete(limits = dataset_limits, expand = c(0,0)) + #Scale of continents
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "KASP marker name", y="Dataset", fill="ALT Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) + #Swivel the x-axis labels 90 degrees
  coord_flip() #nah it's better the other way around
tile_plot

ggsave(
  filename = paste("XP_KASP_tile_plot", "_", lubridate::today(), ".svg", sep = ""),
  tile_plot,
  width = 160,
  height = 120,
  units = "mm",
  device = "svg"
)

dataset_limits <- c("All",        "France",     "Germany",    "UK",         "Czech",      "Exotic",
                    "RAGT_Exome")

tile_plot <-
  KASP_frequencies_long %>%
  ggplot(aes(x=fct_reorder(as.factor(`Variant ID`), gene.name, .desc = TRUE), y=Dataset, fill=Frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high="navy", limits = frequency_limits) + #Colour scale
  scale_y_discrete(limits = dataset_limits, expand = c(0,0)) + #Scale of continents
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "KASP marker name", y="Dataset", fill="ALT Frequency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) + #Swivel the x-axis labels 90 degrees
  coord_flip() #nah it's better the other way around
tile_plot

ggsave(
  filename = paste("XP_KASP_tile_plot", "_RAGTonly_", lubridate::today(), ".svg", sep = ""),
  tile_plot,
  width = 160,
  height = 120,
  units = "mm",
  device = "svg"
)
