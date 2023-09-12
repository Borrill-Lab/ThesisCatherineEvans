# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Run Haplotype Pipeline for RAGT Exome Capture Data
#'
#'27/05/2022
#'practice data
#'
#'10/06/2022
#'Can we use real data??
#'Test with NAM-B2 from Pont et al
#'Remove 'NAME' column as absent
#'
#'11/07/2022
#'Run functions
#'
#'13/07/2022
#'Run on a loop for each gene in the RAGT dataset
#'
#'Built in R 4.2.0 and tidyverse 1.3.1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD PACKAGES
require(tidyverse)
require(readxl)
theme_set(theme_minimal())

#SOURCE haplotype_functions.R
source("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-05-27 Haplotypes/haplotype_functions_2022-07-13.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#OPEN FILES

setwd("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-07-01 RAGT Exome Capture")

data_name <- "NAC.qc.filter.nosingletons"
input <- read_delim("NAC.qc.filter.nosingletons.recode.vcf", delim = "\t", skip = 85)

head(input)

#For gene names
variant_summary_full <- read_csv("variant_summary_full_RAGT_Capture_2022-07-07.csv")

#For KASP marker quality
markers_analysis <- read_excel("markers_analysis_polymarker_NAC.nosingletons_2022-07-11.xlsx")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#DEFINE HAPLOTYPES

#Add columns to identify genes
#unique() to remove bonus rows for transcripts
gene_label <- variant_summary_full %>% select(CHROM_parts, POS_parts, REF:gene.name) %>% unique()

input <- input %>%
  inner_join(gene_label, by = c("#CHROM" = "CHROM_parts", "POS" = "POS_parts", "REF", "ALT")) %>%
  select(`#CHROM`, POS, REF, ALT, Uploaded_variation, Gene, gene.name, everything())

#Define genes
genes <- unlist(unique(gene_label$gene.name))

#Run name_haplotypes()
key_list <- list()
by_variety_list <- list()
for(i in 1:length(genes)){
  gene_input <- input %>%
    filter(gene.name == genes[i]) %>%
    select(-`#CHROM`, -POS, -REF, -ALT, -Gene, -gene.name, -ID, -QUAL, -FILTER, -INFO, -FORMAT)
  
  by_variety <- row_per_variety(gene_input, named_snps = TRUE)
  
  key <- name_haplotypes(by_variety)
  key <- mutate(key, gene.name = genes[i])
  key_list <- c(key_list, list(key))
  
  by_variety <- mutate(by_variety, gene.name = genes[i])
  by_variety_list <- c(by_variety_list, list(by_variety))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Categorise exclusions
#Based on quality of KASP marker

#Given a certain set of present / missing data, which haplotypes could each variety possibly have?
snp_table <- input %>%
  select(Uploaded_variation, gene.name) %>%
  unique() %>%
  left_join(markers_analysis, by = c("Uploaded_variation" = "Marker"))

snp_table_v1 <- snp_table %>%
  mutate(include = TRUE) %>%
  select(Uploaded_variation, gene.name, include)

snp_table_v2 <- snp_table %>%
  mutate(include = replace_na(ifelse(Cat == "Yes" | Cat == "Prob", TRUE, FALSE), FALSE)) %>%
  select(Uploaded_variation, gene.name, include)

snp_table_v3 <- snp_table %>%
  mutate(include = replace_na(ifelse(Comment == "Excellent", TRUE, FALSE), FALSE)) %>%
  select(Uploaded_variation, gene.name, include)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#IMPUTE HAPLOTYPES

#Given a certain set of present / missing data, which haplotypes could each variety possibly have?
output_list <- impute_haplotypes_multi(key_list = key_list, by_variety_list = by_variety_list, genes = genes,
                                       snp_table = snp_table_v1)
key_imputed_list <- output_list[[1]]
by_variety_summary <- output_list[[2]]

output_list <- impute_haplotypes_multi(key_list = key_list, by_variety_list = by_variety_list, genes = genes, 
                                       snp_table = snp_table_v2)
key_imputed_list_v2 <- output_list[[1]]
by_variety_summary_v2 <- output_list[[2]]

output_list <- impute_haplotypes_multi(key_list = key_list, by_variety_list = by_variety_list, genes = genes, 
                                       snp_table = snp_table_v3)
key_imputed_list_v3 <- output_list[[1]]
by_variety_summary_v3 <- output_list[[2]]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#GRAPH

#Simplify categories
by_variety_summary <- by_variety_summary %>%
  mutate(haplist_simplified = fct_collapse(haplist, hap_a = "hap_a", hap_b = "hap_b", hap_c = "hap_c", hap_d = "hap_d",
                                           hap_e = "hap_e", hap_f = "hap_f", other_level = "ambiguous"),
         set = "All")

by_variety_summary_v2 <- by_variety_summary_v2 %>%
  mutate(haplist_simplified = fct_collapse(haplist, hap_a = "hap_a", hap_b = "hap_b", hap_c = "hap_c", hap_d = "hap_d",
                                           hap_e = "hap_e", hap_f = "hap_f", other_level = "ambiguous"),
         set = "Yes_Prob")

by_variety_summary_v3 <- by_variety_summary_v3 %>%
  mutate(haplist_simplified = fct_collapse(haplist, hap_a = "hap_a", hap_b = "hap_b", hap_c = "hap_c", hap_d = "hap_d",
                                           hap_e = "hap_e", hap_f = "hap_f", other_level = "ambiguous"),
         set = "Excellent")

by_variety_combined <- full_join(full_join(by_variety_summary, by_variety_summary_v2), by_variety_summary_v3)


#Plot full summary
ggplot(by_variety_summary, aes(x = gene.name, fill = haplist)) +
  geom_histogram(stat = "count") +
  scale_fill_viridis_d()

ggplot(by_variety_summary, aes(x = gene.name, fill = haplist_simplified)) +
  geom_histogram(stat = "count") +
  scale_fill_viridis_d()

#Compare full vs Excluded
by_variety_combined %>%
  filter(gene.name %in% c("NAM-A1", "NAM-A2", "NAM-B2", "NAC-5A", "NAC-5B", "NAC-5D")) %>%
  ggplot(aes(x = set, fill = haplist)) +
  geom_histogram(stat = "count") +
  scale_fill_viridis_d() +
  facet_wrap(vars(gene.name)) +
  scale_x_discrete(limits = c("All", "Yes_Prob", "Excellent")) +
  labs(x = "Quality of KASP marker", y = "Count of varieties")

by_variety_combined %>%
  ggplot(aes(x = set, fill = haplist_simplified)) +
  geom_histogram(stat = "count") +
  scale_fill_viridis_d() +
  facet_wrap(vars(gene.name)) +
  scale_x_discrete(limits = c("All", "Yes_Prob", "Excellent")) +
  labs(x = "Quality of KASP marker", y = "Count of varieties")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#WRITE

#Write key_list
setwd("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-07-01 RAGT Exome Capture")

for(i in 1:length(genes)){
  write_csv(key_list[[i]], file = paste(paste("haplotype_key", genes[i], lubridate::today(), sep = "_"), ".csv", sep = ""))
}

for(i in 1:length(genes)){
  write_csv(key_imputed_list[[i]], file = paste(paste("haplotype_key_imputed", genes[i], lubridate::today(), sep = "_"), ".csv", sep = ""))
}
