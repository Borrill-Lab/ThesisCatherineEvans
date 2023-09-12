# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Summary Tables
#'
#09/02/21
#Understand specific subset of variants with missense mutations
#Make a combined table with as much information as possible
#
#09/05/22
#Include all variants, not just missense
#Add cultivar and continent info
#Full output table + brief output table with fewer columns, no counts
#Remove plants table (see 09/02/21 for this one)
#'
#'10/06/22
#'Run with data from RAGT Exome Capture
#'This doesn't need allele_freq as it's already been added
#'Also join with existing summary tables copying summary_table_prep_2022-06-17
#'
#'1.	Sort equivalent column names, choosing more informative (if longer) names
#'2.	Each column unique to SNP * Variant source must either be removed or duplicated for each variant source.
#'3.	Join into one big table
#'4.	Remove X columns
#'5.	Pivot to 1 row for each SNP * Transcript
#'6.	Select “full” table with all essential and optional columns
#'7.	Select “intermediate” table with all essential and some optional columns
#'8.	Select “brief” table with essential columns
#
#Built in R 4.0.5
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#LOAD PACKAGES
theme_set(theme_bw())
library(tidyverse)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#OPEN FILES

setwd("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-07-01 RAGT Exome Capture")

#Info on allele frequency etc compiled by PolyMarker Input Prep
SummaryTable <- read_csv("Summary_NAC.qc.filter.nosingletons_2022-07-06.csv")

#Gene names
Gene_summary <- read_delim("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/NAC_target_genes_Catherine.txt")

#Summary table for all publically available data
combined_summary <- read_csv("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-05-27 Haplotypes/summary_SNP_Transcript_Variantsource_2022-07-06.csv")

#Open the VeP file
variant_effect <- read_delim(file = "VeP_output_NAC.qc.filter.nosingletons_2022-07-06.txt", delim = "\t", skip = 0)
head(variant_effect)


# RUN ALL THE VARIANTS

#Combine with VeP table

variant_summary_full <- SummaryTable %>%
  left_join(variant_effect, by = c("Uploaded_variation" = "#Uploaded_variation")) %>%
  select(!starts_with("X")) %>%
  left_join(Gene_summary[c(1,5)], by = c("Gene" = "Gene stable ID")) %>%
  rename(gene.name = `Gene name`) %>%
  select(CHROM, POS, REF, ALT, Uploaded_variation, Gene, gene.name, CDS_position, Consequence, SIFT, everything())


#Filter
variant_effect_brief <- variant_effect %>%
  select(`#Uploaded_variation`, Gene, CDS_position, Consequence, SIFT)

variant_summary_brief <- SummaryTable %>%
  left_join(variant_effect_brief, by = c("Uploaded_variation" = "#Uploaded_variation")) %>%
  left_join(Gene_summary[c(1,5)], by = c("Gene" = "Gene stable ID")) %>%
  rename(gene.name = `Gene name`) %>%
  select(CHROM, POS, REF, ALT, Uploaded_variation, Gene, gene.name, CDS_position, Consequence, SIFT, everything()) %>%
  group_by(CHROM, POS, REF, ALT) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.x)), collapse = "")))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TIDY FILES

#'1.	Sort equivalent column names, choosing more informative (if longer) names
#'2.	Each column unique to SNP * Variant source must either be removed or duplicated for each variant source.

colnames(combined_summary); colnames(variant_summary_full)

variant_summary_full_RAGT_Capture <- variant_summary_full %>%
  mutate(`Variant source` = "RAGT_Capture", `Variant name` = Uploaded_variation) %>%
  mutate(across(REF_FREQ:N_CHR, .fns = NULL, .names = "RAGT_{.col}"), .keep = "unused") %>%
  mutate(across(c(cDNA_position, CDS_position, Protein_position), as.numeric))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#JOIN
#'3.	Join into one big table
#'

combined_summary_v2 <- full_join(combined_summary, variant_summary_full_RAGT_Capture,
          by = c("Variant name",
            "Variant source",
                 "Chromosome/scaffold name" = "CHROM2",
                 "Chromosome/scaffold position start (bp)" = "POS",
                 "Chromosome/scaffold position end (bp)" = "ENDPOS",
                 "REF",
                 "ALT",
                 "Gene stable ID" = "Gene",
                 "Variant consequence" = "Consequence",
                 "Variant start in cDNA (bp)" = "cDNA_position",
                 "Variant start in translation (aa)" = "Protein_position", 
                 "Variant start in CDS (bp)" = "CDS_position",
                 "Transcript stable ID" = "Feature",
                 "Uploaded_variation", "gene.name", "SIFT", "IMPACT", "Feature_type", "BIOTYPE", "EXON", "INTRON", "Amino_acids", "Codons",
                 "Existing_variation", "DISTANCE", "CANONICAL")
)

#DEFINE essential and optional columns
column_names <- colnames(combined_summary)

essential_1 <- c("Variant name", "Variant source",
                 "Chromosome/scaffold name",
                 "Chromosome/scaffold position start (bp)",
                 "Chromosome/scaffold position end (bp)",
                 "REF",
                 "ALT",
                 "Strand",
                 "Uploaded_variation",
                 "Gene stable ID",
                 "gene.name",
                 "Transcript stable ID"
)
essential_2 <- c(column_names[c(31:47)], "RAGT_REF_FREQ", "RAGT_ALT_FREQ", "RAGT_N_CHR")

optional_1 <- c(column_names[c(13:16)])

optional_2 <- c(column_names[c(17:30)])

optional_3 <- c(column_names[c(48:57)])


#4.	Remove X columns
combined_summary_v3 <- combined_summary_v2 %>%
  select(all_of(essential_1), all_of(optional_1), all_of(optional_2), all_of(essential_2), all_of(optional_3)
  )


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#5.	Pivot to 1 row for each SNP * Transcript
#6.	Select “full” table with all essential and optional columns

#Collect separate row for each distinct variable name
combined_summary_v4 <- combined_summary_v3 %>%
  pivot_wider(id_cols = c("Chromosome/scaffold name",
                          "Chromosome/scaffold position start (bp)",
                          "REF",
                          "ALT",
                          "Gene stable ID",
                          "Transcript stable ID"),
              names_from = "Variant source", values_from = "Variant name")

#get all blanks reverted to NA
#group by SNP * Transcript
# everything()
#picks up all columns except those already used for grouping
# unique(na.omit(.x))
#should boil down to 1 value for all columns unique to SNP or SNP * Transcript
# summarise(across(everything(), ~ paste(unique(na.omit(.x)), collapse = ".")))
#Variant source will be concatenated with a . in between each e.g. Exome_Capture_Diversity.Pont_2019
# left_join(combined_summary_v4
#separate row for each distinct variable name above
combined_summary_v5 <- combined_summary_v3 %>%
  mutate(across(everything(), ~na_if(.x, "-"))) %>%
  rename(Axiom = `Synonym name`) %>%
  select(-`Variant name`) %>%
  group_by(`Chromosome/scaffold name`,
           `Chromosome/scaffold position start (bp)`,
           REF,
           ALT,
           `Gene stable ID`,
           `Transcript stable ID`) %>%
  summarise(across(`Variant source`, ~ paste(unique(na.omit(.x)), collapse = ".")),
            across(everything(), ~ paste(unique(na.omit(.x)), collapse = ""))) %>%
  left_join(combined_summary_v4,
            by = c("Chromosome/scaffold name", "Chromosome/scaffold position start (bp)", "REF", "ALT", "Gene stable ID",
                   "Transcript stable ID"))

#7.	Select “intermediate” table with all essential and some optional columns
#new columns
optional_4 <- c(colnames(combined_summary_v5)[c(60:65, 28)])

#45 columns
combined_summary_v6 <- combined_summary_v5 %>%
  select(all_of(essential_1[3:12]), all_of(optional_1), all_of(optional_4), all_of(essential_2), all_of(optional_3)
  )

#8.	Select “brief” table with essential columns
#26 columns
combined_summary_v7 <- combined_summary_v5 %>%
  select(all_of(essential_1[3:12]), all_of(essential_2)
  ) %>%
  select(-`Transcript stable ID`) %>%
  group_by(`Chromosome/scaffold name`,
           `Chromosome/scaffold position start (bp)`,
           REF,
           ALT) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.x)), collapse = "")))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#WRITE
name <- "RAGT_Capture"
write_csv(variant_summary_full,
          file = paste("variant_summary_full_", name, "_", lubridate::today(),".csv", sep = ""),
          col_names = TRUE)
write_csv(variant_summary_brief,
          file = paste("variant_summary_brief_", name, "_", lubridate::today(),".csv", sep = ""),
          col_names = TRUE)


write_csv(combined_summary_v3, paste("summary_SNP_Transcript_Variantsource_", "+_", name, "_", lubridate::today(), ".csv", sep = ""))
write_csv(combined_summary_v5, paste("summary_SNP_Transcript_", "full_", "+_", name, "_", lubridate::today(), ".csv", sep = ""))
write_csv(combined_summary_v6, paste("summary_SNP_Transcript_", "intermediate_", "+_", name, "_", lubridate::today(), ".csv", sep = ""))
write_csv(combined_summary_v7, paste("summary_SNP_Transcript_", "brief_", "+_", name, "_", lubridate::today(), ".csv", sep = ""))
