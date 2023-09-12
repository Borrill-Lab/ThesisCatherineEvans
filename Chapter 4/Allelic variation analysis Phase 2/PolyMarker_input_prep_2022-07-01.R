#'by CATHERINE EVANS
#'Polymarker and VeP input preparation
#'
#'01/07/22
#'Stick together data from the command line into format for PolyMarker and VeP input
#'
#'06/07/22
#'Also make a summary table including allele frequencies
#'Also convert parts to true chromosomes - necessary for VeP to work!!
#'
#'Made in R 4.2.0 and tidyverse


require(tidyverse)

setwd("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-07-01 RAGT Exome Capture/5_flanking_regions")

name = "NAC.qc.filter.nosingletons"

coordinates <- read_delim("NAC.qc.filter.nosingletons.coordinates.bed", col_names = c("CHROM", "STARTPOS", "ENDPOS"))
locations <- read_delim("NAC.qc.filter.nosingletons.locations.txt", col_names = c("CHROM", "POS", "REF", "ALT"))
leftflank <- read_delim("NAC.qc.filter.nosingletons.leftflank.txt", col_names = c("Location", "LEFTFLANK"))
righttflank <- read_delim("NAC.qc.filter.nosingletons.rightflank.txt", col_names = c("Location", "RIGHTFLANK"))
frq <- read_delim("G:/Breeding/Wheat Inbred UK/General/Pathology/Students and work placements/Catherine Evans - JIC/NAC genetic variation 2022/2022-07-01 RAGT Exome Capture/variants_NAC_stats.frq",
                  col_names = c("CHROM",   "POS",     "N_ALLELES",       "N_CHR",   "REF:FREQ", "ALT:FREQ"),
                  skip = 1, delim = "\t")

Chinese_Spring_v1_0_pseudomolecules_parts_to_chr <- read_excel("C:/Users/CEvans/Desktop/CS/161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.xlsx", 
                                                                            col_names = c("CHROM_parts", "S_parts", "E_parts", "CHROM", "S", "E"))
#parsing issues are due to multi-allelic SNPs - i.e. this script fails for multi-allelic SNPs

#TIDY
frq2 <- frq %>%
  filter(N_ALLELES == 2) %>%
  extract(`REF:FREQ`, into = c("REF", "REF_FREQ"), regex = "([AGCT]*)\\:([0-9\\.]*)") %>%
  extract(`ALT:FREQ`, into = c("ALT", "ALT_FREQ"), regex = "([AGCT]*)\\:([0-9\\.]*)") %>%
  mutate(POS = as.numeric(POS))


#JOIN
summary <- data.frame(coordinates, locations, leftflank, righttflank) %>%
  left_join(frq2, by = c("CHROM", "POS", "REF", "ALT")) %>%
  select(CHROM, POS, STARTPOS, ENDPOS, REF, ALT, REF_FREQ, ALT_FREQ, N_CHR, LEFTFLANK, RIGHTFLANK)


#CONVERT to RefSeq v1.0 coordinates
#Parts coordinates now suffixed with _parts
summary_refseq <- summary %>%
  rename_with(~paste(.x, "_parts", sep = ""), c(CHROM, POS, STARTPOS, ENDPOS)) %>%
  left_join(Chinese_Spring_v1_0_pseudomolecules_parts_to_chr, by = "CHROM_parts") %>%
  mutate(POS = POS_parts + S, STARTPOS = STARTPOS_parts + S,ENDPOS = ENDPOS_parts + S) %>%
  select(CHROM, POS, STARTPOS, ENDPOS, REF, ALT, REF_FREQ, ALT_FREQ, N_CHR, LEFTFLANK, RIGHTFLANK, CHROM_parts, POS_parts, STARTPOS_parts, ENDPOS_parts) %>%
  extract(CHROM, into = "CHROM2", regex = "chr([[:alnum:]]*)", remove = FALSE) %>%
  mutate(Uploaded_variation = paste(CHROM2, "_", POS, "_", REF, "/", ALT, sep = ""))


#VEP
#You can tell this is older code
VarList <- summary_refseq

VarListB <- data.frame("CHROM"=rep(0,nrow(VarList)), "POS"=0, "POS2"=0, "BP"=0, "C1"=1)
VarListB["CHROM"] <- str_sub(VarList["CHROM"][[1]],4,5)
VarListB["POS"] <- VarList["POS"]
VarListB["POS2"] <- VarList["ENDPOS"]
VarListB["BP"] <- str_c(VarList["REF"][[1]], "/",VarList["ALT"][[1]])


#POLYMARKER
#Paste REF and ALT between flanking sequences
Polymarker_input <- summary_refseq %>%
  transmute(Uploaded_variation = paste(CHROM2, "_", POS, "_", REF, "/", ALT, sep = ""),
            CHROM2,
            sequence = paste(str_to_lower(LEFTFLANK), "[", REF, "/", ALT, "]", str_to_lower(RIGHTFLANK), sep = "")
  )

write_csv(summary_refseq, paste("Summary_", name, "_", lubridate::today(), ".csv", sep = ""))
write.table(VarListB, paste("VeP_input_", name, "_", lubridate::today(), ".txt", sep = ""), quote=FALSE, sep = " ", row.names=FALSE, col.names=FALSE)
write_csv(Polymarker_input, paste("Polymarker_input_", name, "_", lubridate::today(), ".csv", sep = ""))
