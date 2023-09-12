#!/bin/bash
#Script by CATHERINE EVANS 2022

#Intersect NAC target genes with variants in NIAB founders
#Assume modules are already loaded
#Include header!!

set -e

cd /home/cevans/jbaison_exomes

bedtools intersect -wa -sorted -header -a ./Exome_Capture_RAGT_Capture/Exome_Cap_All.vcf.gz -b ./cevans_analysis/target_genes/NAC_target_genes_pts.bed > ./cevans_analysis/2_variants/variants_NAC.vcf
