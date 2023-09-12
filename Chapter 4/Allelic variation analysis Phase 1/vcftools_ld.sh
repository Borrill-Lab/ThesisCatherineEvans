#!/bin/bash
#Script by CATHERINE EVANS 2022

#Output statistics which might be useful for the subset of variants in target genes
#LD stats
#
## Error: Only one output function may be called.

set -eux

cd /home/cevans/jbaison_exomes/cevans_analysis

input_VCF=4_filtered/NAC.qc.filter.nosingletons.recode.vcf
sample_name=NAC.qc.filter.nosingletons

vcftools --vcf $input_VCF \
--geno-r2 \
--ld-window-bp 50000 \
--out ./6_linkage/${sample_name}

