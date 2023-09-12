#!/bin/bash
#Script by CATHERINE EVANS 2022

#Output statistics which might be useful for the subset of variants in target genes
#All the stats that could possibly be useful

#These are relevant for filtering:
#Allele counts and frequencies
#Coverage depth
#Missingness
#Heterozygosity
#Site quality

#These are relevant for downstream analysis or of general interest
#Allele counts and frequencies
#Singletons
#Haplotype counts per gene
#
## Error: Only one output function may be called.

set -eux

cd /home/cevans/jbaison_exomes/cevans_analysis

input_VCF=./2_variants/variants_NAC.vcf

vcftools --vcf $input_VCF \
--freq \
--out ./3_quality_stats/variants_NAC_stats

vcftools --vcf $input_VCF \
--depth \
--out ./3_quality_stats/variants_NAC_stats

vcftools --vcf $input_VCF \
--site-mean-depth \
--out ./3_quality_stats/variants_NAC_stats

vcftools --vcf $input_VCF \
--missing-indv \
--out ./3_quality_stats/variants_NAC_stats

vcftools --vcf $input_VCF \
--missing-site \
--out ./3_quality_stats/variants_NAC_stats

vcftools --vcf $input_VCF \
--singletons \
--out ./3_quality_stats/variants_NAC_stats

vcftools --vcf $input_VCF \
--get-INFO QD --get-INFO FS --get-INFO MQRankSum --get-INFO ReadPosRankSum --get-INFO SOR --get-INFO MQ \
--out ./3_quality_stats/variants_NAC_stats.qc

vcftools --vcf $input_VCF \
--get-INFO AN --get-INFO DP \
--out ./3_quality_stats/variants_NAC_stats


