#!/bin/bash
#Script by CATHERINE EVANS 2022

#These are relevant for filtering:
#Allele counts and frequencies
#Coverage depth
#Missingness
#Heterozygosity
#Site quality

#6.	Filter SNPs on read depth
#7.	Filter SNPs for quality control characteristics
#8.	Filter SNPs for <20% missing calls 

#based heavily on Scott 2021
#https://github.com/michaelfscott/DIVERSE_MAGIC_WHEAT/blob/147fb04d2da1395f8b60a8c747c340b25e0d692c/SCRIPTS/Variants/4make_posfile_from_vcfs.sh

set -eux

cd /home/cevans/jbaison_exomes/cevans_analysis

input_VCF=./2_variants/variants_NAC.vcf

samplename=NAC
min_mean_DP_param=5
max_mean_DP_param=60
min_DP_param=2
max_DP_param=120
minQ_param=100
max_missing_param=0.8 #between 0 and 1; 0 = any missing; 1 = no missing

#Biallelic sites
#Both SNPs and INDELS
#Passing DP (read depth) filter
#Passing Missing values filter

vcftools --vcf $input_VCF \
--recode --recode-INFO-all \
--min-alleles 2 --max-alleles 2 \
--min-meanDP ${min_mean_DP_param} \
--max-meanDP ${max_mean_DP_param} \
--minDP ${min_DP_param} \
--maxDP ${max_DP_param} \
--max-missing ${max_missing_param} \
--out 4_filtered/${samplename}.qc

#optional extra QC filter

bcftools view 4_filtered/${samplename}.qc.recode.vcf \
--exclude 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
--output-type v \
--output-file 4_filtered/${samplename}.qc.filter.vcf

#optional remove calls with heterozygotes

bcftools view 4_filtered/${samplename}.qc.recode.vcf \
--genotype ^het \
--output-type v \
--output-file 4_filtered/${samplename}.qc.nohet.vcf

#optional extra QC filter

bcftools view 4_filtered/${samplename}.qc.recode.vcf \
--genotype ^het \
--exclude 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
--output-type v \
--output-file 4_filtered/${samplename}.qc.nohet.filter.vcf

#optional remove singletons
vcftools --vcf 4_filtered/${samplename}.qc.recode.vcf \
--recode --recode-INFO-all \
--exclude-positions 3_quality_stats/variants_NAC_stats.singletons \
--out 4_filtered/${samplename}.qc.nosingletons

vcftools --vcf 4_filtered/${samplename}.qc.filter.vcf \
--recode --recode-INFO-all \
--exclude-positions 3_quality_stats/variants_NAC_stats.singletons \
--out 4_filtered/${samplename}.qc.filter.nosingletons

vcftools --vcf 4_filtered/${samplename}.qc.nohet.vcf \
--recode --recode-INFO-all \
--exclude-positions 3_quality_stats/variants_NAC_stats.singletons \
--out 4_filtered/${samplename}.qc.nohet.nosingletons

vcftools --vcf 4_filtered/${samplename}.qc.nohet.filter.vcf \
--recode --recode-INFO-all \
--exclude-positions 3_quality_stats/variants_NAC_stats.singletons \
--out 4_filtered/${samplename}.qc.nohet.filter.nosingletons

#################################################################
#save which sites have been removed as well
vcftools --vcf $input_VCF \
--removed-sites \
--min-alleles 2 --max-alleles 2 \
--out 4_filtered/${samplename}.biallelic

vcftools --vcf $input_VCF \
--removed-sites \
--min-meanDP ${min_mean_DP_param} \
--max-meanDP ${max_mean_DP_param} \
--minDP ${min_DP_param} \
--maxDP ${max_DP_param} \
--out 4_filtered/${samplename}.DP

vcftools --vcf $input_VCF \
--removed-sites \
--max-missing ${max_missing_param} \
--out 4_filtered/${samplename}.lmiss

bcftools view 4_filtered/${samplename}.qc.recode.vcf \
--include 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
--output-type v \
--output-file 4_filtered/${samplename}.qc.filter.removed.sites

bcftools view 4_filtered/${samplename}.qc.recode.vcf \
--genotype het \
--output-type v \
--output-file 4_filtered/${samplename}.qc.nohet.removed.sites

#################################################################
##Original version
##Abbreviated from
##https://github.com/michaelfscott/DIVERSE_MAGIC_WHEAT/blob/147fb04d2da1395f8b60a8c747c340b25e0d692c/SCRIPTS/Variants/4make_posfile_from_vcfs.sh
#
#min_mean_DP_param=5
#max_mean_DP_param=60
#min_DP_param=2
#max_DP_param=120
#minQ_param=100
#max_missing_param=1
#
#vcftools --vcf $input_VCF \
#--recode --recode-INFO-all \
#--remove-indels \
#--min-alleles 2 --max-alleles 2 \
#--min-meanDP ${min_mean_DP_param} \
#--max-meanDP ${max_mean_DP_param} \
#--minDP ${min_DP_param} \
#--maxDP ${max_DP_param} \
#--max-missing ${max_missing_param} \
#--maf 0.01 \
#--out ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs
#
#bcftools view ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf \
#--exclude 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
#--output-type v \
#--output-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.filter.vcf

##remove calls with heterozygotes
#
#bcftools view ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf \
#--genotype ^het \
#--output-type v \
#--output-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.vcf
#
#bcftools view ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.recode.vcf \
#--genotype ^het \
#--exclude 'QD<2 || FS>60.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0 || SOR>3.0 || MQ<40' \
#--output-type v \
#--output-file ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.filter.vcf
#
#grep -v "^#" ${TMPDIR}/${samplename}_chr${chr}${genome}${half}_biallelic_polymorphic.SNPs.nohet.filter.vcf |
#awk -v OFS='\t' '{print $1, $2, $4, $5}' > ${TMPDIR}/chr${chr}${genome}${half}.SNPs.txt
