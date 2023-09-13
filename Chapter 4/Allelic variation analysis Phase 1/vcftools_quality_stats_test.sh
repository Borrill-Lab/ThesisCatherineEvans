#!/bin/bash
#Script by CATHERINE EVANS 2022

#Output statistics which might be useful for the subset of variants in target genes
#maximum output options
set -eux

cd /home/cevans/jbaison_exomes/cevans_analysis

input_VCF=2_variants/variants_NAC.vcf
echo $input_VCF

vcftools --vcf $input_VCF