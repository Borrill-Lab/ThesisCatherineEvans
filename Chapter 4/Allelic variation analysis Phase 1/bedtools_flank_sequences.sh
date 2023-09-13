#!/bin/bash
#Script by CATHERINE EVANS 2022

#Find 100bp flanking sequence left and right of each SNP location

#Assume modules are already loaded
#Include header!!

set -eux

cd /home/cevans/jbaison_exomes/cevans_analysis

input_VCF=4_filtered/NAC.qc.filter.nosingletons.recode.vcf
sample_name=NAC.qc.filter.nosingletons
flank=100
genome_FASTA=/home/cevans/jbaison_exomes/cevans_analysis/target_genes/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta

#Convert VCF to BED (yes this might not be the best method)
bedtools intersect -sorted -a target_genes/NAC_target_genes_pts.bed -b $input_VCF  > 5_flanking_regions/${sample_name}.coordinates.bed

#Extract REF and ALT alleles
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $input_VCF > 5_flanking_regions/${sample_name}.locations.txt

#Define flanking region
bedtools flank -i 5_flanking_regions/${sample_name}.coordinates.bed -g target_genes/chinese_spring_genome_parts.txt -l $flank -r 0 > 5_flanking_regions/${sample_name}.leftflank.bed
bedtools flank -i 5_flanking_regions/${sample_name}.coordinates.bed -g target_genes/chinese_spring_genome_parts.txt -r $flank -l 0 > 5_flanking_regions/${sample_name}.rightflank.bed

#Get fasta sequence for flanking region from reference genome
bedtools getfasta -fi $genome_FASTA -bed 5_flanking_regions/${sample_name}.leftflank.bed -tab > 5_flanking_regions/${sample_name}.leftflank.txt
bedtools getfasta -fi $genome_FASTA -bed 5_flanking_regions/${sample_name}.rightflank.bed -tab > 5_flanking_regions/${sample_name}.rightflank.txt
