# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#'CATHERINE EVANS
#'Functions for Haplotype Pipeline
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
#'Make functions
#'
#'13/07/2022
#'Split function script from run script
#'make impute_haplotypes_multi() to deal with sets of multiple independent genes or regions.
#'
#'Built in R 4.2.0 and tidyverse 1.3.1
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#'so
#'we've somehow selected a subset of variants
#'and we want to see what the gene-wide haplotype is in each plant
#'NOT phased by parents, but assuming a homozygous call


#'Takes a table with
#'1 row per SNP, expressed as CHROM, POS, REF, ALT *OR* as Uploaded_variation if named_snps = TRUE
#'1 col per variety
#'Returns a table with
#'1 row per variety
#'1 col per SNP, expressed as CHROM_POS_REF_ALT *OR* as Uploaded_variation if named_snps = TRUE
#'Prints a summary of the genotype frequency, considering Hom calls only
#'@param input table with 1 row per SNP, e.g. VCF format
#'@return table with 1 row per variety
#'CATHERINE EVANS
row_per_variety <- function(input, named_snps = FALSE){
  
  require(dplyr); require(tidyr)
  if(isTRUE(named_snps)){
    long_data = input %>%
      pivot_longer(cols = !(Uploaded_variation), names_to = "Variety", values_to = "Genotype") %>%
      mutate(Genotype_hom = ifelse(substr(Genotype, 1, 3) %in% c("0|0","0/0"), 0, ifelse(substr(Genotype, 1, 3) %in% c("1|1","1/1"), 1, NA))
      )
    
    #Check allele frequencies
    #doesn't work unless coded as 0 / 1
    allele_summary <- long_data %>%
      
      group_by(Uploaded_variation) %>%
      summarise(alt_freq_hom = mean(na.omit(Genotype_hom)))
    
    print("Genotype summary")
    print(allele_summary)
    
    #haplotype collapse
    
    #these are the columns we want
    snps <- unlist(input$Uploaded_variation)
    
    by_variety <- long_data %>%
      pivot_wider(id_cols = "Variety", names_from = c("Uploaded_variation"), names_sep = "_", values_from = Genotype_hom) %>%
      unite(haplotype, all_of(snps), remove = FALSE, sep = "_", na.rm = FALSE)
    
    return(by_variety)
  }else{
  long_data = input %>%
    pivot_longer(cols = !(c(CHROM, POS, REF, ALT)), names_to = "Variety", values_to = "Genotype") %>%
    mutate(Genotype_hom = ifelse(startsWith(Genotype, "0|0"), 0, ifelse(startsWith(Genotype, "1|1"), 1, NA))
    )
  
  #Check allele frequencies
  #doesn't work unless coded as 0 / 1
  allele_summary <- long_data %>%
    
    group_by(CHROM, POS, REF, ALT) %>%
    summarise(alt_freq_hom = mean(na.omit(Genotype_hom)))
  
  print("Genotype summary")
  print(allele_summary)
  
  #haplotype collapse
  
  #these are the columns we want
  snps <- paste(input$CHROM, input$POS, input$REF, input$ALT, sep  ="_") 
  
  by_variety <- long_data %>%
    pivot_wider(id_cols = "Variety", names_from = c("CHROM", "POS", "REF", "ALT"), names_sep = "_", values_from = Genotype_hom) %>%
    unite(haplotype, all_of(snps), remove = FALSE, sep = "_", na.rm = FALSE)
  
  return(by_variety)
  }
}

#'Name haplotypes
#'Based on rows without any missing values
#'Order so that hap_a is the most common
#'@param by_variety table with 1 row per variety
#'@return key with 1 row per haplotype, and haplotype names
#'CATHERINE EVANS
name_haplotypes <- function(by_variety){
  key_1 <- by_variety %>%
    group_by(haplotype) %>%
    summarise(n = n(), p = n/nrow(by_variety)) 
  
  key_3 <- by_variety %>%
    select(-Variety) %>%
    unique
  
  key_3 <- key_3 %>%
    filter(!if_any(1:ncol(key_3), is.na))
  
  key_4 <- key_1 %>%
    inner_join(key_3, by = "haplotype") %>%
    arrange(-n) %>%
    cbind(haplotype_name = paste("hap_", letters[seq(1, nrow(key_3), 1)], sep = "")) %>%
    select(haplotype, haplotype_name, n, p, everything())
  
  return(key_4)
}

#'Impute possible haplotypes
#'Work out which haplotypes rows WITH MISSING VALUES could represent
#'To compare the impact of lacking certain genotypes, replace with missing values before this step
#'@param by_variety table with 1 row per variety
#'@param key key with 1 row per haplotype, and haplotype names
#'@return full key including all combinations of missing values AND by_variety table with a "haplist" column
#'CATHERINE EVANS
impute_haplotypes <- function(key, by_variety){
  
  key_1 <- by_variety %>%
    group_by(haplotype) %>%
    summarise(n = n(), p = n/nrow(by_variety)) 
  
  key_2 <- by_variety %>%
    select(-Variety) %>%
    unique
  
  print("key and key_2")
  print(key)
  print(key_2)
  
  mat <- matrix(nrow = length(1:nrow(key_2)), ncol = length(1:nrow(key)))
  
  #impute possible haplotypes for sets with missing values
  for(i in 1:nrow(key_2)){
    for(j in 1:nrow(key)){
      #print(c(i, j))
      #print(all(key_2[i,-1] == key[j,-1:-4]))
      mat[i,j] <- !(isFALSE(all(key_2[i, -1] == key[j, -1:-4])))
    }
  }
  colnames(mat) <- key$haplotype_name
  
  #cur_column() returns column name
  
  #add stuff together
  key_5 <- cbind(key_2, mat) %>%
    mutate(across(starts_with("hap_"), ~ ifelse(.x == TRUE, cur_column(), NA))
    ) %>%
    unite(haplist, all_of(key$haplotype_name), remove = FALSE, sep = ".", na.rm = TRUE) %>%
    left_join(key_1) %>%
    select(haplotype, haplist, n, p, everything())
  
  print("Key imputed")
  print(key_5)
  
  #report
  by_variety_2 <- by_variety %>%
    left_join(key_5[c("haplotype", "haplist")], by = "haplotype") %>%
    mutate(haplist = as.factor(haplist))
  
  print("Haplotype summary")
  print(summary(by_variety_2$haplist))
  
  return(list(key_5, by_variety_2))
}


#'A wrapper for impute_haplotypes()
#'Runs impute_haplotypes in a loop for every gene in genes
#'@param by_variety_list list of tables with 1 row per variety
#'@param key_list list of keys with 1 row per haplotype, and haplotype names
#'@param genes vector of genes in gene.name column of by_variety_list and key_list
#'@param snp_table table with Uploaded_variation, gene.name and include
#'@return list of full keys including all combinations of missing values AND by_variety table with a "haplist" column and a "gene.name" column
#'CATHERINE EVANS
impute_haplotypes_multi <- function(key_list, by_variety_list, genes, snp_table){
  
  #Set empty variables
  key_imputed_list <- list()
  by_variety_summary <- c()
  
  #Loop through genes
  for(i in 1:length(genes)){
    
    by_variety <- by_variety_list[[i]] %>%
      select(-gene.name)
    
    #exclude a subset of snps to see what happens
    gene_snp_table <- snp_table %>%
      filter(gene.name == genes[i])
    
    snps <- gene_snp_table %>%
      select(Uploaded_variation) %>% unlist() %>% as.vector()
    
    exclusions <- gene_snp_table %>%
      filter(!include) %>%
      select(Uploaded_variation) %>% unlist() %>% as.vector()
    
    by_variety <- by_variety %>%
      mutate(across(all_of(exclusions), ~ NA)) %>%
      unite(haplotype, all_of(snps), remove = FALSE, sep = "_", na.rm = FALSE) %>%
      select(Variety, haplotype, everything())
    
    #run impute_haplotypes()
    output_list <- impute_haplotypes(key_list[[i]][1:ncol(key_list[[i]])-1], by_variety)
    
    key_imputed <- output_list[[1]]
    by_variety_imputed <- output_list[[2]]
    
    #Add gene name
    key_imputed <- mutate(key_imputed, gene.name = genes[i])
    by_variety_imputed <- mutate(by_variety_imputed, gene.name = genes[i])
    
    #Compile lists
    key_imputed_list <- c(key_imputed_list, list(key_imputed))
    
    by_variety_imputed_v2 <- by_variety_imputed %>%
      select(Variety, haplist, haplotype, gene.name)
    by_variety_summary <- rbind(by_variety_summary, by_variety_imputed_v2)
  }
  
  return(list(key_imputed_list, by_variety_summary))
}

