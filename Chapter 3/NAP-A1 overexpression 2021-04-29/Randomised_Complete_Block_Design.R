#'25.05.21
#'Randomised design
#'Generate a randomised design for NAP-A greenhouse experiment
#'Randomised complete block design = design where treatments are assigned to blocks;
#'each block contains constant conditions and equal number of samples of each treatment;
#'samples within a block are randomised;
#'conditions can vary between blocks;
#'and 'block' is added as a factor in the ANOVA, with 'blocks-1' degrees of freedom.
#'
#'Here I will use 6 blocks with 2 plants of each block
library(tidyverse)

genotypes = c("CTA7.con1","CTA10.con1","CTA7.5","CTA7.19","CTA10.3","CTA10.5","CTA10.9","CTA7.7","CTA7.8","CTA7.9","CTA10.1","CTA10.8")
n_plants = 144
n_blocks = 6
plants_per_block = n_plants/n_blocks
plants_per_genotype = n_plants/length(genotypes)
genotypes_per_block <- plants_per_block/length(genotypes)

design <- data.frame(unrandomised_number = 1:n_plants, block = sort(rep(1:n_blocks, plants_per_block)), genotype = rep(genotypes, plants_per_genotype), plant_num = sort(rep(1:plants_per_genotype, length(genotypes))))

set.seed(25052021) #for reproducible random-ness
randomised_design <- design %>%
  group_by(block) %>% #each group is treated like a separate table
  mutate(randomised_number = sample(unrandomised_number, replace = FALSE)) #replace=FALSE ensures each number is also picked once

#also randomise plant number as this affected the layout of the plants at seedling stage
randomised_plant_number <- randomised_design %>%
  group_by(genotype) %>%
  mutate(randomised_plant_num = sample(plant_num, replace = FALSE)) %>%
  select(randomised_number, block, genotype, randomised_plant_num) %>%
  arrange(randomised_number)

setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics")
write.csv(randomised_plant_number, "NAPA_randomised_complete_block_design_20210525.csv")
