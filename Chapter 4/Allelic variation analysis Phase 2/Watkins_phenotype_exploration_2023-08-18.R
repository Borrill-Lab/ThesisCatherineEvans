
#exploratory watkins phenotyping

#Compare NAM-A1 haplotypes to phenotypes in Watkins core collection

#https://wisplandracepillar.jic.ac.uk/results_resources.htm

#cannot be published without permission from Simon Griffiths

#Grain and senescence data from different years. No correction for population structure

require(ggpubr)
require(tidyverse)
require(readxl)

theme_set(theme_bw())

Watkins_by_variety_2022_06_15 <- read_csv("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC genetic variation 2022/2022-06-15 Watkins data/Watkins_by_variety_2022-06-15.csv")

WISP_WGIN_Watkins_JIC_2011_grain <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC genetic variation 2022/2022-06-15 Watkins data/WISP_WGIN_Watkins_JIC_2006_2010_2011_Field_Data.xlsx", 
                                                              sheet = "2011 Watkins Core phenotypes", 
                                                              col_types = c("text", "numeric", "numeric", 
                                                                            "numeric", "numeric", "date", 
                                                                            "numeric"))

WISP_WGIN_Watkins_JIC_2010_senescence <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC genetic variation 2022/2022-06-15 Watkins data/WISP_WGIN_Watkins_JIC_2006_2010_2011_Field_Data.xlsx", 
                                                              sheet = "2010 JIC Church Farm phenotypes", 
                                                              col_types = c("text", "date", "numeric", 
                                                                            "numeric", "numeric", "text"))

JIC_Wheat_Pre_Breeding_WISP_Watkins_Core_Collection <- read_excel("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/NAC genetic variation 2022/2022-06-15 Watkins data/JIC_Wheat_Pre-Breeding_WISP_Watkins_Core_Collection.xlsx", 
                                                                  sheet = "Watkins Core Set",
                                                                  col_types = c("text", "text", "text", 
                                                                                "text", "text", "text",
                                                                                "numeric"))

Watkins_genotypes <- Watkins_by_variety_2022_06_15 %>%
  extract(Sample, "ID", "([0-9]*)")

Watkins_table <- Watkins_genotypes %>%
  left_join(WISP_WGIN_Watkins_JIC_2011_grain, by = c("ID" = "ACCESSION number")) %>%
  left_join(WISP_WGIN_Watkins_JIC_2010_senescence, by = c("ID" = "ACCESSION number")) %>%
  left_join(JIC_Wheat_Pre_Breeding_WISP_Watkins_Core_Collection, by = c("ID" = "Watkins Accession"))

Watkins_table <- Watkins_table %>%
  mutate(NAMA1_hap = fct_collapse(
    interaction(`6A_77098609`, `6A_77099397`),
    NAMA1a = "TG|TG.G|G",
    NAMA1b = "G|G.G|G",
    NAMA1c = "TG|TG.A|A",
    NAMA1d = "G|G.A|A",
    other_level = "Other")
  )

#PLOT PHENOTYPES

variable_list = c(
  "TGRWT",
  "GRSA",
  "GRWD",
  "GRLG",
  "EM",
  "PlGRYLD",
  "`heading date`",
  "`senesce flag leaf 22/7/10`",
  "`senesce stem 22/7/10`"#,
  #"lodging"
)

variable_labels <- c(
  TGRWT = "Thousand grain weight 2011",
  GRSA = "Grain area 2011",
  GRWD = "Grain width 2011",
  GRLG = "Grain length 2011",
  EM = "Ear emergence 2011",
  PlGRYLD = "Plot yield 2011",
  `heading date` = "Ear emergence 2010",
  `senesce flag leaf 22/7/10` = "senesce flag leaf 22/7/10",
  `senesce stem 22/7/10` = "senesce stem 22/7/10"#,
  #"lodging"
)


#no way...
plot_list <- list()
for(i in 1:length(variable_list)){
 boxplot <- 
   Watkins_table %>%
   filter(!is.na(NAMA1_hap)) %>%
   ggplot(aes_string(x = "NAMA1_hap", y = variable_list[i])) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0) +
   scale_x_discrete(limits = c("NAMA1a", "NAMA1b", "NAMA1c", "NAMA1d", "Other")) +
   labs(x = "", y = variable_labels[i]) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
 plot_list[[i]] <- boxplot
}

ggarrange(plotlist = plot_list, nrow = 3, ncol = 3)

save_pdf_1(ggarrange(plotlist = plot_list, nrow = 3, ncol = 3), "haplotypes_Watkins_Core_Collection.pdf", gene = "NAMA1")
save_A4_svg(ggarrange(plotlist = plot_list, nrow = 3, ncol = 3), "haplotypes_Watkins_Core_Collection.pdf", gene = "NAMA1")

x <- paste("`", colnames(Watkins_table)[2], "`", sep = "")

for(j in 2:22){
  plot_list <- list()
  for(i in 1:length(variable_list)){
    boxplot <- 
      Watkins_table %>%
      ggplot(aes_string(x = paste("`", colnames(Watkins_table)[j], "`", sep = ""),
                        y = variable_list[i])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(height = 0) +
      labs(x = "", y = variable_labels[i])
    plot_list[[i]] <- boxplot
  }
  print(ggarrange(plotlist = plot_list, nrow = 3, ncol = 3))
}


#PLOT GEOGRAPHIC DISTRIBUTION

#It's not coming up with anything super exciting

#Possible to use curated metadata from He et al as they also use Watkins lines - but did they use the same ones? Possibly not

stack_plot1 <- ggplot(Watkins_table, aes(x=fct_infreq(Origin), fill=NAMA1_hap))+
  geom_bar(position = "stack") +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(type = "seq", palette = "YlGnBu") +
  theme(axis.text.x = element_text(angle = 90))
stack_plot1

save_A5_svg(stack_plot1, paste("Watkins_haplotype_stackplot_origin_", lubridate::today(), ".svg"), gene = "NAMA1")

stack_plot1 <- ggplot(Watkins_table, aes(x=fct_infreq(`Ancestral Group`), fill=NAMA1_hap))+
  geom_bar(position = "stack") +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(type = "seq", palette = "YlGnBu") +
  theme(axis.text.x = element_text(angle = 90))
stack_plot1

save_A5_svg(stack_plot1, paste("Watkins_haplotype_stackplot_", lubridate::today(), ".svg"), gene = "NAMA1")

