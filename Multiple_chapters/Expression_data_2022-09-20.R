#'Expression data
#'
#'20/09/2022
#'
#'06/10/2022
#'Add poster theme
#'
#'27/10/2022
#'Add AAB seminar theme
#'
#'
#'Built under R 4.0.4
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#'SOURCE
#'Packages and source files
require(readxl); require(readr); require(lubridate)
require(tidyverse)
theme_set(theme_bw())

#LOAD
setwd("C:/Users/evansc/OneDrive - Norwich BioScience Institutes/DAP-seq/2022-10-06 NAC TF background")
expression <- read_delim("heatmap_20_9_2022_tpm_raw.tsv", delim = "\t", skip = 13)
target_genes <- read_excel("NAC_target_genes.xlsx")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#THEME
theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 32), plot.title = element_text(size = 32), legend.title = element_text(size=32),
                  axis.text = element_text(size = 20), legend.text = element_text(size = 20), strip.text = element_text(size = 20),
                  panel.grid = element_blank(), legend.position = "right")

subgenome_scale <- scale_colour_manual(values=c("A"="#3182bdff","B"="#62c2acff","D"="#e5f58dff"))

subgenome_scale_2 <- scale_colour_manual(values=c("A"="#0868acff","B"="#43a2caff","D"="#a8ddb5"))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#TIDY
#Make long data
expression_long <- expression %>%
  pivot_longer(cols = 13:29, names_to = "Gene stable ID", values_to = "tpm") %>%
  left_join(target_genes, by = "Gene stable ID")

#Select focus dataset
senescence <- expression_long %>%
  filter(study == "PRJNA497810") %>%
  select(Experiments, study, Age, Variety, `Gene stable ID`, `Gene name`, Triad, Subgenome, tpm) %>%
  tidyr::extract(col = Age, into = "DPA", regex = "([0-9]+)dpa", remove = FALSE)

fewer_triads <- c("NAM-1", "NAC-5", "NAP")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#plot

plot <- ggplot(senescence, aes(x = as.numeric(DPA), y = tpm, col = Subgenome, group = `Gene name`)) +
  #geom_point() +
  facet_wrap(vars(Triad)) +
  stat_summary(fun = mean, geom = "line") +
  subgenome_scale +
  labs(x = "Days after Anthesis", y = "mRNA in flag leaf (tpm)")

plot
plot + facet_wrap(vars(Triad), scale = "free_y")

ggsave(filename = paste("NAC_all_tpm_", today(), ".pdf", sep = ""), plot,
       height = 180, width = 320, units = "mm",
       device = "pdf")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#POSTER PLOT
plot <- senescence %>%
  filter(Triad %in% fewer_triads) %>%
  mutate(Triad = factor(Triad, levels =c("NAC-5", "NAP", "NAM-1"))) %>%
  ggplot(aes(x = as.numeric(DPA), y = tpm, col = Subgenome, group = `Gene name`)) +
  #geom_point() +
  facet_wrap(vars(Triad), nrow = 1) +
  stat_summary(fun = mean, geom = "line", size = 1.5) +
  subgenome_scale +
  scale_x_continuous(breaks = c(3, 7, 14, 21, 26)) +
  labs(x = "Days after Anthesis", y = str_wrap("Relative expression in flag leaf (tpm)", width = 20))

plot
plot + my_theme

ggsave(filename = paste("NAC_poster_tpm_", today(), ".pdf", sep = ""), plot + my_theme,
       height = 120, width = 370, units = "mm",
       device = "pdf")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#AAB seminar plot
plot <- senescence %>%
  #filter(Triad %in% c("NAC-5", "NAP", "NAM-1", "NAM-2")) %>%
  mutate(Triad = factor(Triad, levels =c("NAM-2", "NAC-5", "NAM-1", "NAP", "NAC3-1", "NACS"))) %>%
  ggplot(aes(x = as.numeric(DPA), y = tpm, col = Subgenome, lty = Subgenome, group = `Gene name`)) +
  #geom_point() +
  facet_wrap(vars(Triad), nrow = 2) +
  stat_summary(fun = mean, geom = "line", size = 1.5) +
  subgenome_scale_2 +
  scale_x_continuous(breaks = c(3, 7, 10, 15, 21, 26)) +
  labs(x = "Days after Anthesis", y = str_wrap("Relative expression in flag leaf (tpm)", width = 20))

plot <- plot + theme_bw(base_size = 16) + theme(panel.grid = element_blank(), legend.position = "none")
plot
plot + theme(legend.position = "right")

ggsave(filename = paste("NAC_AAB_seminar_tpm_", today(), ".svg", sep = ""), plot,
       height = 150, width = 250, units = "mm",
       device = "svg")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#AAB seminar plot individual
for(triad in unique(senescence$Triad)){
  plot <- senescence %>%
    filter(Triad == triad) %>%
    ggplot(aes(x = as.numeric(DPA), y = tpm, col = Subgenome, lty = Subgenome, group = `Gene name`)) +
    facet_wrap(vars(Triad)) +
    stat_summary(fun = mean, geom = "line", size = 1.5) +
    subgenome_scale_2 +
    scale_x_continuous(breaks = seq(0, 26, 5)) +
    labs(x = "Days after Anthesis", y = str_wrap("Relative expression in flag leaf (tpm)", width = 20))
  
  plot <- plot + theme_bw(base_size = 20) + theme(panel.grid = element_blank(), legend.position = "none")
  plot
  
  ggsave(filename = paste("NAC_AAB_seminar_tpm_", triad, "_", today(), ".svg", sep = ""), plot,
         height = 75, width = 125, units = "mm",
         device = "svg")
}
