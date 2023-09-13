#'23/03/21
#'Make neat looking plots of NAC-5A data for poster
#'Use "light" data only and generate new linear models
#'Some of the code taken from NAC5A_detailed_analysis.rmd
#'24/03/21
#'Added Tukey post-hoc tests as a more appropriate statistical test for multiple comparison
#'Replace line names with names without "CTA".

library(readr); library(ggplot2); library(reshape2); library(dplyr)
theme_set(theme_bw())
setwd("/Users/u1984449/OneDrive - University of Warwick/NAC transgenics/NAC phenotyping")
Table <- read_csv(file = "NAC5A_phenotype_yday_190321.csv", col_names = TRUE)

#Create table with only lines in the LIGHT
#Remove missing plants
Table <- Table %>%
  filter(!is.na(Height) & Light %in% c("L", "LL"))

#Make a table for SPAD values (copied from NAC5A_graphs_050321.R)
four_lines <- c("CTA5.con1", "CTA8.con1", "CTA8.2", "CTA8.3")
Table_melted <- melt(Table, measure.vars = c("SPAD_1820DAH", "SPAD_2527DAH","SPAD_3234DAH","SPAD_3941DAH","SPAD_4648DAH"))
excl_date=53 #Exclude any SPAD values taken at 53 days after 1st Jan
Table_melted <- Table_melted %>%
  mutate(value=ifelse((variable=="SPAD_1820DAH" & First_score_date == excl_date)|(variable=="SPAD_2527DAH" & First_score_date == excl_date-7)|(variable=="SPAD_3234DAH" & First_score_date == excl_date-14)|(variable=="SPAD_3941DAH" & First_score_date == excl_date-21)|(variable=="SPAD_4648DAH" & First_score_date == excl_date-28), NA, value))
Table_SPAD <- Table_melted %>%
  filter(Line_name %in% four_lines) %>%
  dcast(formula = Line_name + Plant_no + Plant_name + Parental_genotype + Heading_date + Light + DAH + Flag_leaf_25_DAH + Peduncle_100_DAH + variable ~ "SPAD")

#Organise theme and font sizes
theme_set(theme_bw())
my_theme <- theme(axis.title = element_text(size = 24), plot.title = element_text(size = 32), axis.text = element_text(size = 16), 
                  panel.grid = element_blank(), legend.position = "right", legend.text = element_text(size = 16),
                  legend.title = element_text(size=24))

#Labels etc.
four_lines <- c("CTA5.con1", "CTA8.con1", "CTA8.2", "CTA8.3")
SPAD_labs <- c(SPAD_1820DAH="18-20", SPAD_2527DAH="25-27",SPAD_3234DAH="32-34",SPAD_3941DAH="39-41",SPAD_4648DAH="46-48")
genotype_labs <- c(CTA5.con1="5.con", CTA8.con1="8.con", CTA5.4="5.4", CTA8.9="8.9", CTA8.26="8.26", CTA8.2="8.2", CTA8.3="8.3", CTA8.1 = "8.1", CTA8.10 = "8.10", CTA8.23 = "8.23", CTA8.25 = "8.25", CTA5.5 = "5.5", CTA5.1 = "5.1", CTA5.2="5.2")
genotype_order <- c("CTA5.con1","CTA8.con1", "CTA8.1","CTA8.2","CTA8.3","CTA8.10")
genotype_scale <- scale_fill_manual(values=c("control"="#ca0020","single"="#92c5de","multi"="#0571b0"))
line_scale <- scale_color_manual(values=c("CTA5.con1"="#ca0020","CTA8.con1" = "#f4a582", "CTA8.2" ="#92c5de","CTA8.3"="#0571b0"), limits = c("CTA5.con1", "CTA8.con1", "CTA8.2", "CTA8.3"), labels = genotype_labs)
Heading_labs<- c("4"="04/01","8"="08/01","12"="12/01","16"="16/01")
#This colour scale is based on ColorBrewer palette RdBu with 4 classes 
#https://colorbrewer2.org/?type=diverging&scheme=RdBu&n=4

par(mfrow=c(2,2))
#Plot heading date
heading_plot <- ggplot(Table, aes(x=Line_name, y=Heading_date, fill =Parental_genotype)) +
  geom_boxplot() +
  scale_x_discrete(limits = genotype_order, labels = genotype_labs) +
  scale_y_continuous(breaks = seq(4,16, 4), labels = Heading_labs) +
  my_theme +
  genotype_scale +
  theme(legend.position = "none") +
  labs(x = "Line", y = "Heading date")
heading_plot

#Plot leaf senescence
plot_F <- ggplot(Table, aes(x=Line_name, y=Flag_leaf_25_DAH, fill = Parental_genotype)) +
  geom_boxplot() +
  scale_x_discrete(limits = genotype_order, labels = genotype_labs) +
  my_theme +
  genotype_scale +
  theme(legend.position = "none") +
  labs(x = "Line", y = "25% Flag leaf senescence (DAH)")
plot_F
#Plot peduncle senescence
plot_P <- ggplot(Table, aes(x=Line_name, y=Peduncle_100_DAH, fill = Parental_genotype)) +
  geom_boxplot() +
  scale_x_discrete(limits = genotype_order, labels = genotype_labs) +
  my_theme +
  genotype_scale +
  theme(legend.position = "none") +
  labs(x = "Line", y = "100% Peduncle senescence (DAH)")
plot_P

#Do the statistical tests
lm_Heading_date <- lm(Heading_date ~ Line_name, data = Table)
summary(lm_Heading_date)
anova(lm_Heading_date)
par(mfrow=c(1,2))
plot(lm_Heading_date, which=c(1,2)) #Hmm, it doesn't look as robust as the previous, probably because it's all near 0....

lm_Flag_leaf_25_DAH <- lm(Flag_leaf_25_DAH ~ Line_name, data = Table)
summary(lm_Flag_leaf_25_DAH)
anova(lm_Flag_leaf_25_DAH)
par(mfrow=c(1,2))
plot(lm_Flag_leaf_25_DAH, which=c(1,2))

lm_Peduncle_100_DAH <- lm(Peduncle_100_DAH ~ Line_name, data = Table)
summary(lm_Peduncle_100_DAH)
anova(lm_Peduncle_100_DAH)
par(mfrow=c(1,2))
plot(lm_Peduncle_100_DAH, which=c(1,2))

#Plot SPAD
SPAD_plot <- ggplot(filter(Table_SPAD, variable != "SPAD_4648DAH"), aes(x=variable, y=SPAD, group = Line_name, shape=Line_name, color=Line_name))
SPAD_plot <- SPAD_plot +
  stat_summary(fun.data = "mean_se", geom=("errorbar"), width=0.1) +
  stat_summary(fun.data = "mean_se", geom="line") +
  stat_summary(fun.data = "mean_se", geom="point", size=3) +
  labs(x="Days after Heading", y="SPAD", shape = "Line", color="Line") +
  scale_x_discrete(labels=SPAD_labs) +
  my_theme
SPAD_plot + line_scale +
  scale_shape(limits = c("CTA5.con1", "CTA8.con1", "CTA8.2", "CTA8.3"), labels = genotype_labs) +
  theme(legend.position = "none")

#Can we do any statistical tests for the SPAD values?
lm_SPAD_together <- lm(SPAD ~ variable * Parental_genotype, data = Table_SPAD)
summary(lm_SPAD_together)
anova(lm_SPAD_together)
plot(lm_SPAD_together, which=c(1,2))


#I'm really not sure about the heading date comparisons so I'm going to validate them with a non-parametric test
#Since they are really on an arbitrary scale from 1st January
#Taken from Stats4_MixedEffects.R
pValues.wilcox   <- numeric(5)
pValues.wilcox[] <- NA

#Comparing all against CTA5.con1, find those heading dates
CTA5con1 <- Table %>%
  filter(Line_name == "CTA5.con1") %>%
  select(Heading_date)

for(i in 1:5){
  compare <- Table %>%
    filter(Line_name == genotype_order[(i+1)]) %>%
    select(Heading_date)
  outputObject      <- wilcox.test(unlist(compare), unlist(CTA5con1), exact=FALSE) #unlist to make numeric
  pValues.wilcox[i] <- outputObject$p.value
}
cbind(genotype_order[2:6], pValues.wilcox*5) #Bonferroni correction: multiply by 5 'cause we made 5 comparisons


#TUKEY test for post-hoc comparsion with ANOVA.
aov_H <- aov(lm_Heading_date)
aov_F <- aov(lm_Flag_leaf_25_DAH)
aov_P <- aov(lm_Peduncle_100_DAH)

summary(Heading_aov) #This gives us the same output as calling anova()
TUKEY_H <- TukeyHSD(x=aov_H, conf.level=0.95)
plot(TUKEY_H , las=1 , col="brown")
TUKEY_F <- TukeyHSD(x=aov_F, conf.level=0.95)
plot(TUKEY_F , las=1 , col="brown")
TUKEY_P <- TukeyHSD(x=aov_P, conf.level=0.95)
plot(TUKEY_P , las=1 , col="brown")

#Function credits to R GRAPH GALLERY
# I need to group the treatments that are not different each other together.
library(multcompView)
generate_label_df <- function(TUKEY, variable){
  require(multcompView)
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  print(Tukey.levels)
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY_H , "Line_name")
#Add these letters to the top of the boxplot
heading_plot + geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = 16, fill = "white", size = 6, label.size = 0)

LABELS <- generate_label_df(TUKEY_F , "Line_name")
plot_F + geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = 40, fill = "white", size = 6, label.size = 0)

LABELS <- generate_label_df(TUKEY_P , "Line_name")
plot_P + geom_label(data = LABELS, label = LABELS[,1], x = LABELS[,2], y = 50, fill = "white", size = 6, label.size = 0) +
  ylim(28, 50)
