#################################################

# Total Gene Counts -- GE-normalized ARGs

#################################################
### Load libraries

library(tidyverse)
library(ggplot2)
library(ggpubr)

### Dataset Preparation
args <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/Resistome/campylobacter_casecontrol_RPKG_resistome_gene.csv',
                 header = TRUE)

args <- args[, colSums(args != 0) > 0]
args <- args[rowSums(args != 0) >0, ]


meta <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/campylobacter_metadata_casecontrol.csv', header = TRUE)

health <- meta %>%
  select(ER_ID, Case.status)%>%
  filter(Case.status != 'FollowUp')

all <- left_join(health, args, by = "ER_ID")

#####################################################
# Determining Actual Abundance of ARGs among Samples
#####################################################

sub_all <- all %>%
  mutate(Total_Count = rowSums(all[,-c(1,2)])) %>%
  mutate(Avg_Count = Total_Count/(ncol(all)-3)) %>%
  select(ER_ID, Case.status, Total_Count, Avg_Count)

group_means <- sub_all %>%
  group_by(Case.status)%>%
  summarise(AvgTotalCount = mean(Total_Count), AllTotalCount = sum(Total_Count))

compare_means(AvgTotalCount ~ Case.status, sub_all, method = 'wilcox.test')
compare_means(AllTotalCount ~ Case.status, sub_all, method = 'wilcox.test')


### Plot alpha diversity of total ARG abundances
my_comparisons = list(c('Case','Control')) 
colors_cc = c('cyan4', 'darkorange3')

ggplot(data=sub_all, aes(x=Case.status, y=Total_Count))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
  #  geom_point(aes(group=Family.ID, fill=Case.status),size=2,        
  #             position = position_dodge(0.2), shape=21)+
  #  geom_line(aes(group=Family.ID), position = position_dodge(0.2), color='gray44')+    # These comments are relevant to the family analysis
  #  facet_wrap(~variable, scales="free_y") +
  scale_color_manual(values=colors_cc)+
  annotate("label", x = 2.25, y =200, label = " Means: \n Case: 64.2 \n Control: 12.0")+
  theme_bw(base_size = 11) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+
  labs(
    x = '\nHealth Status',
    y = 'Average Total Counts Assigned to ARGs'
  )+
  stat_compare_means(comparisons = my_comparisons,
                     label = 'p.format')


#################################################
# Investigating Specific ARG Counts
#################################################

class <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/Resistome/campylobacter_Casecontrol_GEnorm_resistome_class.csv',
                  header = TRUE)

all_class <- left_join(health, class, by='ER_ID')

class_means <- all_class %>%
  group_by(Case.status) %>%
  summarise(AvgAminocoumarins=mean(Aminocoumarins),
            AvgAminoglycosides=mean(Aminoglycosides),
            AvgBacitracin=mean(Bacitracin),
            AvgBetalactams=mean(betalactams),
            AvgCAT=mean(Cationic_antimicrobial_peptides),
            AvgElfamycins=mean(Elfamycins),
            AvgFluoroquinolones=mean(Fluoroquinolones),
            AvgFosfomycin=mean(Fosfomycin),
            AvgGlycopeptides=mean(Glycopeptides),
            AvgLipopeptides=mean(Lipopeptides),
            AvgMetronidazole=mean(Metronidazole),
            AvgMLS=mean(MLS),
            AvgMDR=mean(Multi.drug_resistance),
            AvgPhenicol=mean(Phenicol),
            AvgRifampin=mean(Rifampin),
            AvgSulfonamide=mean(Sulfonamides),
            AvgTetracycline=mean(Tetracyclines),
            AvgTrimethoprim=mean(Trimethoprim))

View(class_means)

class_sum <- all_class %>%
  group_by(Case.status) %>%
  summarise(AvgAminocoumarins=sum(Aminocoumarins),
            AvgAminoglycosides=sum(Aminoglycosides),
            AvgBacitracin=sum(Bacitracin),
            AvgBetalactams=sum(betalactams),
            AvgCAT=sum(Cationic_antimicrobial_peptides),
            AvgElfamycins=sum(Elfamycins),
            AvgFluoroquinolones=sum(Fluoroquinolones),
            AvgFosfomycin=sum(Fosfomycin),
            AvgGlycopeptides=sum(Glycopeptides),
            AvgLipopeptides=sum(Lipopeptides),
            AvgMetronidazole=sum(Metronidazole),
            AvgMLS=sum(MLS),
            AvgMDR=sum(Multi.drug_resistance),
            AvgPhenicol=sum(Phenicol),
            AvgRifampin=sum(Rifampin),
            AvgSulfonamide=sum(Sulfonamides),
            AvgTetracycline=sum(Tetracyclines),
            AvgTrimethoprim=sum(Trimethoprim))

View(class_sum)

