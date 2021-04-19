#################################################

# Sequencing Metrics of Case and Control Samples

#################################################
# Load libraries
library(Nonpareil)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(outliers)

#################################################
# Assessing Coverage with Nonpareil
#################################################
# Dataset Preparation

# Nonpareil output
non_out <- read.csv('D://HPCC/Nonpareil_output_all_FINAL2.csv', sep = ',', header = TRUE)
non_out$ID <- as.character(non_out$ID)


# Metadata file used to extract relevant sample information
campy_meta <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/ScientificReports_DataFiles/campylobacter_metadata_casecontrol_Hansen2021.csv',
                 header = TRUE)


# Extract the sample IDs to match with Nonpareil output
campy_IDs <- campy_meta %>%
  select(ID)
campy_IDs$ID <- as.character(campy_IDs$ID)

df.campy <- left_join(campy_IDs, non_out, by = 'ID')


# Plot Nonpareil coverage curves
# (Figure S1)
nps<- Nonpareil.curve.batch(df.campy$File, r=df.campy$R, g=df.campy$G, b=df.campy$B, labels='', plot=TRUE, plot.opts=list(plot.observed=FALSE), star = 95)

# Show the estimated values of Nonpareil
print(nps)

# Show current coverage of the dataset (as %)
summary(nps)[,"C"]*100

# Extract Nonpareil sequence diversity (Nd) index
summary(nps)[,"diversity"]

# Extract actual sequencing effort (in Gbp)
summary(nps)[,'LR']/1e9

# Extract sequencing effort for nearly complete coverage (in Gbp)
summary(nps)[,"LRstar"]/1e9

# Predict coverage for a sequencing effort of 10Gbp
sapply(nps$np.curves, predict, 10e9)

#################################################
# Exploring Average Genome Size (AGS) and Genome Equivalents (GE)
#################################################

campy_meta <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/campylobacter_metadata_casecontrol.csv',
                 header = TRUE)

# Use the "campy_meta" datagframe from above to extract the Average Genome Size (AGS) and Genome Equivalents (GE)
campy_seq <- campy_meta %>%
  select(ER_ID, Case.status, Avg_GS, GenomeEquivalents)

# Convert to "long" format
smc_long <- melt(campy_seq, id.vars=c('ER_ID','Case.status'))

#Determine means
c_means <- campy_seq %>%
  group_by(Case.status)%>%
  summarise(MeanAGS = mean(Avg_GS), MeanGE = mean(GenomeEquivalents))

c_minmax <- campy_seq %>%
  group_by(Case.status)%>%
  summarise(MinAGS=min(Avg_GS), MaxAGS=max(Avg_GS), MinGE=min(GenomeEquivalents), MaxGE=max(GenomeEquivalents))


compare_means(Avg_GS ~ Case.status, campy_seq, method='wilcox.test')
compare_means(GenomeEquivalents ~ Case.status, campy_seq, method= 'wilcox.test')


campy_comparisons <- list(c('Case', 'Control'))


#Generate the Plot
# (Figure S2)

colors_campy = c('cyan4','darkorange3')
labels <- c(Avg_GS = 'Average Genome Size', GenomeEquivalents = "Genome Equivalents")

ggplot(data=smc_long, aes(x=Case.status, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
  facet_wrap(~variable, scales="free_y",labeller=labeller(variable=labels)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size=12),
        panel.spacing.x = unit(2, 'lines'),
        plot.margin=unit(c(1,1,0.5,0.5), 'cm'))+
  labs(
    x = 'Health Status\n',
    y = 'Alpha Diversity Value\n'
  )+
  stat_compare_means(comparisons = campy_comparisons,
                     method = 'wilcox.test',
                     label = 'p.format')


