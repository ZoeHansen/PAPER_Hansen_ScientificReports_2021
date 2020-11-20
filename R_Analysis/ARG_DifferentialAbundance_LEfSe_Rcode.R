########################################################

# Re-plotting Differential Abundance Data (LEfSe)

########################################################
# Load libraries and data

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)


diff.abd <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/Resistome_Data/LefseAnalysis/LDA_scores_CLASS_CaseControl.csv', 
                  header=TRUE, na.strings = c("",'NA'))

########################################################
# Prepare data for plotting
########################################################

diff.abd <- diff.abd %>%
  drop_na()
#diff.abd$LDA.Score <- as.numeric(levels(diff.abd$LDA.Score))[diff.abd$LDA.Score]
diff.abd$LDA.Score[diff.abd$Class=='Case']= -(diff.abd$LDA.Score)
diff.abd$Feature = factor(diff.abd$Feature, levels=diff.abd$Feature[order(diff.abd$LDA.Score)])

########################################################
# Plot the LDA Scores (Differential Abundance)
########################################################

### Plot the graph
### (Figure S4)
ggplot(diff.abd, aes(x=Feature, y=LDA.Score, fill=LDA.Score>0))+
  geom_col()+
  geom_hline(yintercept=0, size=1.25)+
  coord_flip()+
  scale_fill_manual(values=c("cyan4","darkorange3"),
                    labels=c("Case","Control"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = 'black',linetype = 'dashed'),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'))+
  labs(
    x = '\nARG Class\n',
    y = 'LDA Score (log10)\n',
    fill = 'Health Status')

