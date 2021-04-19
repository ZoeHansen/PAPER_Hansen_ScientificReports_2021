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
  drop_na()%>%
  filter(LDA_score >4.0) #for 'group' level only 

diff.abd$LDA_score <- as.numeric(diff.abd$LDA_score)
diff.abd$LDA_score <- with(diff.abd, ifelse(Class == 'Case', -LDA_score, LDA_score))
diff.abd1 <- diff.abd[order(diff.abd$LDA_score, decreasing=TRUE),]
diff.abd1$Feature=factor(diff.abd1$Feature, levels=diff.abd1$Feature)

########################################################
# Plot the LDA Scores (Differential Abundance)
########################################################

### Plot the graph
### (Figures S3 and S4)
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


########################################################
# Plot LEfSe for Case Clusters
########################################################
#### Case clustering preparation 

diff.abd <- diff.abd[order(diff.abd$Class, diff.abd$LDA_score, decreasing=TRUE),]
diff.abd$Feature=factor(diff.abd$Feature, levels=diff.abd$Feature)

### Plot the graph
### (Figures S6 and S7)

ggplot(diff.abd, aes(x=Feature, y=LDA_score, fill=Class))+
  geom_col()+
  geom_hline(yintercept=0, size=1.25)+
  coord_flip()+
  scale_fill_manual(values=c("red2","blue2","green2"),
                    labels=c("Cluster 1","Cluster 2", "Cluster 3"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = 'black',linetype = 'dashed'),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'))+
  labs(
    x = '\nARG Group\n',
    y = 'LDA Score (log10)\n',
    fill = 'Case Cluster')
