########################################################

# Performing Ecological Statistics (Case vs Control)

########################################################
# Load libraries and data

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(vegan)
library(calibrate)
library(indicspecies)

pheno <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/Data_files_Hansen_2020/campylobacter_metadata_Hansen_2020.csv',
                 header = TRUE, na.strings = c("", 'NA'))

data = read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/Resistome_Data/Campy_fullgene_AGS_normal.csv', 
                header = TRUE)

########################################################
# Selecting relevant metadata for analysis
########################################################

#Extract the Cases and Controls from our dataset and remove IDs of excluded samples
health <- pheno %>%
  dplyr::select(ID, Case.status) %>%   
  drop_na()

# Merge metadata with AGS-normalized abundances
#Full Case_Control Campylobacter dataset with metadata of interest
data_cc <- left_join(health, data, by = "ID")

# Remove any lingering zeroes to avoid downstream errors
data_cc <- data_cc[, colSums(data_cc != 0) > 0] 
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 
data_cc$Case.status <- factor(data_cc$Case.status)

########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Richness (of genes) across our samples
r.c = specnumber(data_cc[,-c(1:2)])

# Shannon diversity
h.c = diversity(data_cc[,-c(1:2)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

# Combine alpha diversity data and Case status information
div.c=tibble(health$ID, health$Case.status, r.c, h.c, pielou.c)
colnames(div.c)=c("ID", "Case.status", "Richness", "Shannon", "Pielou")

### Plot alpha diversity ###
# (Figure 1)

# Reshape the data to "long" format
div.c.long=melt(div.c, id.vars=c("ID", "Case.status"))

# Determine the Mean, Min, and Max for these measures
div.means <- div.c %>%
  group_by(Case.status) %>%
  summarise(MeanRich = mean(Richness), MeanShannon = mean(Shannon), MeanPielou = mean(Pielou)) 

div.minmax <- div.c %>%
  group_by(Case.status) %>%
  summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
            MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))

my_comp_cc  <- list(c('Case','Control'))
colors_cc = c('cyan4', 'darkorange3')

# Plot alpha diversity values for Case Status
ggplot(data=div.c.long, aes(x=Case.status, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
  facet_wrap(~variable, scales="free_y") +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+
  labs(
    x = '\nHealth Status',
    y = 'Alpha Diversity Value\n'
  )+
  stat_compare_means(comparisons = my_comp_cc,
                     label = 'p.format')

########################################################
# Calculate and plot across-sample (Beta) diversity
########################################################

#Calculate Bray-curtis dissimilarity. 
bc.d.c=vegdist(data_cc[,-c(1:4)], method="bray")


#Map desired colors to the Case statuses to make our future legend for the PCoA.  
Class=rep(1 ,nrow(data_cc))
Class[data_cc$Case.status=="Case"]= as.integer(21) #'cyan4'
Class[data_cc$Case.status=="Control"]= as.integer(24) #'darkorange3'

# Identify the controls with Abx treatment:
Abx=rep(1, nrow(data_cc))
Abx[data_cc$Abx.use == 'Yes']=24
Abx[data_cc$Abx.use == 'No']=21

# For PCoA indicating Family ID (16 total families):
fam = rep('black', nrow(data_cc))
fam[data_cc$New.Family.ID == '1'] = 'blue'
fam[data_cc$New.Family.ID == '2'] = 'brown1'
fam[data_cc$New.Family.ID == '3'] = 'blueviolet'
fam[data_cc$New.Family.ID == '4'] = 'chartreuse'
fam[data_cc$New.Family.ID == '5'] = 'chartreuse4'
fam[data_cc$New.Family.ID == '6'] = 'coral1'
fam[data_cc$New.Family.ID == '7'] = 'cyan2'
fam[data_cc$New.Family.ID == '8'] = 'darkgoldenrod1'
fam[data_cc$New.Family.ID == '9'] = 'darkorchid4'
fam[data_cc$New.Family.ID == '10'] = 'darkslategrey'
fam[data_cc$New.Family.ID == '11'] = 'deeppink'
fam[data_cc$New.Family.ID == '12'] = 'gold4'
fam[data_cc$New.Family.ID == '13'] = 'gray'
fam[data_cc$New.Family.ID == '14'] = 'darkorange1'
fam[data_cc$New.Family.ID == '15'] = 'azure4'
fam[data_cc$New.Family.ID == '16'] = 'chocolate4'


### Perform PCoA 
bc.c.pcoa=cmdscale(bc.d.c, eig=TRUE) # for data_cc

# Calculate percent variance explained by each axes 1 and 2
ax1.v.bc.c=bc.c.pcoa$eig[1]/sum(bc.c.pcoa$eig)
ax2.v.bc.c=bc.c.pcoa$eig[2]/sum(bc.c.pcoa$eig)

# Plot the ordination.
plot(bc.c.pcoa$points[,1],bc.c.pcoa$points[,2],
     cex=1.5,pch=Class,
     bg=Res,
     xlab= paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""),
     xlim=c(-0.7,0.4),
     ylim=c(-0.5, 0.5))

legend('topright', c('Case','Control','Received Antibiotics'), pch = c(21, 21, 2), 
       pt.bg=c('cyan4', 'darkorange3', 'black'), lty=0)

#Add ID labels to points (if desired)
textxy(X=bc.c.pcoa$points[,1], Y=bc.c.pcoa$points[,2],labs=data_cc$New.Family.ID, cex=0.6)


########################################################
# PERMANOVA
########################################################
# Centroid
a.c = adonis(bc.d.c~Class, distance=TRUE, permutations=1000)

# Disperson
b.c=betadisper(bc.d.c, group=Class)

# Post-hoc Tukey's Honestly Significant Difference 
TukeyHSD(b.c, which = "group", ordered = FALSE,conf.level = 0.95)


# Performing a pairwise PERMANOVA:
# (from P. Martinez Arbizu: https://github.com/pmartinezarbizu/pairwiseAdonis)
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] 
                ~ factors[factors %in% c(as.character(co[1,elem]),(co[2,elem]))] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

PW.Adonis= pairwise.adonis(x, factors, sim.method = 'bray', p.adjust.m='bonferroni')


########################################################
# Indicator Species Analysis
########################################################

indval = multipatt(data_cc[,-c(1:6)], cluster = data_cc$Case.status, control=how(nperm=999))
summary(indval, indvalcomp = TRUE)

