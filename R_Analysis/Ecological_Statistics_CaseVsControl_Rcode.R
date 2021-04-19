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

pheno <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/campylobacter_metadata_casecontrol.csv',
                 header = TRUE)

data = read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/Resistome/campylobacter_casecontrol_RPKG_resistome_gene.csv', 
                header = TRUE)

########################################################
# Selecting relevant metadata for analysis
########################################################

#Extract the Cases and Controls from our dataset and remove IDs of excluded samples
health <- pheno %>%
  dplyr::select(ER_ID, Case.status, New.Family.ID, FamID.wStatus) %>% 
  drop_na()

enviro <- pheno %>%
  dplyr::select(ER_ID, Days.since.exposed.to.case, Age.years,Gender, Residence.type)

meta <- left_join(health, enviro, by = 'ER_ID')
meta[is.na(meta)] <- 0

# Merge metadata with AGS-normalized abundances
#Full Case_Control Campylobacter dataset with metadata of interest
data_cc <- left_join(meta, data, by = "ER_ID")

# Remove any lingering zeroes to avoid downstream errors
data_cc <- data_cc[, colSums(data_cc != 0) > 0] 
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 
data_cc$Case.status <- factor(data_cc$Case.status)
data_cc$New.Family.ID <-factor(data_cc$New.Family.ID)
data_cc$FamID.wStatus <- factor(data_cc$FamID.wStatus)

########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Richness (of genes) across our samples
r.c = specnumber(data_cc[,-c(1:4)])

# Shannon diversity
h.c = diversity(data_cc[,-c(1:4)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

# Combine alpha diversity data and Case status information
div.c=tibble(health$ER_ID, health$Case.status, health$New.Family.ID, r.c, h.c, pielou.c)
colnames(div.c)=c("ER_ID", "Case.status", "FamilyID",'Richness', "Shannon", "Pielou")

### Plot alpha diversity ###
# (Figure 1)

# Reshape the data to "long" format
div.c.long=melt(div.c, id.vars=c("ER_ID", "Case.status", 'FamilyID'))

# Determine the Mean, Min, and Max for these measures
# Case status
div.means <- div.c %>%
  group_by(Case.status) %>%
  summarise(MeanRich = mean(Richness), MeanShannon = mean(Shannon), MeanPielou = mean(Pielou)) 

div.minmax <- div.c %>%
  group_by(Case.status) %>%
  summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
            MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))

rich.comp.cc <- compare_means(Richness ~ Case.status, div.c, method = 'wilcox.test')
shan.comp.cc <- compare_means(Shannon ~ Case.status, div.c, method = 'wilcox.test')
even.comp.cc <- compare_means(Evenness ~ Case.status, div.c, method = 'wilcox.test')

# Family ID
fam.div.means <- div.c %>%
  group_by(FamilyID) %>%
  summarise(MeanRich = mean(Richness), MeanShannon = mean(Shannon), MeanPielou = mean(Pielou)) #%>%

fam.div.minmax <- div.c %>%
  group_by(FamilyID)%>%
  summarise(MinRich=min(Richness), MaxRich=max(Richness), MinShannon=min(Shannon), MaxShannon=max(Shannon),
            MinPielou=min(Pielou), MaxPielou=max(Pielou))

# Statistical Comparisons
rich.comp.fam <- compare_means(Richness ~ Family.ID, div.c, method = 'wilcox.test')
shan.comp.fam <- compare_means(Shannon ~ Family.ID, div.c, method = 'wilcox.test')
even.comp.fam <- compare_means(Evenness ~ Family.ID, div.c, method = 'wilcox.test')

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
bc.d.c=vegdist(data_cc[,-c(1:8)], method="bray")


#Map desired colors to the Case statuses to make our future legend for the PCoA.  
Class=rep(1 ,nrow(data_cc))
Class[data_cc$Case.status=="Case"]= as.integer(21) #'cyan4'
Class[data_cc$Case.status=="Control"]= as.integer(24) #'darkorange3'

# Identify cases vs. controls and those who received Abx treatment:
status=rep(1,nrow(data_cc))
status[data_cc$Case.status=='Case']=21
status[data_cc$Case.status=='Control']=24
status[data_cc$Abx.use == 'Yes']=22

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


# We can add other environmental factors, but for now just leave them out
env.sub= data_cc[, c('Gender', 'Days.since.exposed.to.case','Residence.type','Age.years')]
envEF.bc=envfit(bc.c.pcoa, env.sub, permutations = 999, na.rm = TRUE)
envEF.bc

### Note: Change "rownames(envEF.bc$factors$centroids)" to alter the factor levels in plot

# Plot the ordination.
plot(bc.c.pcoa$points[,1],bc.c.pcoa$points[,2],
     cex=1.5,pch=Class,
     bg=fam,
     xlab= paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""),
     xlim=c(-0.75,0.65),
     ylim=c(-0.45, 0.45))

plot(envEF.bc, p.max=0.05, col='black', labels=list(vectors=paste('Days Since Exposed to Case')))

legend('topright', c('Case','Control','Received Antibiotics'), pch = c(21, 24, 22), 
       pt.bg=c('cyan4', 'darkorange3', 'black'), lty=0)

#Add ID labels to points (if desired)
textxy(X=bc.c.pcoa$points[,1], Y=bc.c.pcoa$points[,2],labs=data_cc$FamID.wStatus, cex=0.6)


########################################################
# PERMANOVA
########################################################
# Centroid
a.c = adonis(bc.d.c~data.cc$Case.status, distance=TRUE, permutations=1000) # Repeat for Family ID

# Disperson
b.c=betadisper(bc.d.c, group=data_cc$Case.status) #Repeat for Family ID
permutest(b.c)

# Post-hoc Tukey's Honestly Significant Difference 
TukeyHSD(b.c, which = "group", ordered = FALSE,conf.level = 0.95)

