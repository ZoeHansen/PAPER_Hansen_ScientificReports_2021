#################################################

# Case Clustering Analysis 

#################################################
# Load libraries and data

library(tidyverse)
library(ape)
library(RColorBrewer)
library(vegan)

gene = read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/Resistome_Data/Campy_fullgene_AGS_normal.csv', 
                header = TRUE)

meta <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/campylobacter_metadata.csv', 
                  header=TRUE, na.strings = c("",'NA'))

meta <- meta %>%
  filter(!grepl('\\<23\\>', ID))%>%  
  filter(!grepl('\\<66\\>', ID)) %>%
  filter(!grepl('\\<85\\>', ID))%>%
  filter(!grepl('\\<86\\>', ID)) %>%
  dplyr::select(ID,Case.status)%>% 
  filter(!grepl('FollowUp', Case.status))%>%
  drop_na() 

gene.data <- left_join(meta, gene, by = "ID")

gene.data <- gene.data[, colSums(gene.data != 0) > 0] 
gene.data <- gene.data[rowSums(gene.data != 0) > 0,] 
gene.data$Case.status <- factor(gene.data$Case.status)
gene.data$ID <- as.character(gene.data$ID)

#################################################
# Exract the data for cases to perform clustering analysis
#################################################

# The "gene.data" dataframe will be used to include controls in the the comparison of case clusters

case.data <- gene.data %>%
  filter(grepl('Case', Case.status))

# Rownames will be used as labels
case.data1 <- case.data %>%
  remove_rownames(.)%>%
  column_to_rownames(var = 'ID')

bc.dist.case <- vegdist(case.data1[,-1], method = 'bray')
case.clust <- hclust(bc.dist.case, method = 'average')


colors = c("red", "blue" , "green")
clus3 = cutree(case.clust, 3)

# Create the dendrogram
# (Figure S6)

plot(as.phylo(case.clust),
     tip.color = colors[clus3], show.tip.label = TRUE, label.offset = 0.01, cex = 1.0)

legend('topleft', legend = c('Cluster 1', 'Cluster 2','Cluster 3'), fill = c('red', 'blue','green'))


#################################################
# Add Case Cluster assignments to data
#################################################

clus3 <- as.data.frame(clus3)

clus3.2 <- rownames_to_column(clus3)
colnames(clus3.2) = c('ID','Cluster')

case.clusters <- left_join(clus3.2, gene.data, by = 'ID')

### Perform alpha-, beta- diversity analyses for case clusters with this dataframe 
### See 'Ecological_Statistics_CaseVsControl_Rcode' 



# Create a separate dataframe that includes controls for comparison

case.control.clust <- full_join(clus3.2, gene.data, by = "ID")
case.control.clust <- case.control.clust %>%
  replace_na(list(Cluster = 'control' ))

### Perform alpha-, beta- diversity analyses for case clusters vs. controls with this dataframe 
### See 'Ecological_Statistics_CaseVsControl_Rcode'
### Figure 6 - PCoA of case clusters shown with controls

