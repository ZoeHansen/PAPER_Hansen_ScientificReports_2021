#################################################

# Heatmap of ARG Abundances 

#################################################
# Load libraries and data

library(tidyverse)
library(metagenomeSeq)

arg_data <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/Resistome/campylobacter_Casecontrol_GEnorm_resistome_group.csv',
               header = TRUE)

meta <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/ScientificReports_DataFiles/campylobacter_metadata_casecontrol_Hansen2021.csv',
                 header = TRUE)

meta <- meta %>%
  dplyr::select(ID, Case.status)%>%
  drop_na()%>%
  arrange(Case.status, ID)

arg_data <- left_join(meta, arg_data, by = 'ID')


# Remove all samples/groups in which there are exactly 0.00 counts overall
arg_data <- arg_data[, colSums(arg_data != 0) > 0]
arg_data <- arg_data[rowSums(arg_data != 0) >0, ]

# If desired: create cutoff to cleanup the heatmap and reduce the high number of zero counts
# To do this, set a cutoff value for the colSums
arg_data_cut <- arg_data[, colSums(arg_data != 0) >5]
arg_data_cut <- as.data.frame(arg_data_cut)

#################################################
# Preparing data for MetagenomeSeq
#################################################

# Making the dataframe of read counts (Samples = columns, groups = rows)
count <- arg_data_cut %>%
  remove_rownames(.) %>%
  select(-Case.status)%>% 
  column_to_rownames(., 'ID') %>%
  t(.) 

count <- as.data.frame(count)

# Using a separate file for the classifications (samples = rows)
pheno <- arg_data_cut %>%
  remove_rownames(.) %>%
  column_to_rownames(., 'ID') %>%
  select(Case.status)

pheno <- as.data.frame(pheno)


# Making a separate dataframe with the group names (match order given in the 'count' dataframe)
group.id <- count %>%
  rownames_to_column(., var = 'Group') %>%
  select(Group)%>%
  column_to_rownames(., 'Group')

# Make sure to remove the rownames in "count" AFTER you have made the group.id object. 
count <- remove_rownames(count)

#Make the MetagenomeSeq annotated dataframes:
c.pheno <- AnnotatedDataFrame(pheno)

c.group = AnnotatedDataFrame(group.id)

c.data = newMRexperiment(count, phenoData = c.pheno, featureData = c.group)

#################################################
# Constructing the heatmap
#################################################

c.trials = pData(c.data)$Case.status 
hmColColors = brewer.pal(8, "Dark2")[as.integer(factor(c.trials))] 
hmCols = colorRampPalette(brewer.pal(9, "BuPu"))(16) 

# Create the heatmap
# (Figure 4)

plotMRheatmap(obj = c.data, n = 70,
              norm = TRUE,
              log = TRUE,
              margins = c(2,4),
              xlab = 'Samples',
              ylab = 'ARG Groups',
              cexRow = 0.5, 
              labCol = FALSE,
              Colv = TRUE,
              trace = "none",
              col = hmCols,
              ColSideColors = hmColColors)

legend('topleft', inset = c(0,0.15), c('Case','Control'), pch = 22, pt.bg=c('cyan4','darkorange3'))


#################################################
# Construct a correlation plot
#################################################

# Perform the correlation test using MetagenomeSeq
cors = correlationTest(c.data, norm = TRUE, log = TRUE)

plotCorr(obj = c.data, n = 70, 
         norm = TRUE,
         log = TRUE,
         cexRow = 0.4, 
         cexCol = 0.5, 
         trace = "none", 
         dendrogram = "col", 
         col = hmCols)




