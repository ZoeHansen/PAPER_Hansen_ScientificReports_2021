###############################################

# MaAsLin2 - Linear Model Construction

###############################################
### Load libraries

library(Maaslin2)

#################################################
# Loading the data 
#################################################

setwd('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Third_Analysis_ScientificReports_Submission/Resistome/MaAsLin2/')

df_input_data = read.table(file = 'campylobacter_Casecontrol_GEnorm_resistome_class.txt', header = TRUE, sep = "\t",
                           row.names = 1,
                           stringsAsFactors = FALSE)

df_input_metadata = read.table(file = 'campylobacter_metadata_casecontrol.txt', header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)


#################################################
# Subsetting the data (as needed)
#################################################

##### Subset to isolate only families 
# Metadata
df_input_metadata_fam <- subset(df_input_metadata, !is.na(New.Family.ID))
df_input_metadata_fam$New.Family.ID <- factor(df_input_metadata_fam$New.Family.ID)

# Resistome features
remove_these <- c('ER0210','ER0385','ER0557','ER0576','ER0599','ER0730','ER0775','ER0785','ER0964',
                  'ER1005','ER0236','ER0237','ER0238','ER0217','ER0218','ER0219','ER0220')

#Now we find the indicies of the rows that need to be removed
rows_to_remove <- which(row.names(df_input_data) %in% remove_these)

#And remove these corresponding rows.
df_input_data_fam <- df_input_data[-rows_to_remove,]



###### Subset Cases Only -- Case Clustering Analysis 
c.clust.meta <- read.table(file='CaseClustering/metadata_CaseClustersOnly.txt',
                           header=TRUE, sep ='\t', row.names=1, stringsAsFactors = FALSE)
c.clust.meta$Cluster <- factor(c.clust.meta$Cluster)
c.clust.meta$Residence.type <- factor(c.clust.meta$Residence.type)

keep_these <- row.names(c.clust.meta)

df_input_data_clust <- df_input_data[row.names(df_input_data) %in% keep_these,]

### Remove "Unknown" residence Type
c.clust.meta.new <- subset(c.clust.meta, Residence.type != 'Unknown')
keep_these_new <- row.names(c.clust.meta.new)

df_input_data_clust_new <- df_input_data_clust[row.names(df_input_data_clust) %in% keep_these_new,]



##### Subset Controls Only 
df_input_metadata_control <- subset(df_input_metadata, Case.status == 'Control')

keep_these <- row.names(df_input_metadata_control)

df_input_data_control <- df_input_data[row.names(df_input_data) %in% keep_these,]

### If only considering family controls
df_input_metadata_control_fam <- subset(df_input_metadata_control, New.Family.ID > 0)
keep_these_fam <- row.names(df_input_metadata_control_fam)
df_input_data_control_fam <- df_input_data_control[row.names(df_input_data_control) %in% keep_these_fam,]

df_input_metadata_control_fam$New.Family.ID <- factor(df_input_metadata_control_fam$New.Family.ID)



#################################################
# Running MaAsLin2 
#################################################

fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  output = "class_CaseStatusfixed_GenderAgeResTyperandom",  
  fixed_effects = 'Case.status',                            # include fixed effects as needed 
  random_effects = c('Gender','Age.group','Residence.type'), # include random effects as needed
  reference = 'Case.status, Control')


