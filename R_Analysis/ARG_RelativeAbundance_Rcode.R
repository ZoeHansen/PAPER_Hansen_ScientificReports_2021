#################################################

# Relative Abundance of ARGs

#################################################
# Load libraries and data

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(scales)

class_data <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/Resistome_Data/Campy_fullclass_AGS_normalized.csv',
                 header = TRUE)
class_data <- class_data[, colSums(class_data != 0) > 0]


# Import metadata
meta <- read.csv('D://Manning_ERIN/CampylobacterSubset_AIM_ONE/Second_Analysis/campylobacter_metadata.csv',
                 header = TRUE)
meta <- meta %>%
  select(ID, Case.status)%>%
  filter(!grepl('\\<23\\>', ID))%>% 
  filter(!grepl('\\<66\\>', ID)) %>%
  filter(!grepl('\\<85\\>', ID))%>%
  filter(!grepl('\\<86\\>', ID))%>%
  filter(!grepl('FollowUp', Case.status))%>%
  drop_na()

#################################################
# Calculate relative abundance of ARGs in samples
#################################################

# Add a column with the total reads/genome equivalents for each sample
class_data1 <- class_data %>%
  mutate(., SampTotal = rowSums(class_data[,-1]))

# Create a new dataframe with the relative abundance information
ra <- cbind(class_data1$ID, class_data1[, -c(1,19)] / class_data1$SampTotal)
colnames(ra)[1] <- 'ID'


# Merge our relative abundance dataframe with metadata and order by Case.status
ra_df <- left_join(meta, ra, by = "ID") %>%
  arrange(., Case.status)

# Rename a few of the columns for formatting (Class level)
ra_df <- rename(ra_df, "MDR" = 'Multi.drug.resistance')
ra_df <- rename(ra_df, 'Beta-lactams' = 'betalactams')
ra_df <- rename(ra_df, 'CAP' = 'Cationic.antimicrobial.peptides')

#################################################
# Prepare data for plotting
#################################################

# Extract Case status data to append later (the melt function requires only one variable to collapse upon)
ra_df_ordered <- ra_df %>%
  arrange(.,Case.status)

Health <- ra_df_ordered$Case.status

# Remove the case status data from the dataframe to prepare for melting 
ra_df.cc <- ra_df_ordered %>%
  select(-Case.status)

# Create an arbitrary sequential numbering system to avoid spacing in the x-axis of the plot
id.num <- seq(1,70, 1)

# Melt our ra_df.cc variable, and attach the Health and Num variables to the dataframe
ra_df.long <- melt(ra_df.cc, id.vars = 'ID', variable.name = 'ARG_Class')
ra_df.long$Case.status = rep(Health, times = (ncol(ra_df.cc)-1))
ra_df.long$Num = rep(id.num, times = (ncol(ra_df.cc)-1))

#################################################
# Plot relative abundance 
#################################################

# Designate a color palette (used later in our ggplot function)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Create the plot
# (Figure 3)

ggplot(data = ra_df.long, aes(x = Num, y = value, fill = ARG_Class))+
  geom_bar(stat = 'identity', width = 0.9, position = 'stack')+ 
  scale_fill_manual(values = getPalette(17), 
  guide = guide_legend(nrow=4))+
  scale_x_discrete('Num', name = 'Case Status')+
  scale_y_continuous(expand = c(0.01,0))+
  facet_wrap( ~ Case.status, strip.position = 'bottom', scales = 'free_x')+ 
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'ARG Class')+
  xlab('\nCase Status\n')+
  ylab('Relative Abundance per Sample\n')























