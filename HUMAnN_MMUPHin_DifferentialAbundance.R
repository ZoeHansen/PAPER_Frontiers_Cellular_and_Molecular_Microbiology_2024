#################################################

# Aim 3: MMUPHin - Exploring Differentially Abundant HUMAnN Pathways

#################################################

# Load libraries
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(viridis)
library(vegan)

# Load data (MMUPHin requires relative abundance/proportions)

# Note: metadata should contain samples as rows and feature data should contain samples as columns

# Pathway abundances
abd.data <- read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/ERIN_pathabundance_relabd_merged_CaseFollowPairs.csv',
                     header=TRUE)

abd.data <- read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_PHYLUM_results/ERIN_pathabundance_phylum_relabd_CaseFollowPairs.csv',
                     header=TRUE)

abd.data <- abd.data %>%
  dplyr::filter(!grepl('Case.status', ER_ID))%>%
  dplyr::filter(!grepl('Case.Follow_ID', ER_ID)) %>%
  dplyr::filter(!grepl('Pathogen', ER_ID))# %>%
#  dplyr::filter(!grepl('UNMAPPED', ER_ID)) %>% # If we want to remove "UNMAPPED" to actually observe pathway differences
#  dplyr::filter(!grepl('UNINTEGRATED', ER_ID)) # Remove at PHYLUM level
  
# Since we are interested in extracting community total abundances for pathways (not necessarily those
# associated with a specific taxon), we will only use rows that do not contain a " | taxon-id" delimiter

abd.data.com <- abd.data %>%
  dplyr::filter(!grepl('g_', ER_ID))%>%
  dplyr::filter(!grepl('unclassified',ER_ID))

com.paths <- abd.data.com$ER_ID

abd.data.com <- as.data.frame(lapply(abd.data.com[,c(2:ncol(abd.data.com))], as.numeric))

rownames(abd.data.com) <- com.paths

abd.data.com <- abd.data.com[,colSums(abd.data.com !=0)>0]
abd.data.com <- abd.data.com[rowSums(abd.data.com !=0) >0,]

# Make sure to sort ER_ID by descending before importing (*insert eye roll here*)
meta.data <- read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/ERIN_MetaboliteAnalysis_Metadata_CaseFollowPairs_clean.csv',
                      header=TRUE)

meta.data <- meta.data %>%
  dplyr::select('ER_ID','Case.status','Pathogen','Run','Case.Follow_ID',
                'Followup.days','Residence.type','Age.years','Hospital',
                'Avg_GS','GenomeEquivalents','Antibiotics', 'Gender','Year')

rownames(meta.data) <- meta.data$ER_ID

meta.data$Antibiotics[meta.data$Antibiotics == ""] <- "No"
meta.data$Residence.type[meta.data$Residence.type == ''] <- 'Unknown'
meta.data$Hospital[meta.data$Hospital == ''] <- "Unknown"

meta.data$Run <- factor(meta.data$Run, levels=c('1','2','3','4'))
meta.data$Year <- as.numeric(meta.data$Year)
meta.data$Followup.days <- as.integer(meta.data$Followup.days)
meta.data$Case.status <- factor(meta.data$Case.status, levels=c('Case','FollowUp'))
meta.data$Antibiotics <- factor(meta.data$Antibiotics, levels=c('No','Yes'))


#### Batch Correction ####

# First, we will attempt to correct for batch effects due to Sequencing Run (Microbiome Data)

fit_adjust_batch <- adjust_batch(feature_abd=abd.data.com,    # abundance matrix
                                 batch='Run',             # the batch variable requiring correction
                                 covariates='Case.status',# covariates to adjust for (i.e. Case vs. Follow-up status)
                                 data=meta.data,          # metadata file
                                 control=list(verbose=FALSE)) # other specifiers, if desired

# Isolate adjusted abundance table
path_abd_adj <- fit_adjust_batch$feature_abd_adj

# Perform PERMANOVA to examine variability due to Sequencing Run before and after adjustment

dist_before <- vegdist(t(abd.data.com), method='bray')
dist_after <- vegdist(t(path_abd_adj), method='bray')

fit_adonis_before <- adonis(dist_before ~ Run, meta.data)
fit_adonis_after <- adonis(dist_after ~ Run, meta.data)

print(fit_adonis_before)
print(fit_adonis_after)

# Note: compare amount of variability explained using the R2 measure; check for significance
# This may depend on the resolution of taxa (e.g. phylum vs. genus vs. species)



#### Meta-analytical Differential Abundance ####

# 'lm_meta()' combines the functionality of MaAsLin2 and vegan to retrieve differentially abundant
# features between groups

# General setup of model
fit_lm_meta <- lm_meta(feature_abd=abd.data.com,       # designate abundance data
                       batch='Run',                # list the variable to correct for (batch)
                       exposure='Case.status',     # group variable
                       covariates=c(''),           # Other environmental variables to account for in the model
                       data=meta.data,             # metadata file
                       control=list(verbose=FALSE))# other specifications (if desired)

#meta.data$Run <- as.numeric(meta.data$Run)

# Run the model
fit_lm_meta <- lm_meta(feature_abd=abd.data.com,       
                       batch='Run',                
                       exposure='Case.status',     
                       covariates=c('GenomeEquivalents','Gender',
                                    'Antibiotics','Age.years'),           
                       data=meta.data,             
                       control=list(verbose=FALSE))


# View the 'meta_fits' output
meta_fits <- fit_lm_meta$meta_fits

meta_fits_signif <- meta_fits %>%
  filter(qval.fdr < 0.05)         # pull out significant p-adjusted values

write.csv(meta_fits_signif, 'D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_PHYLUM_results/MMUPHin_DifferentialAbundance/MMUPHin_DifferentialAbundance_MetaCycPathways_UNMAPPED_UNINTEGRATEDremoved_CaseStatus_CaseFollowPairs.csv',
          row.names=FALSE)

# Plot the significant findings (differential abundance)
meta_fits_signif %>%
  filter(abs(coef) > 0.035) %>%
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature, fill=coef>0)) +
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_manual(values=c("cyan4","darkorchid3"),
                    labels=c("Case","FollowUp"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = 'black',linetype = 'dashed'),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'))+
  labs(
    x = '\nMetaCyc Pathway\n',
    y = '\nCoefficient\n',
    fill = 'Health Status')


##### Investigating Discrete Population Structure #####

# For this portion, we will investigate cases and follow-ups separately

case.meta <- subset(meta.data, Case.status == 'Case')
case.abd.adj <- micro_abd_adj[, rownames(case.meta)]

D_case <- vegdist(t(case.abd.adj), method='bray')


follow.meta <- subset(meta.data, Case.status == 'FollowUp')
follow.abd.adj <- micro_abd_adj[, rownames(follow.meta)]

D_follow <- vegdist(t(follow.abd.adj), method='bray')

##### Investigating Continuous Population Structure #####

# Identify variables driving continuous structure

###### Cases & Follow-Ups Separate ######
### Cases ### 

case.fit_continuous <- continuous_discover(feature_abd = case.abd.adj,
                                           batch = "Run",
                                           data = case.meta,
                                           control = list(var_perc_cutoff = 0.5,
                                                          verbose = FALSE))

case.loading <- data.frame(feature = rownames(case.fit_continuous$consensus_loadings),
                           loading1 = case.fit_continuous$consensus_loadings[, 1])

case.loading %>%
  arrange(-abs(loading1)) %>%
  slice(1:20) %>%
  arrange(loading1) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  ggplot(aes(x = feature, y = loading1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Features with top loadings")

case.mds <- cmdscale(d = D_case)
colnames(case.mds) <- c("Axis1", "Axis2")
as.data.frame(case.mds) %>% 
  mutate(score1 = case.fit_continuous$consensus_scores[, 1]) %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = score1)) +
  geom_point() +
  coord_fixed()



### Follow-Ups ###
follow.fit_continuous <- continuous_discover(feature_abd = follow.abd.adj,
                                             batch = "Run",
                                             data = follow.meta,
                                             control = list(var_perc_cutoff = 0.5,
                                                            verbose = FALSE))

follow.loading <- data.frame(feature = rownames(follow.fit_continuous$consensus_loadings),
                             loading1 = follow.fit_continuous$consensus_loadings[, 1])

follow.loading %>%
  arrange(-abs(loading1)) %>%
  slice(1:20) %>%
  arrange(loading1) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  ggplot(aes(x = feature, y = loading1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Features with top loadings")

follow.mds <- cmdscale(d = D_follow)
colnames(follow.mds) <- c("Axis1", "Axis2")
as.data.frame(follow.mds) %>% 
  mutate(score1 = follow.fit_continuous$consensus_scores[, 1]) %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = score1)) +
  geom_point() +
  coord_fixed()



##### Cases + Follow-Ups Together ######
cf.fit_continuous <- continuous_discover(feature_abd = path_abd_adj,
                                         batch = "Run",
                                         data = meta.data,
                                         control = list(var_perc_cutoff = 0.5,
                                                        verbose = FALSE))

cf.loading <- data.frame(feature = rownames(cf.fit_continuous$consensus_loadings),
                         loading1 = cf.fit_continuous$consensus_loadings[, 1])

shape=rep(1, nrow(meta.data))
shape[meta.data$Case.status=='Case']='Case'
shape[meta.data$Case.status=='FollowUp']='FollowUp'
shape[meta.data$Antibiotics=='Yes']='Received Antibiotics'

cf.loading.data <- cf.loading %>%
  arrange(-abs(loading1)) %>%
  slice(1:30) %>%
  arrange(loading1) %>%
  mutate(feature = factor(feature, levels = feature))

write.csv(cf.loading.data, 'D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_PHYLUM_results/MMUPHin_DifferentialAbundance/MMUPHin_MetaCycPathways_Top30LoadingData_UNMAPPED_UNINTEGRATEDremoved_CaseStatus_CaseFollowPairs.csv',
          row.names=FALSE)

ggplot(data=cf.loading.data, aes(x = feature, y = loading1, fill=loading1>0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values=c("azure4","azure3"),
                    labels=c("Case-like","FollowUp-like"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = 'black',linetype = 'dashed'),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.5,0.25,0,0.25), 'cm'),
        plot.title=element_text(face='bold',hjust=0.5, vjust=2))+
  ggtitle("Features with Top Loadings")+
  labs(
    x = '\nMetaCyc Pathway\n',
    y = '\nLoading 1\n',
    fill = 'Health Status')

cf.mds <- cmdscale(d = dist_after)
colnames(cf.mds) <- c("Axis1", "Axis2")
cf.mds.data <- as.data.frame(cf.mds) %>% 
  mutate(score1 = cf.fit_continuous$consensus_scores[, 1])

write.csv(cf.mds.data, 'D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_PHYLUM_results/MMUPHin_DifferentialAbundance/MMUPHin_MetaCycPathways_MDS_ScoringData_UNMAPPED_UNINTEGRATEDremoved_CaseStatus_CaseFollowPairs.csv',
          row.names=FALSE)


ggplot(data=cf.mds.data, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color=score1, fill=score1, shape=shape),
             size=3) +
  scale_shape_manual(values=c(16, 15, 17))+
  guides(fill=guide_colorbar(title = "Score",
                             label.position='right',
                             title.position='top',
                             title.vjust=1,
                             frame.color='black',
                             barwidth=1.5,
                             barheight=15),
         shape=guide_legend(override.aes=list(size=5)),
         color='none',
         size='none')+
  scale_color_viridis(discrete = FALSE, option = 'H')+
  scale_fill_viridis(discrete = FALSE, option = 'H')+
  coord_fixed()+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        legend.key= element_rect(fill=NA),
        legend.key.size = unit(0.75, "cm"),
        legend.spacing.y = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'))+
  labs(fill = "Score",
       shape = "Health Status")
