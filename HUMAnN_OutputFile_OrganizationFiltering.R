#########################################################

# HUMAnN 3.0 Output Organization and Filtering 

#########################################################

# Author : Zoe Hansen
# Last modified : 2022.02.09

# This code will use the 'pathcoverage.tsv' file to determine the average pathway coverage among samples
# Samples will then be filtered out of the 'pathabundance.tsv' file based on these filtering cutoffs. 

# A subset will be made for Case-Follow pairs included in the other metabolome analyses. The entire 
# sample-set will still be stored but not included in downstream analyses. 

############## Loading Libraries and Data #########################

library(tidyverse)
library(plyr)

path.abd <- read.delim('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/ERIN_humann_pathabundance_relabd_merged.tsv',
                       header=TRUE, sep='\t')

path.cov <- read.delim('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_GENUS_results/ERIN_pathcoverage_genus_merged.tsv',
                       header=TRUE, sep='\t')

metadata <- read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/ERIN_MetaboliteAnalysis_Metadata_CaseFollowPairs_clean.csv',
                     header=TRUE)

########### Isolate Case/Follow Pairs ######################

meta_sub <- metadata %>%
  dplyr::select(ER_ID, Case.status, Case.Follow_ID, Pathogen)

path.abd.t <- path.abd %>%
  gather(key = key, value = value, 2:ncol(path.abd)) %>%
  spread(key=names(path.abd)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., ER_ID=key)

path.abd.cf <- left_join(meta_sub, path.abd.t, by='ER_ID')

path.abd.cf <- path.abd.cf[, colSums(path.abd.cf[,-c(1:4)] != 0) > 0] 
path.abd.cf <- path.abd.cf[rowSums(path.abd.cf[,-c(1:4)] != 0) > 0,] 

path.abd.cf <-left_join(meta_sub, path.abd.cf, by=c('ER_ID','Pathogen'))


path.cov.t <- path.cov %>%
  gather(key = key, value = value, 2:ncol(path.cov)) %>%
  spread(key=names(path.cov)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., ER_ID=key)

path.cov.cf <- left_join(meta_sub, path.cov.t, by='ER_ID')

path.cov.cf <- path.cov.cf %>%
  dplyr::select(which(!colSums(path.cov.cf[,-c(1:4)], na.rm=TRUE) %in% 0))
path.cov.cf <- path.cov.cf[rowSums(path.cov.cf[,-c(1:4)] !=0) >0, ]


########### Determining Average Pathway Coverage ###########

# Since we have varied pathway coverage among our samples (due to variability that we may not want to 
# overlook...) I am going to hold off on filtering based on pathway coverage. However, here is code
# designed to perform filtering based on the average coverage estimate per pathway across all 118 samples.

# According to Eric Franzosa who helped develop HUMAnN 3.0:

#   To clarify this: We are leaning toward replacing the pathway coverage output file 
#   (which is a holdover from the original HUMAnN) with an intermediate abundance file, 
#   probably reaction abundances. Modern HUMAnN is tuned such that a non-zero pathway abundance 
#   is a good indicator of our confidence in the pathway's presence, independent of the separate 
#   probability estimate provided by the coverage file.

path.cov.cf <- path.cov.cf %>%
  dplyr::select(-Case.status, -Case.Follow_ID, -Pathogen)

path.cov.cf.t <- path.cov.cf %>%
  gather(key = key, value = value, 2:ncol(path.cov.cf)) %>%
  spread(key=names(path.cov.cf)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., ER_ID=key)

path.cov.avg <- path.cov.cf.t %>%
  mutate(., AvgCov = rowSums(path.cov.cf.t[,-1])/ncol(path.cov.cf.t[,-1]))

path.cov.cutoff <- path.cov.avg %>%
  filter(AvgCov > 0.1)


########### Determining Average Pathway Abundance ##################

# Instead of filtering based on coverage, I will compute the average relative abundances of these pathways
# and exclude all pathways that register an average equal to 0. This will remove all irrelevant/absent 
# pathways from our dataset. 

path.abd.cf <- path.abd.cf %>%
  dplyr::select(-Case.status, -Case.Follow_ID, -Pathogen)

path.abd.cf.t <- path.abd.cf %>%
  gather(key = key, value = value, 2:ncol(path.abd.cf)) %>%
  spread(key=names(path.abd.cf)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., ER_ID=key)

path.abd.avg <- path.abd.cf.t %>%
  mutate(., AvgAbd = rowSums(path.abd.cf.t[,-1])/ncol(path.abd.cf.t[,-1]))


# Isolate pathways with average relative abundance > 0 (excluding non-relevant pathways)
path.abd.cutoff <- path.abd.avg %>%
  filter(., AvgAbd > 0)

write.csv(path.abd.cutoff, 'D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/ERIN_pathabundance_relabd_merged_CaseFollowPairs.csv',
          row.names=FALSE)

