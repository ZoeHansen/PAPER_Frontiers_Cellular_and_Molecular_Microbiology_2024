########################################################
# Aim 3: HUMAnN 3.0 - Statistics for Case/Follow Pairs
########################################################

# Load libraries and data

library(tidyverse)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(vegan)
library(calibrate)

meta <- read.csv('D://ERIN_MetaboliteAnalysis_Metadata_CaseFollowPairs_clean.csv', header = TRUE) 
comm.level = read.csv('D://ERIN_community_abundance_MetaCycPathways_GENUS_CaseFollow.csv', header = TRUE) # Relative abundances from HUMAnN 3.0

### Sample names should be in rows, features in columns (transpose if needed)
data <- comm.level %>%
  gather(key = key, value = value, 2:ncol(comm.level)) %>%
  spread(key=names(comm.level)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

data <- data[, colSums(data[, -1] != 0) > 0] 
data <- data[rowSums(data[,-1] != 0) > 0,] 
data[is.na(data)] <- 0

########################################################
# Selecting relevant metadata for analysis
########################################################

# Select the relevant metadata fields to include in analysis
health <- meta %>%
  dplyr::select(ER_ID, Case.status, Pathogen, Antibiotics, Case.Follow_ID) 

health$Antibiotics[health$Antibiotics == ""] <- "No"

IDs_use <- health %>%
  dplyr::select(ER_ID)

# Merge metadata with abundances
data_2 <- left_join(IDs_use, data, by = 'ER_ID')

# Remove any lingering zeroes to avoid downstream errors
data_2 <- data_2[, colSums(data_2[,c(2:ncol(data_2))] != 0) > 0]
data_2 <- data_2[rowSums(data_2[,c(2:ncol(data_2))] != 0) > 0,] 

# Create comprehensive dataset with features and metadata of interest combined
data_cc <- left_join(health, data_2, by = "ER_ID")

# Change variables of interest to factors (where appropriate)
data_cc$ER_ID <- factor(data_cc$ER_ID)
data_cc$Case.status <- factor(data_cc$Case.status)
data_cc$Case.Follow_ID <- factor(data_cc$Case.Follow_ID)
data_cc$Pathogen <- factor(data_cc$Pathogen)
data_cc$Antibiotics <- factor(data_cc$Antibiotics)

data_cc[is.na(data_cc)] <- 0
data_cc <- data_cc[, colSums(data_cc != 0) > 0]
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 

########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Richness (of genes) across our samples
r.c = specnumber(data_cc[,-c(1:5)])

# Shannon diversity
h.c = diversity(data_cc[,-c(1:5)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

### Case-FollowUp Pairs
# Combine alpha diversity data and Case status information
div.c=tibble(data_cc$ER_ID, data_cc$Case.status, data_cc$Pathogen, data_cc$Case.Follow_ID, r.c, h.c, pielou.c)
colnames(div.c)=c("ER_ID", "Case.status", "Pathogen","Pair","Richness", "Shannon", "Pielou")

# Determine outliers #
out.rich <- boxplot.stats(div.c$Richness)$out
out_ind.rich <- which(div.c$Richness %in% c(out.rich))
outliers.rich <- div.c[out_ind.rich, ]

out.shan <- boxplot.stats(div.c$Shannon)$out
out_ind.shan <- which(div.c$Shannon %in% c(out.shan))
outliers.shan <- div.c[out_ind.shan, ]

out.even <- boxplot.stats(div.c$Pielou)$out
out_ind.even <- which(div.c$Pielou %in% c(out.even))
outliers.even <- div.c[out_ind.even, ]

# ER0194 (Case) and ER0611 (FollowUp) were found to be outliers in each of the three measures of Alpha Diversity
# These samples were excluded from analysis
# For consistency, we will will also remove their paired samples (ER0235 and ER0535, respectively) from analysis
# This results in a final count of 118 samples (59 pairs)

### COMMUNITY LEVEL - Abundances #
# Even after removal of these four problematic samples, ER0106 appeared as an outlier for the measure
# of richness. However, we will include it in our analysis 

# Determine the Mean, Min, and Max for these measures
# Health Status
div.means <- div.c %>%
  dplyr::group_by(Case.status) %>%
  dplyr::summarise(MeanRich = mean(Richness), MeanShannon = mean(Shannon), MeanPielou = mean(Pielou)) 
div.minmax <- div.c %>%
  dplyr::group_by(Case.status) %>%
  dplyr::summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
                   MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))

# Pathogen
div.bact.means <- div.c %>%
  dplyr::group_by(Case.status,Pathogen) %>%
  dplyr::summarise(MeanRich=mean(Richness),MeanShannon=mean(Shannon),MeanPielou=mean(Pielou))
div.bact.minmax <- div.c %>%
  dplyr::group_by(Case.status,Pathogen) %>%
  dplyr::summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
                   MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))


### Plot alpha diversity ###
hansen_theme = theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14, face='bold', hjust=0.5),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+

# Reshape the data to "long" format
div.c.long=melt(div.c, id.vars=c("ER_ID", "Case.status", 'Pathogen','Pair'))
div.c.long$variable_f = factor(div.c.long$variable, levels=c('Richness','Shannon','Pielou'))

# Isolate Alpha Diversity measures for Different Pathogens (Case, Control, Follow)
div.bact.richness <- div.c.long %>%
  dplyr::filter(variable == 'Richness') 
div.bact.shannon <- div.c.long %>%
  dplyr::filter(variable == 'Shannon') 
div.bact.even <- div.c.long %>%
  dplyr::filter(variable == 'Pielou')

# Comparisons (for stat_compare_means)
case.follow.comp <- list(c('Case','FollowUp'))
colors.cf <- c('cyan4','darkorchid3')

bact.comp <- list(c('Salmonella (SA)', 'Campylobacter (CA)'))
colors.bact <- c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2')

labels <- c(Richness='Richness', Shannon='Shannon Index', Pielou='Pielou Evenness')

# Plot alpha diversity values for Health Status
alpha.div.plot<- ggplot(data=div.c.long, aes(x=Case.status, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(color=Case.status, shape=Case.status))+
  facet_wrap(~factor(variable, levels=c('Richness','Shannon','Pielou')),
             scales="free", labeller=labeller(variable=labels)) +
  scale_color_manual(values=colors.cf) +
  scale_fill_manual(values=colors.cf)+
  theme_bw(base_size = 10) +
  hansen_theme+
  labs(
    x = '\nHealth Status',
    y = 'Alpha Diversity Value\n'
  )+
  stat_pvalue_manual(means.to.plot, label = "p.adj", tip.length = 0.01)

### Wilcoxon Signed-Rank Test (for Case-Follow Pairs Only) 
wilcox.test(Richness ~ Case.status, data = div.c, paired = TRUE)
wilcox.test(Shannon ~ Case.status, data = div.c, paired = TRUE)
wilcox.test(Pielou ~ Case.status, data = div.c, paired = TRUE)

# For adding to plot
means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",    1.212e-07,    325,       "Richness",
  "Case",    "FollowUp",    7.49e-10,    3.5,         "Shannon",
  "Case",    "FollowUp",    1.67e-09,    0.63,       "Pielou"
)

########################################################
# Calculate and plot across-sample (Beta) diversity
########################################################

#Calculate Bray-curtis dissimilarity of HUMAnN pathways  
bc.d.c=vegdist(data_cc[,-c(1:5)], method="bray")

### Perform PCoA 
bc.c.pcoa=cmdscale(bc.d.c, eig=TRUE, x.ret=TRUE, k=3) 

# Calculate percent variance explained by each axes 1 and 2
ax1.v.bc.c=bc.c.pcoa$eig[1]/sum(bc.c.pcoa$eig)
ax2.v.bc.c=bc.c.pcoa$eig[2]/sum(bc.c.pcoa$eig)
ax3.v.bc.c=bc.c.pcoa$eig[3]/sum(bc.c.pcoa$eig)

#### Plot Beta Div ####
data_cc_out <- as.data.frame(bc.c.pcoa$points)

beta.div.plot <- ggscatter(data_cc_out, x = 'V1', y = 'V2')+
  geom_point(aes(fill=data_cc$Case.status, color=data_cc$Case.status, 
                 shape=data_cc$Case.status),
             size=4)+  
  geom_text(aes(label=data_cc$Abx_key), size=10)+
  theme(legend.position = c(0.02, 0.99), 
        legend.justification = c(0.02, 0.99),
        legend.box.background = element_rect(colour = "black", size=1),
        panel.background = element_rect(colour = "black", size=1, fill=NA))+
  scale_color_manual(values=c('cyan4', 'darkorchid3'))+
  xlab(paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""))+
  ylab(paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""))+
  labs(fill='Health Status',
       color='Health Status',
       shape='Health Status')

############## COMBINE PLOTS (Figure 1) #####################

alpha.beta.comb<- plot_grid(alpha.div.plot, beta.div.plot, 
                            rel_heights = c(1, 0.8),
                            ncol=1, nrow=2,
                            labels = c("A", "B"))

ggsave('D://AlphaBetaDiv_MetaCycPathways_CommLevel_CaseFollow_CombinedPlot_20231204.png',
       plot=alpha.beta.comb,
       width=25, height=38, 
       units='cm',dpi=300, device='png')

########################################################
# PERMANOVA
########################################################
# Centroid
a.c = adonis(bc.d.c~data_cc$Case.status, distance=TRUE, permutations=1000)

# Disperson
b.c=betadisper(bc.d.c, group=data_cc$Case.status)
permutest(b.c)
plot(b.c)
boxplot(b.c)

# Post-hoc Tukey's Honestly Significant Difference 
TukeyHSD(b.c, which = "group", ordered = FALSE,conf.level = 0.95)

