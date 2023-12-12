#########################################################
# Metabolome: GNPS Output Organization, Normalization, and Random Forest
#########################################################
### Credit ###
# Code for random forest was written by Douglas V. Guzior and modified by Zoe A. Hansen
# Code for alpha diversity, beta diversity, and heatmap generation was written by Zoe A. Hansen

#########################################################

#### Loading libraries and data ####
library(tidyverse)
library(randomForest) # Had to install version 4.6-14 since the newest version is incompatible with my R version
library(wordspace)
library(data.table)
library(summarytools)
library(reshape2)
library(janitor)
library(ggpubr) 
library(ggfortify)
library(ggpmisc)
library(viridis)
library(FactoMineR)
library(factoextra)
library(MASS)
library(cowplot)
library(vegan)
library(pheatmap)
library(RColorBrewer)

# Case-Follow colors
cf.col <- c('cyan4','darkorchid3')

# Set a function that is opposite of %in% to exclude things from a datatable
`%ni%` <- Negate(`%in%`)

set.seed(808)
dat0 <- read.csv('bucket.csv') #metabolites as columns, "filename" as very first row with individual files as rows
metadata <- read.csv('metadata.csv') #"filenames" as the first row with matches to the output table from MZmine
metadata$HealthStatus <- as.factor(metadata$HealthStatus)
cluster.info <- read.csv('cluster-info.csv') #pulled from GNPS outputs
  
dat1 <- merge(metadata,dat0, by = 'filename')

#Blank Signal Removal & Normalization - 
#Removing all metabolites present in the Blank samples
datB <- data.frame(t(dat0))
blank.names <- dat0[str_detect(dat0$filename, c("Blank")),]; blank.names <- unlist(blank.names$filename)
names(datB) <- as.matrix(datB[1, ]); datB <- datB[-1,]
datB <- setDT(datB, keep.rownames = T); names(datB)[1] <- "Cluster"
datB.name <- datB[,1] %>% mutate(rn = row_number())

datBa <- lapply(datB[,-1], function(x) type.convert(as.numeric(x))) %>% 
  as.data.frame(); datBb <- datBa %>% 
  mutate(rn = row_number())

datBc <- merge(datB.name, datBb, by = "rn"); datBd <- datBc[,-1]; datBd <- datBd %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "Cluster")

datBd$BlankSum <- rowSums(datBd[ ,blank.names])

datBe <- datBd %>% 
  filter(BlankSum < 1) %>% 
  t() %>% 
  as.data.frame()

datBf <- setDT(datBe, keep.rownames = T); names(datBf)[1] <- "filename"
datBf$filename <- gsub("X","",as.character(datBf$filename)); datBf$filename <- gsub("mzML","mzXML",as.character(datBf$filename))

datNB <- slice(datBf, 1:(n() - 1)) #Removes BlankSum row

# Remove unneeded variables from workspace
remove(blank.names,datB,datB.name,datBa,datBb,datBc,datBd,datBe,datBf)

#Normalizing to total signal
raw_pivot0 <- t(datNB)
raw_pivot <- as.data.frame(raw_pivot0)
names(raw_pivot) <- as.matrix(raw_pivot[1, ])
raw_pivot <- raw_pivot[-1,]
raw_pivot[] <- lapply(raw_pivot, function(x) type.convert(as.numeric(x))) #warning: na introduced by coercion is fine
raw_pivot <- raw_pivot %>% drop_na() #solves warning, there must be two files without anything (probably controls)
quant_normT <- scale(raw_pivot, center = F, scale = colSums(raw_pivot)) # sum scaling
quant_norm <- t(quant_normT); quant_norm <- as.data.frame(quant_norm)
quant_norm <- setDT(quant_norm, keep.rownames = TRUE)[]; names(quant_norm)[1] <- "filename"

dat2 <- merge(metadata, quant_norm, by = "filename")
names(dat2) <- make.names(names(dat2)) # Can write this dataframe to CSV for other statistical analysis

# Remove unneeded variables from workspace
remove(raw_pivot0,raw_pivot,quant_norm, quant_normT,datNB)

##########################################################################
############### Performing Statistics on Metabolome ######################

# Use 'dat2' variable and filter out unwanted samples
working <- dat2 %>% filter(SampleType %ni% c("BLANK","POOL")) %>% droplevels()

########## Alpha Diversity ##########
##### Health Status #####
# Richness (of metabolites) across our samples
r.c = specnumber(working[(ncol(metadata)+1):ncol(working)])

# Shannon diversity
h.c = diversity(working[(ncol(metadata)+1):ncol(working)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

# Combine alpha diversity data and Case status information
div.c=tibble(working$PairID, working$HealthStatus, working$Pathogen,r.c, h.c, pielou.c)
colnames(div.c)=c("PairID", "HealthStatus", "Pathogen","Richness", "Shannon", "Pielou")

# Determine the Mean, Min, and Max for these measures
# Health Status
div.means <- div.c %>%
  dplyr::group_by(HealthStatus) %>%
  dplyr::summarise(MeanRich = mean(Richness), MeanShannon = mean(Shannon), MeanPielou = mean(Pielou)) 

div.minmax <- div.c %>%
  dplyr::group_by(HealthStatus) %>%
  dplyr::summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
                   MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))


### Wilcoxon Signed-Rank Test (for Case-Follow Pairs Only) 

wilcox.test(Richness ~ HealthStatus, data = div.c, paired = TRUE)
wilcox.test(Shannon ~ HealthStatus, data = div.c, paired = TRUE)
wilcox.test(Pielou ~ HealthStatus, data = div.c, paired = TRUE)

# Reshape the data to "long" format
div.c.long=melt(div.c, id.vars=c("PairID", "HealthStatus", 'Pathogen'))

cf.compare=list(c('CASE','FOLLOW'))

# Plot alpha diversity - health status #
metab.alpha.div.plot<- ggplot(div.c.long, aes(x=HealthStatus, y=value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=HealthStatus, shape=HealthStatus))+
  facet_wrap(~variable, scales="free") +
  scale_color_manual(values=cf.col) +
  scale_fill_manual(values=cf.col)+
  stat_compare_means(paired = T, label = 'p.format',
                     position='identity',
                     comparisons=cf.compare)+
  theme_bw(base_size = 10) +
  hansen_theme+
  labs(
    x = '\nHealth Status',
    y = 'Alpha Diversity Value\n'
  )

##### Pathogen ######

# Pathogen
div.bact.means <- div.c %>%
  dplyr::group_by(HealthStatus,Pathogen) %>%
  dplyr::summarise(MeanRich=mean(Richness),MeanShannon=mean(Shannon),MeanPielou=mean(Pielou))

div.bact.minmax <- div.c %>%
  dplyr::group_by(HealthStatus,Pathogen) %>%
  dplyr::summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
                   MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))

### Wilcoxon Signed-Rank Test (for Case-Follow Pairs Only) 
div.c.campy <- div.c %>%
  dplyr::filter(Pathogen == 'CAMPYLOBACTER')
div.c.salm<- div.c %>%
  dplyr::filter(Pathogen == 'SALMONELLA')
div.c.shig<- div.c %>%
  dplyr::filter(Pathogen == 'SHIGELLA')
div.c.stec<- div.c %>%
  dplyr::filter(Pathogen == 'STEC')

wilcox.test(Richness ~ HealthStatus, data = div.c.shig, paired = TRUE) # Perform separate Wilcoxon tests for each pathogen
wilcox.test(Shannon ~ HealthStatus, data = div.c.shig, paired = TRUE)
wilcox.test(Pielou ~ HealthStatus, data = div.c.shig, paired = TRUE)


# Isolate Alpha Diversity measures for Different Pathogens (Case, Follow)
div.bact.richness <- div.c.long %>%
  dplyr::filter(variable == 'Richness') 
div.bact.shannon <- div.c.long %>%
  dplyr::filter(variable == 'Shannon') 
div.bact.even <- div.c.long %>%
  dplyr::filter(variable == 'Pielou')

bact.comp <- list(c('SALMONELLA', 'CAMPYLOBACTER'))
colors.bact <- c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2')

# Plot alpha diversity - by pathogen #
metab.even.plot<- ggplot(div.bact.even, aes(x=Pathogen, y=value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Pathogen, shape=Pathogen))+
  facet_wrap(~HealthStatus) +
  scale_color_manual(values=colors.bact) +
  scale_fill_manual(values=colors.bact)+
  stat_compare_means(paired = F, label = 'p.format',
                     position='identity',
                     comparisons=bact.comp)+
  theme_bw(base_size = 10) +
  hansen_theme+
  labs(
    x = '\nPathogen',
    y = 'Evenness\n'
  )

########## Beta Diversity - Health Status ##########
#Calculate Bray-curtis dissimilarity. 
bc.d.c=vegdist(working[(ncol(metadata)+1):ncol(working)], method="bray")
# Perform PCoA
bc.c.pcoa=cmdscale(bc.d.c, eig=TRUE, x.ret=TRUE, k=3) 

# Calculate percent variance explained by each axes 1 and 2
ax1.v.bc.c=bc.c.pcoa$eig[1]/sum(bc.c.pcoa$eig)
ax2.v.bc.c=bc.c.pcoa$eig[2]/sum(bc.c.pcoa$eig)
ax3.v.bc.c=bc.c.pcoa$eig[3]/sum(bc.c.pcoa$eig)

#### Plot beta diversity ####
working_out <- as.data.frame(bc.c.pcoa$points)

metab.beta.div.plot <- ggscatter(working_out, x = 'V1', y = 'V2')+
  geom_point(aes(fill=working$HealthStatus, color=working$HealthStatus, 
                 shape=working$HealthStatus),
             size=4)+  
  hansen_theme+
  theme(legend.position = c(0.02, 0.05), 
        legend.justification = c(0.02, 0.05),
        legend.box.background = element_rect(colour = "black", size=1),
        panel.background = element_rect(colour = "black", size=1, fill=NA))+
  scale_color_manual(values=c('cyan4', 'darkorchid3'))+
  xlab(paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""))+
  ylab(paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""))+
  labs(fill='Health Status',
       color='Health Status',
       shape='Health Status')


### Combine Plots (Figure 1) ###
alpha.beta.comb<- plot_grid(metab.alpha.div.plot, metab.beta.div.plot, 
                            rel_heights = c(1, 0.7),
                            ncol=1, nrow=2,
                            labels = c("A", "B"))

ggsave('D://AlphaBetaDiv_polar_metab_CaseFollow_CombinedPlot_20231206.png',
       plot=alpha.beta.comb,
       width=25, height=38, 
       units='cm',dpi=300, device='png')

### PERMANOVA - by health status
# Centroid
a.c = adonis(bc.d.c~working$HealthStatus, distance=TRUE, permutations=1000)

# Disperson
b.c=betadisper(bc.d.c, group=working$HealthStatus)
permutest(b.c)
plot(b.c)
boxplot(b.c)

# Post-hoc Tukey's Honestly Significant Difference 
TukeyHSD(b.c, which = "group", ordered = FALSE,conf.level = 0.95)

################## Beta Diversity - Pathogen ######################
#Calculate Bray-curtis dissimilarity. 
bc.d.c=vegdist(working[(ncol(metadata)+1):ncol(working)], method="bray")

### Perform PCoA 
bc.c.pcoa=cmdscale(bc.d.c, eig=TRUE, x.ret=TRUE, k=3) 

# Calculate percent variance explained by each axes 1 and 2
ax1.v.bc.c=bc.c.pcoa$eig[1]/sum(bc.c.pcoa$eig)
ax2.v.bc.c=bc.c.pcoa$eig[2]/sum(bc.c.pcoa$eig)
ax3.v.bc.c=bc.c.pcoa$eig[3]/sum(bc.c.pcoa$eig)

# Plot beta diversity by pathogen #
working_out <- as.data.frame(bc.c.pcoa$points)

path.metab.beta.div.plot <- ggscatter(working_out, x = 'V1', y = 'V2')+
  geom_point(aes(fill=working$Pathogen, color=working$Pathogen, 
                 shape=working$HealthStatus),
             size=4)+  
  hansen_theme+
  theme(legend.position = c(0, 0), 
        legend.justification = c(0, 0),
        legend.box.background = element_rect(colour = "black", size=1),
        legend.text = element_text(size=8),
        panel.background = element_rect(colour = "black", size=1, fill=NA))+
  scale_color_manual(values=c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2'))+
  xlab(paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""))+
  ylab(paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""))+
  labs(fill='Pathogen',
       color='Pathogen',
       shape='Health Status')

### Combine plots ###
metab.alpha.comb <- plot_grid(metab.rich.plot + rremove('xlab'), 
                            metab.shan.plot + rremove('xlab'),
                            metab.even.plot + rremove('xlab'),
                            rel_heights = c(0.7,0.7,0.7),
                            ncol=1, nrow=3,
                            labels="C",
                            common.legend=TRUE)

ggsave('D://alpha_div_nonpolar_metab_PATHOGEN_CaseFollow_20231206.png',
       plot=metab.alpha.comb,
       width=25, height=35, 
       units='cm',dpi=300, device='png')

metab.beta.comb <- plot_grid(path.metab.beta.div.plot,
                              labels="D")

ggsave('D://beta_div_nonpolar_metab_PATHOGEN_CaseFollow_20231206.png',
       plot=metab.beta.comb,
       width=25, height=22, 
       units='cm',dpi=300, device='png')

### PERMANOVA - by pathogen
# Centroid
a.c = adonis(bc.d.c~working$Pathogen, distance=TRUE, permutations=1000)

# Disperson
b.c=betadisper(bc.d.c, group=working$Pathogen)

permutest(b.c)

plot(b.c)
boxplot(b.c)

# Post-hoc Tukey's Honestly Significant Difference 
TukeyHSD(b.c, which = "group", ordered = FALSE,conf.level = 0.95)

#############################################################################################
                      
############## Performing Random Forest - Idenfity Distinguishing Metabolites ###############
# Health Status #
rf.healthstatus  <- randomForest(working$HealthStatus ~., data = working[(ncol(metadata)+1):ncol(working)],
                                 na.action = na.omit, importance = T, ntree = 5000)
rf.healthstatus

vip.healthstatus <- varImpPlot(rf.healthstatus, type = 1)
vip <- setDT(data.frame(vip.healthstatus), keep.rownames = T)
names(vip)[1]<-'Cluster'
vip.wIDs <- merge(cluster.info, vip, by= 'Cluster')

t30.healthstatus <- top_n(vip, 30); names(t30.healthstatus)[1] <- 'Cluster'
t30.healthstatus <- merge(cluster.info, t30.healthstatus, by = 'Cluster')

# Pathogen #
working$Pathogen <- as.factor(working$Pathogen)
rf.pathogen  <- randomForest(working$Pathogen ~., data = working[(ncol(metadata)+1):ncol(working)],
                                 na.action = na.omit, importance = T, ntree = 5000)
rf.pathogen

vip.pathogen <- varImpPlot(rf.pathogen, type = 1)
vip.p <- setDT(data.frame(vip.pathogen), keep.rownames = T)
names(vip.p)[1]<-'Cluster'
vip.p.wIDs <- merge(cluster.info, vip.p, by= 'Cluster')

t30.pathogen <- top_n(vip.p, 30); names(t30.pathogen)[1] <- 'Cluster'
t30.pathogen <- merge(cluster.info, t30.pathogen, by = 'Cluster')


# Plot clusters of importance - facet by pathogen #
cluster_plot8369<- ggplot(working, aes(HealthStatus, Cluster8369)) +  # Designate which cluster to plot here #
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=HealthStatus, shape=HealthStatus))+
  scale_color_manual(values=cf.col) +
  scale_fill_manual(values=cf.col)+
  scale_y_continuous(trans = 'pseudo_log', breaks=seq(0, 0.008, by=0.0020)) + 
  stat_compare_means(paired = T, label = 'p.signif',
                     position='identity',
                     label.x=1.5,
                     label.y=0.004) +
  labs(x = '', y = 'Normalized Abundance\n',
       title ='Cluster 8369') + 
  facet_wrap(Pathogen ~., nrow = 1)+
  theme_bw(base_size = 10) +
  hansen_theme

# Plot clusters of importance - facet by health status #
cluster_plot6581<- ggplot(working, aes(Pathogen, Cluster6581)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=Pathogen, shape=Pathogen))+
  scale_color_manual(values=colors.bact) +
  scale_fill_manual(values=colors.bact)+
  scale_y_continuous(trans = 'pseudo_log', breaks=seq(0, 0.006, by=0.0015)) + 
  stat_compare_means(paired = T, label = 'p.signif',
                     position='identity',
                     label.x=1.5,
                     label.y=0.003,
                     size=6) +
  labs(x = '', y = 'Normalized Abundance\n',
       title ='Cluster 6581') + 
  facet_wrap(HealthStatus ~., nrow = 1)+
  theme_bw(base_size = 10) +
  hansen_theme

clusters.comb <- plot_grid(cluster_plot2964,
                           cluster_plot6581,
                           cluster_plot8369,
                           ncol=1, nrow=3,
                           hjust=-0.1,
                           common.legend=TRUE)

ggsave('D://clusters6581_8369_2964_RF_norm_abd_PATHOGEN_health_status_facet_CaseFollow_20231206.png',
       plot=clusters.comb,
       width=26, height=42, 
       units='cm',dpi=300, device='png')

################### Heatmaps - Visualizing Random Forest Output #################

# Create/modify metabolite table
met_table = working[,c(1,4,6,17:ncol(working))]
met_table$filename <- str_extract(met_table$filename, "[A-Z][A-Z]\\d+")

met.meta = met_table[,1:3]
met.cluster = met_table[,c(1,4:ncol(met_table))]

met.cluster <- met.cluster %>%
  gather(key = key, value = value, 2:ncol(met.cluster)) %>%
  spread(key=names(met.cluster)[1], value = 'value') %>%
  dplyr::rename(., filename=key)

# Take the top-30 metabolites identified by Random Forest 
met.cluster.30 <- met.cluster %>%
  filter(filename %in% t30.pathogen$Cluster)

met.cluster.30 <- met.cluster.30 %>% 
  remove_rownames %>% 
  column_to_rownames(var="filename")

met.cluster.30=as.matrix(met.cluster.30)
met.cluster.30 <- met.cluster.30*10^9
mode(met.cluster.30) <- "integer"

# Designate annotation (metadata)
annotation_samples <- met.meta %>% 
  remove_rownames %>% column_to_rownames(var="filename")

anno_color <- list(Pathogen = c(SALMONELLA = "#EE2C2C", SHIGELLA= "#EEB422", 
                               CAMPYLOBACTER="#1C86EE", STEC="#9F79EE"),
                   HealthStatus = c(CASE = "cyan4", FOLLOW = "darkorchid3")
                   )

#HEATMAP
heatmap <- pheatmap(
  mat               = log(met.cluster.30+1),
  border_color      = NA,
  show_colnames     = F,
  show_rownames     = T,
  angle_col = 90,
  #  drop_levels       = TRUE,
  #fontsize_col = 6,
  fontsize_row = 8,
  #  fontsize          = 14,
  # color             = brewer.pal(9,"RdYlBu"),
  #color = inferno(100),  
  #number_color = NA,
  annotation_col = annotation_samples,
  annotation_colors = anno_color,
  annotation_names_col = F,
  annotation_names_row = F,
  cluster_cols = T,
  cluster_rows = T,
  clustering_method = "ward.D2",
  gaps_row = FALSE
)
heatmap

ggsave(plot=heatmap, 
       "heatmap_top30_RF_PATHOGEN_nonpolar_metabolites_20231207.png", 
       width = 12, height = 10)
