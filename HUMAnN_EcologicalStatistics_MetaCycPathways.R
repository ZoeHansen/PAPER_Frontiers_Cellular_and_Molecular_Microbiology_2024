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

meta <- read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/ERIN_MetaboliteAnalysis_Metadata_CaseFollowPairs_clean.csv',
                 header = TRUE) 

comm.level = read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_GENUS_results/ERIN_community_abundance_MetaCycPathways_GENUS_CaseFollow.csv', 
                header = TRUE)

genus.data = read.csv('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_GENUS_results/ERIN_pathabundance_genus_relabd_CaseFollowPairs.csv',
                      header=TRUE)

genus.data <- genus.data %>%
  dplyr::filter(!grepl('Case.status', ER_ID))%>%
  dplyr::filter(!grepl('Case.Follow_ID', ER_ID)) %>%
  dplyr::filter(!grepl('Pathogen', ER_ID))

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

enviro <- meta %>%
  dplyr::select(ER_ID, Case.status, Followup.days, Age.years,Gender, Residence.type, Case.Follow_ID,
                Run, Avg_GS, GenomeEquivalents, Year, Hospital, County, Stool.type) 

enviro$Stool.type[enviro$Stool.type == ''] <- 'Unknown'
enviro$Hospital[enviro$Hospital == ''] <- "Unknown"
enviro$Residence.type[enviro$Residence.type == ''] <- "Unknown"

meta_env <- left_join(health, enviro, by = c('ER_ID','Case.Follow_ID', 'Case.status'))

meta[is.na(meta)] <- 0


# Merge metadata with AGS-normalized abundances
#Full Case_Control Campylobacter dataset with metadata of interest
data_2 <- left_join(IDs_use, data, by = 'ER_ID')


# Remove any lingering zeroes to avoid downstream errors
#data_2 <- data_2 %>%
#  dplyr::select(which(!colSums(data_2[,-1], na.rm=TRUE) %in% 0))
data_2 <- data_2[, colSums(data_2[,c(2:ncol(data_2))] != 0) > 0]
data_2 <- data_2[rowSums(data_2[,c(2:ncol(data_2))] != 0) > 0,] 

# Create comprehensive dataset with features and metadata of interest combined
data_cc <- left_join(health, data_2, by = "ER_ID")

data_cc$ER_ID <- factor(data_cc$ER_ID)
data_cc$Case.status <- factor(data_cc$Case.status)
data_cc$Case.Follow_ID <- factor(data_cc$Case.Follow_ID)
data_cc$Pathogen <- factor(data_cc$Pathogen)
data_cc$Antibiotics <- factor(data_cc$Antibiotics)

data_cc[is.na(data_cc)] <- 0
data_cc <- data_cc[, colSums(data_cc != 0) > 0]
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 

data.met <- data_cc[,c(1:5)]

data.val <- as.data.frame(lapply(data_cc[,c(6:ncol(data_cc))], as.numeric))

data_cc <- cbind(data.met, data.val)

########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Richness (of genes) across our samples
r.c = specnumber(data_cc[,-c(1:5)])

# Shannon diversity
h.c = diversity(data_cc[,-c(1:5)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

### Case Control Follow
# Combine alpha diversity data and Case status information
div.c=tibble(data_cc$ER_ID, data_cc$Case.status, data_cc$Pathogen, data_cc$Case.Follow_ID, r.c, h.c, pielou.c)
colnames(div.c)=c("ER_ID", "Case.status", "Pathogen","Pair","Richness", "Shannon", "Pielou")


# Determine outliers

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
# These samples will be excluded from analysis :( 
# For consistency, we will will also remove their paired samples (ER0235 and ER0535, respectively) from analysis
# This results in a final count of 118 samples (59 pairs)

# COMMUNITY LEVEL
# Even after removal of these four problematic samples, ER0106 appeared as an outlier for the measure
# of richness. However, we will include it in our analysis 

div.c.campy <- div.c %>%
  dplyr::filter(Pathogen == 'Campylobacter (CA)')
div.c.salm<- div.c %>%
  dplyr::filter(Pathogen == 'Salmonella (SA)')
div.c.shig<- div.c %>%
  dplyr::filter(Pathogen == 'Shigella (SH)')
div.c.stec<- div.c %>%
  dplyr::filter(Pathogen == 'STEC (EC)')



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

# Comparisons
case.follow.comp <- list(c('Case','FollowUp'))
colors.cf <- c('cyan4','darkorchid3')

bact.comp <- list(c('Salmonella (SA)', 'Campylobacter (CA)'))
colors.bact <- c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2')

labels <- c(Richness='Richness', Shannon='Shannon Index', Pielou='Pielou Evenness')

# Plot alpha diversity values for Case Status
alpha.div.plot<- ggplot(data=div.c.long, aes(x=Case.status, y=value))+
  geom_boxplot() + 
  #  geom_point(aes(group=Pair, fill=Case.status),size=2,
  #             position = position_dodge(0.2), shape=21)+
  geom_jitter(aes(color=Case.status, shape=Case.status))+
  #  geom_line(aes(group=Pair), position = position_dodge(0.2), color='gray44')+ # For connecting paired data
  facet_wrap(~factor(variable, levels=c('Richness','Shannon','Pielou')),
             scales="free", labeller=labeller(variable=labels)) +
  scale_color_manual(values=colors.cf) +
  scale_fill_manual(values=colors.cf)+
  theme_bw(base_size = 10) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14, face='bold', hjust=0.5),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+
  labs(
    x = '\nHealth Status',
    y = 'Alpha Diversity Value\n'
    #title = 'Metabolic Diversity - MetaCyc Pathways'
  )+
  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)

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

#Calculate Bray-curtis dissimilarity. 
bc.d.c=vegdist(data_cc[,-c(1:5)], method="bray")


#Map desired colors to the Case statuses to make our future legend for the PCoA.  
Class=rep('black' ,nrow(data_cc))
Class[data_cc$Case.status=="Case"]= 'cyan4'
Class[data_cc$Case.status=='FollowUp']= 'darkorchid3'

shape=rep(1, nrow(data_cc))
shape[data_cc$Case.status=='Case']=21
shape[data_cc$Case.status=='FollowUp']=22
shape[data_cc$Antibiotics=='Yes']=24


# Pathogen
bact=rep('black',nrow(data_cc))
bact[data_cc$Pathogen=='Campylobacter (CA)']='dodgerblue2'
bact[data_cc$Pathogen=='Salmonella (SA)']='firebrick2'
bact[data_cc$Pathogen=='Shigella (SH)']='goldenrod2'
bact[data_cc$Pathogen=='STEC (EC)']='mediumpurple2'


######## PCoA #########
### Perform PCoA 
bc.c.pcoa=cmdscale(bc.d.c, eig=TRUE, x.ret=TRUE, k=3) 

# Calculate percent variance explained by each axes 1 and 2
ax1.v.bc.c=bc.c.pcoa$eig[1]/sum(bc.c.pcoa$eig)
ax2.v.bc.c=bc.c.pcoa$eig[2]/sum(bc.c.pcoa$eig)
ax3.v.bc.c=bc.c.pcoa$eig[3]/sum(bc.c.pcoa$eig)

# Add environmental factors if desired
env.sub=meta_env 
envEF.bc=envfit(bc.c.pcoa, env.sub, permutations = 999, na.rm = TRUE)
envEF.bc


### PLOTTING
#### Plot the ordination for Case Status for Axes 1&2 ####
plot(bc.c.pcoa$points[,1],bc.c.pcoa$points[,2],
     cex=1.5,
     pch=shape,
     bg=Class,
     xlab= paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""))#,
     #main = 'Metabolic Beta Diversity - MetaCyc Pathways',
     #xlim=c(-0.6,0.4),
     #ylim=c(-0.6, 0.3))

# Add environmental variables
plot(envEF.bc, p.max=0.05, col='black', cex = 1)#,
#     labels=list(vectors=paste(c('Follow-up Days','Age (years)'))))


### Ellipses
# Case Follow
ordiellipse(bc.c.pcoa, groups=data_cc$Case.status, 
            col=c('cyan4','darkorchid3'),
            draw = 'lines', kind = 'sd', conf = 0.95, alpha = 0.05, label = TRUE, lwd = 2)

# Bacterial Pathogen 
ordiellipse(bc.c.pcoa, groups=data_cc$Pathogen, 
            col=c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2'),
            draw = 'lines', kind = 'sd', conf = 0.95, alpha = 0.05, label = TRUE, lwd = 2)


### Legend
# Case Follow
legend('topright', c('Case','FollowUp','Antibiotics'), pch = c(21,22,24), 
       pt.bg=c('cyan4', 'darkorchid3','black'), lty=0)

# Pathogen
legend('bottomright',c('Campylobacter (CA)','Salmonella (SA)','Shigella (SH)','STEC (EC)','Case','Follow-Up','Antibiotics'), 
       pch=c(21,21,21,21,21,22,24), pt.bg=c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2','black','black','black'))

### Text Labels
#Add ID labels to points (if desired)
textxy(X=bc.c.pcoa$points[,1], Y=bc.c.pcoa$points[,2],labs=data_cc$Reassigned.Pair.ID, 
       cex=0.6, pos=2)

#### Beta Div with GGPLOT ####

data_cc$Abx_key <- ifelse(data_cc$Antibiotics == 'Yes', '*','')

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
  #geom_text(label=data_cc$Reassigned.Pair.ID, hjust=1.5, vjust=0.2)+
  #  stat_ellipse(aes(group=data_cc$Case.status))+
  #geom_line(aes(group=data_cc$Reassigned.Pair.ID), color='gray60')+
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

ggsave('D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_GENUS_results/EcologicalStatistics/CommunityLevel/AlphaBetaDiv_MetaCycPathways_CommLevel_CaseFollow_CombinedPlot_20231204.png',
       plot=alpha.beta.comb,
       width=25, height=38, 
       units='cm',dpi=300, device='png')



#### Plot with ggplot2() to connect lines

data_cc_out <- as.data.frame(bc.c.pcoa$points)

ggscatter(data_cc_out, x = 'V1', y = 'V2', shape = 21)+
  geom_point(fill=Class, size=4, color=Class)+
  geom_text(label=data_cc$Case.Follow_ID, hjust=1.8, vjust=0.2)+
  #  stat_ellipse(aes(group=data_cc$Case.status))+
  geom_line(aes(group=data_cc$Case.Follow_ID), color='gray60')+
  xlab(paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""))+
  ylab(paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""))



##### Plotting Axes 2&3 #####
plot(bc.c.pcoa$points[,2],bc.c.pcoa$points[,3],
     cex=1.5,
     pch=shape,
     bg=Class,
     xlab= paste("PCoA2: ",100*round(ax2.v.bc.c,3),"% var. explained",sep=""), 
     ylab= paste("PCoA3: ",100*round(ax3.v.bc.c,3),"%var. explained",sep=""))#,

plot(envEF.bc, p.max=0.05, col='black', cex = 1)

textxy(X=bc.c.pcoa$points[,2], Y=bc.c.pcoa$points[,3],labs=data_cc$Case.Follow_ID, 
       cex=0.6, pos=2)

# Case Follow
legend('bottomleft', c('Case','FollowUp','Antibiotics'), pch = c(21,22,24), 
       pt.bg=c('cyan4', 'darkorchid3','black'), lty=0)


data_cc_out23 <- as.data.frame(bc.c.pcoa$points)
ggscatter(data_cc_out23, x = 'V2', y = 'V3', shape = 21)+
  geom_point(fill=Class, size=4, color=Class)+
  geom_text(label=data_cc$Case.Follow_ID, hjust=1.5, vjust=0.2)+
  geom_line(aes(group=data_cc$Case.Follow_ID), color='gray60')+
  xlab(paste("PCoA2: ",100*round(ax2.v.bc.c,3),"% var. explained",sep=""))+
  ylab(paste("PCoA3: ",100*round(ax3.v.bc.c,3),"%var. explained",sep=""))

##### PCoA EnvFit - Investigate species or genus ####

# Intrinsic variables (species or genus)
pcoa.mod.fit <- envfit(bc.c.pcoa, data_cc[,c(6:ncol(data_cc))], permutations=999)
head(pcoa.mod.fit)

pcoa.mod.df <- data.frame((pcoa.mod.fit$vectors)$arrows, (pcoa.mod.fit$vectors)$r,
                          (pcoa.mod.fit$vectors)$pvals)


write.csv(pcoa.mod.df, 'D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/HUMAnN/humann3_GENUS_results/EcologicalStatistics/GenusLevel/ERIN_PCoA_EnvFit_intrinsic_VECTORS_MetaCyc_GenusLevel_CaseFollowPairs.csv',
          row.names=TRUE)

# Plot Intrinsic Variables

# First, subset based on R2 value (otherwise, way too many arrows)
#__FUNCTION: select.envfit__# (found from https://www.researchgate.net/post/How_do_I_set_cutoff_r_values_for_plotting_arrows_from_function_envfit_in_R)
# function (select.envit) filters the resulting list of function (envfit) based on their p values. This allows to display only significant values in the final plot.

select.envfit<-function(fit, r.select){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(fit$vectors$r)) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    if (fit$vectors$r[i]<r.select) { #Check wether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    } #close if-loop
  } #close for-loop
  return(fit) #return fit as the result of the function
} #close the function

pcoa.mod.fit2 <-select.envfit(pcoa.mod.fit, r.select=0.85) # Perform select.envfit on dataset


plot(bc.c.pcoa$points[,1],bc.c.pcoa$points[,2],
     cex=1.5,
     pch=shape,
     bg=Class,
     xlab= paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""),
     main = 'Metabolic Beta Diversity - MetaCyc Pathways',
     xlim=c(-0.9,0.3),
     ylim=c(-0.5,0.3))


legend('topleft', c('Case','FollowUp','Antibiotics'), pch = c(21,22,24), 
       pt.bg=c('cyan4', 'darkorchid3','black'), lty=0)

plot(pcoa.mod.fit2, p.max=0.001, col = "black", cex = 0.7, arrow.mul=0.4)


# Extrinsic variables (Environmental factors)
pcoa.env.fit <- envfit(bc.c.pcoa, env.sub, permutations = 999, na.rm = TRUE)
pcoa.vect.df <- data.frame((pcoa.env.fit$vectors)$arrows, (pcoa.env.fit$vectors)$r,
                           (pcoa.env.fit$vectors)$pvals)
pcoa.fact.df <- data.frame((pcoa.env.fit$factors)$r,
                           (pcoa.env.fit$factors)$pvals)

write.csv(pcoa.fact.df,'D://Manning_ERIN/ERIN_Metabolomics_AIM_THREE/Anvio_Pipeline/Frequency/EcologicalStatistics/ERIN_PCoa_EnvFit_extrinsic_FACTORS_KEGGmodules_CaseFollowPairs.csv',
          row.names = TRUE)

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

