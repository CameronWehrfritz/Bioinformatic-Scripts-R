#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#October 28, 2020


# CRUK ECM
# Analysis of cancer-adjacent-normal to true-normal 

######################
#### Begin Program ###
######################

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/R_workspace") # VPN mac
setwd("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/R_workspace") # VPN windows
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("ggplot2", "tidyr", "reshape2", "hablar", "VennDiagram", "dplyr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# LOAD DATA #

# Cancer-Adjacent-Normal vs. True Normal:

# First data set (Analysis v1.  E03Nadj, E06Nadj, E14Nadj.  compared to EN01, EN02, EN03)
# df <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Candidates_EN_Eadj_E03E06E14_1003_2020_v1.tsv", sep="\t", stringsAsFactors = F) #VPN mac
df.one <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Candidates_EN_Eadj_E03E06E14_1003_2020_v1.tsv", sep="\t", stringsAsFactors = F) #VPN windows

# Cancer-Adjacent-Normal vs. True Normal
# Second data set (Analysis v2  E03Nadj, E06Nadj, E07Nadj.  compared to EN01, EN02, EN03)
# df.two <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Candidates_EN_Eadj_E03E06E07_1003_2020_v2.tsv", sep="\t", stringsAsFactors = F) #VPN windows
df.two <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Candidates_EN_Eadj_E03E06E07_1003_2020_v2.tsv", sep="\t", stringsAsFactors = F) #VPN windows



# CRUK ECM Panel:

# core panel - all proteins
# df.panel <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Table_Core_Panel_All_1029_2020_v4.csv", sep=",", stringsAsFactors = F) #VPN windows
#version 5 includes 4 more proteins that Birgit wanted to include
# df.panel <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Table_Core_Panel_All_1029_2020_v5.csv", sep=",", stringsAsFactors = F) #VPN windows

# tissue specific
# df.panel <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Table_Core_Panel_Tissue-Specific_1029_2020_v4.csv", sep=",", stringsAsFactors = F) #VPN windows

# universal
# df.panel <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Table_Core_Panel_Universal_1028_2020_v4.csv", sep=",", stringsAsFactors = F) #VPN windows=
#version 5 includes 4 more proteins that Birgit wanted to include
df.panel <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1023_CRUK_ECM_Early_Markers/Input/Table_Core_Panel_Universal_1028_2020_v5.csv", sep=",", stringsAsFactors = F) #VPN windows
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# early marker genes - Birgit wanted to add these in
eary.markers <- c("DSP", "JUP", "PPL", "PKP3")
#------------------------------------------------------------------------------------
  

#--------------------------------------------------------------------------------------------
# Replicates

esoph.replicates <- c("E01", "E03", "E04", "E04", "E06", "E07", "E14")
lung.replicates <- c("L01", "L02", "L03", "L04", "L05", "L06", "L07", "L08", "L10")
gastric.replicates <- c("G01", "G02", "G04", "G05")

all.replicates <- c(esoph.replicates, lung.replicates, gastric.replicates)
#--------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------
# calculate average Qvalues for each tissue
df <- df.panel %>%
  mutate(Esophagus.Avg.Qvalue = rowMeans(select(., starts_with("E") & contains("Qvalue")), na.rm = TRUE)) %>% # calc esophagus avg qvalue
  mutate(Lung.Avg.Qvalue = rowMeans(select(., starts_with("L") & contains("Qvalue")), na.rm = TRUE)) %>% # calc lung avg qvalue
  mutate(Gastric.Avg.Qvalue = rowMeans(select(., starts_with("G") & contains("Qvalue")), na.rm = TRUE)) %>% # calc gastric avg qvalue
  select(UniProtIds, Genes, ProteinNames, ProteinDescriptions, Roles.in.Literature, Tissue.Specific.or.Core, 
         Division, Category, GO.Biological.Process, GO.Cellular.Component, GO.Molecular.Function, Organisms, 
         Esophagus, Lung, Gastric, Esophagus.Avg.Qvalue, Lung.Avg.Qvalue, Gastric.Avg.Qvalue, 
         everything()) # reordering columns
#--------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Join Cancer-Adjacent-Normal vs. True Normal w/ core panel

df <- df %>%
  left_join(select(df.one, E.CAN.vs.TN.1=AVG.Log2.Ratio, E.CAN.vs.TN.1.Qvalue=Qvalue, Genes), by="Genes") %>%
  left_join(select(df.two, E.CAN.vs.TN.2=AVG.Log2.Ratio, E.CAN.vs.TN.2.Qvalue=Qvalue, Genes), by="Genes")
#------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------
# data with only individual biological replicates for each tissue type
df.m.ind <- df %>%
  select(-c(Esophagus, Lung, Gastric, Esophagus.Avg.Qvalue, Lung.Avg.Qvalue, Gastric.Avg.Qvalue)) %>% # drop these columns
  melt()

df.m.ind.qvals <- df.m.ind %>%
  filter(grepl("Qvalue", variable)) %>% # retain only Qvalue rows
  rename(Tissue=variable, Qvalue=value) # rename columns

df.ind <- df.m.ind %>%
  filter(! grepl("Qvalue", variable)) %>% # filter out Qvalue rows
  rename(Tissue=variable, LogRatio=value) %>% # rename columns
  cbind(select(df.m.ind.qvals, Qvalue))


# Heatmap
# Individual Replicates - All
library(scales) # squish and limits are used in controlling color saturation
ggplot(df.ind, aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=20))


# Heatmap
# Individual Replicates - Esophagus only
library(scales) # squish and limits are used in controlling color saturation
df.ind %>%
  filter(grepl("E", Tissue) | Tissue=="E.CAN.VS.TN.1" | Tissue=="E.CAN.VS.TN.2") %>% # keep only esophagus samples
  ggplot(aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=20))
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# Row Clustering and Factor Releveling

# make data frame with only individual replicate columns for clustering
# remove early.marker proteins - these will be added at the bottom
df.ind.tocluster <- df %>%
  filter(! Genes %in% eary.markers) %>%
  select(all_of(esoph.replicates))

row.names(df.ind.tocluster) <- df %>%
  filter(! Genes %in% eary.markers) %>%
  pull(Genes)

# early markers
df.early.markers <- df %>%
  filter(Genes %in% eary.markers)

row.names(df.ind.tocluster) <- df$Genes

hc <- hclust(dist(df.ind.tocluster), "ward.D") # hierarchical clustering

# genes in clustered order
genes.hc.order <- hc$labels[hc$order] ## THIS IS THE ORDER WE WANT THE ROWS IN HEATMAP


# set levels of df before melting, then make plots
df$Genes <- factor(df$Genes, levels = c(eary.markers, genes.hc.order)) # mucho importante
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# melt data - after factor releveling
df.m <- df %>%
  melt()
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# data with only individual biological replicates for each tissue type
df.m.ind <- df %>%
  select(-c(Esophagus, Lung, Gastric, Esophagus.Avg.Qvalue, Lung.Avg.Qvalue, Gastric.Avg.Qvalue)) %>% # drop these columns
  melt()

df.m.ind.qvals <- df.m.ind %>%
  filter(grepl("Qvalue", variable)) %>% # retain only Qvalue rows
  rename(Tissue=variable, Qvalue=value) # rename columns

df.ind <- df.m.ind %>%
  filter(! grepl("Qvalue", variable)) %>% # filter out Qvalue rows
  rename(Tissue=variable, LogRatio=value) %>% # rename columns
  cbind(select(df.m.ind.qvals, Qvalue))


# Heatmap
# Individual Replicates - All
library(scales) # squish and limits are used in controlling color saturation
ggplot(df.ind, aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=16))


# Heatmap
# Individual Replicates - Esophagus only
library(scales) # squish and limits are used in controlling color saturation
df.ind %>%
  filter(grepl("E", Tissue) | Tissue=="E.CAN.VS.TN.1" | Tissue=="E.CAN.VS.TN.2") %>% # keep only esophagus samples
  ggplot(aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=16))

# Heatmap
# Individual Replicates - Lung only
library(scales) # squish and limits are used in controlling color saturation
df.ind %>%
  filter(grepl("L", Tissue)) %>% # keep only esophagus samples
  ggplot(aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=16))
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# Averages across biological replicates
# and with the CAN.vs.TN data

test <- df %>%
  select(UniProtIds, Genes, ProteinNames, ProteinDescriptions, Roles.in.Literature, Tissue.Specific.or.Core, 
         Division, Category, GO.Biological.Process, GO.Cellular.Component, GO.Molecular.Function, Organisms, 
         Esophagus, Lung, Gastric, Esophagus.Avg.Qvalue, Lung.Avg.Qvalue, Gastric.Avg.Qvalue,
         E.CAN.vs.TN.1, E.CAN.vs.TN.1.Qvalue,
         E.CAN.vs.TN.2, E.CAN.vs.TN.2.Qvalue) %>%
  melt() 

# data with only average of biological replicates for each tissue type
df.m.abundances <- test %>%
  filter(! grepl("Qvalue", variable)) %>% # variable must be one of these three strings
  rename(Tissue=variable, LogRatio=value) 


# data with only average qvalues of biological replicates for each tissue type
df.m.qvals <- test %>%
  filter(grepl("Qvalue", variable)) %>% # variable must be one of these three strings
  rename(Tissue.Qvalue=variable, Qvalue=value)

# join abundances and qvalues together
df.m.all <- df.m.abundances %>%
  cbind(select(df.m.qvals, Qvalue)) # add Qvalue column 


# Heatmap
# All
library(scales) # squish and limits are used in controlling color saturation
ggplot(df.m.all, aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=16))


# Heatmap
# Esophagus only
library(scales) # squish and limits are used in controlling color saturation
df.m.all %>%
  filter(grepl("E", Tissue) | Tissue=="E.CAN.VS.TN.1" | Tissue=="E.CAN.VS.TN.2") %>% # keep only esophagus samples
  ggplot(aes(x=Tissue, y=Genes, fill=LogRatio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'white', shape = 21) +
  # geom_vline(xintercept = c(3.5, 4.5), linetype = "dashed") +
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "gray") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) +
  xlab("") + 
  ylab("") + 
  # theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size=16))
#--------------------------------------------------------------------------------------------




# END