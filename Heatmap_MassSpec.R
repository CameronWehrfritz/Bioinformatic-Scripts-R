#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute
#March 22, 2021

# Create Heatmap of Mass Spec Data

#### Begin Program ###

#--------------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2020_0701_CRUK_STORMing_Cancer_MASTER/2021_0310_E20_CRUK/R_workspace") # mac
setwd("//bigrock/GibsonLab/users/Cameron/2020_0701_CRUK_STORMing_Cancer_MASTER/2021_0310_E20_CRUK/R_workspace") # pc
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# PACKAGES 
packages = c("dplyr", "openxlsx", "readxl", "writexl", "reshape2", "ggplot2", "VennDiagram",
             "scales", "forcats", "tidyr")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#--------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------
# LOAD DATA

# CRUK .csv files are stored in the repository MASTER_Processed_Input_CSV_Files
# these .csv files have been processed from the original raw .tsv files

df.input <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_0701_CRUK_STORMing_Cancer_MASTER/MASTER_Processed_Input_CSV_Files/Esophagus/Candidates_PG20_E20_0305_2021_v2.csv", sep=",", encoding = "UTF-8-BOM")
names(df.input)[1] <- "UniProtIds" # fix column name
#--------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# CLEAN THE DATA

# check out the comparisons
unique(df.input$Comparison..group1.group2.)

df <- df.input %>%
  filter(grepl("E20_bT", Comparison..group1.group2.)) %>% # keep only ratios containing Tumor
  mutate(Comparison..group1.group2. = gsub(" / ", "_vs_", Comparison..group1.group2.)) %>% # tweak comparison name
  arrange(Qvalue) %>% # order by ascending Qvalue
  mutate(Uni_Gene = paste(UniProtIds, Genes, sep="_")) # create new variable incorporating uniprotid and gene name together

# check out the comparisons - these are the ones which will be plotted in the heatmap
unique(df$Comparison..group1.group2.)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# FILTER FOR SIGNIFICANCE
# by default we are using Q<0.05 and absolute.log2.ratio>0.58 

df <- df %>%
  mutate(AVG.Log2.Ratio=ifelse(Qvalue<0.05 & Absolute.AVG.Log2.Ratio>0.58, AVG.Log2.Ratio, NA)) # replace non-significant observations with NA
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# MELT

df.m <- df %>%
  melt( measure.vars=c("AVG.Log2.Ratio", "Qvalue")) # melt so each row is one of these observations
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# CAST WIDE
# cast avg.log2.ratio and qvalue wide

df.wide <- df.m %>%
  # first cast avg.log2.ratio wide
  filter( grepl("AVG.Log2.Ratio", variable)) %>%
  dcast(formula = Uni_Gene + UniProtIds + Genes + Comparison..group1.group2. ~ variable, value.var=c("value")) %>%
  # join with Qvalues
  inner_join(df.m %>% 
               filter( grepl("Qvalue", variable)) %>%
               dcast(formula = Uni_Gene + UniProtIds + Genes + Comparison..group1.group2. ~ variable, value.var=c("value")), 
             by = c("Uni_Gene", "UniProtIds", "Genes", "Comparison..group1.group2.")) # join by these variables
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# CREATE HEATMAP DATAFRAME

# determine proteins to show in heatmap
df.heatmap <- df.wide %>%
  # determine how many proteins to display in heatmap
  group_by(UniProtIds, Genes) %>% # group by UniProtId and Gene
  mutate(Number.Sig.Obs=sum(!is.na(AVG.Log2.Ratio))) %>% # calculate number of significant observations per uniprot/gene
  ungroup %>%
  filter(Number.Sig.Obs==4) # adjust this number to see more or less proteins in the heatmap
  
# find protein order 
# set up for doing hierarchical clustering
df.to.cluster <- df.m %>%
  filter(Genes %in% df.heatmap$Genes) %>% # only keep proteins in df.heatmap
  filter(variable=="AVG.Log2.Ratio") %>% # only keep Avg.Log2.Ratio for clustering
  dcast(formula = Uni_Gene + UniProtIds + Genes ~ Comparison..group1.group2., value.var = "value")

# hierarchical clustering
hierarchical.cluster <- df.to.cluster %>% 
  select(-UniProtIds, -Genes, -Uni_Gene) %>%
  dist() %>%
  hclust(method = "ward.D")

# protein order after hierarchical cluster 
protein.hc.order <- df.to.cluster[hierarchical.cluster$order, ] %>%
  pull(Uni_Gene)

# rearrange heatmap according to hierarchical clustered ordering
df.heatmap <- df.heatmap %>% 
  mutate(Uni_Gene = factor(Uni_Gene, levels = protein.hc.order)) %>% # relevel factor
  arrange(Uni_Gene) # rearrange rows
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# PLOT HEATMAP

df.heatmap %>%
  ggplot(aes(x=Comparison..group1.group2., y=Uni_Gene, fill=AVG.Log2.Ratio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'black', shape = 21) + # depict Qvalue with a black circle
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  xlab("") +
  ylab("") +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, size=8), axis.text.y = element_text(size=7))
#-----------------------------------------------------------------------------------------------------

# END
