#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute
#March 22, 2021

# Create Heatmap of Mass Spec Data

#### Begin Program ###

#-----------------------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2020_0701_CRUK_STORMing_Cancer_MASTER/2021_0310_E20_CRUK/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2020_0701_CRUK_STORMing_Cancer_MASTER/2021_0310_E20_CRUK/R_workspace") # PC
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# PACKAGES 
packages = c("dplyr", "openxlsx", "readxl", "writexl", "reshape2", "ggplot2", "VennDiagram",
             "scales", "forcats", "tidyr")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# LOAD DATA

# note about CRUK STORMing Cancer Data:
# CRUK .csv files are stored in the repository MASTER_Processed_Input_CSV_Files
# these .csv files have been processed from the original raw .tsv files

df.input <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_0701_CRUK_STORMing_Cancer_MASTER/MASTER_Processed_Input_CSV_Files/Esophagus/Candidates_PG20_E20_0305_2021_v2.csv", sep=",", encoding = "UTF-8-BOM")
names(df.input)[1] <- "UniProtIds" # fix column name
#-----------------------------------------------------------------------------------------------------


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
# typically we are using Q<0.05 and absolute.log2.ratio>0.58 
# slightly more stringent is Q<0.01 and absolute.log2.ratio>0.58 

# set cutoff values for filtering
Qvalue.cutoff <- 0.01
Log.Ratio.cutoff <- 0.58

df <- df %>%
  mutate(AVG.Log2.Ratio=ifelse(Qvalue<Qvalue.cutoff & Absolute.AVG.Log2.Ratio>Log.Ratio.cutoff, AVG.Log2.Ratio, NA)) # replace non-significant observations with NA
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

# rearrange heatmap rows according to hierarchical clustered ordering
df.heatmap <- df.heatmap %>% 
  mutate(Uni_Gene = factor(Uni_Gene, levels = protein.hc.order)) %>% # relevel factor
  arrange(Uni_Gene) # rearrange rows

# # do you want to rearrange the columns (which show each comparison)? if not, skip this step
# df.heatmap <- df.heatmap %>%
#   mutate(Comparison..group1.group2. =  factor(Comparison..group1.group2., levels = c("ENTER COMPARISON HERE", "ENTER COMPARISON HERE", "ETC")))
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# PLOT HEATMAP WITH GGPLOT

# create title
title <- paste("Significance criteria: \n Qvalue < ", Qvalue.cutoff, sep="") %>%
  paste(" \n |Log2.Ratio| > ", Log.Ratio.cutoff, sep="")

# for shorter protein lists
df.heatmap %>%
  ggplot(aes(x=Comparison..group1.group2., y=Uni_Gene, fill=AVG.Log2.Ratio)) +
  coord_equal() +
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'black', shape = 21) + # depict Qvalue with a black circle
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  xlab("Comparisons") +
  ylab("Proteins") +
  ggtitle(title) +
  theme(title = element_text(size=10), text = element_text(size = 20), axis.text.x = element_text(angle = 90, size=8), axis.text.y = element_text(size=7))

# for very long protein list, very small protein font and wide aspect ratio
df.heatmap %>%
  ggplot(aes(x=Comparison..group1.group2., y=Uni_Gene, fill=AVG.Log2.Ratio)) +
  coord_equal(ratio = 1/3) + # aspect ratio y/x
  geom_tile(color = "white") +
  geom_point(aes(size = -log10(Qvalue)), colour = 'black', shape = 21) + # depict Qvalue with a black circle
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-3,3), oob=squish) + # squish and limits are used in controlling color saturation # library(scales)
  xlab("Comparisons") +
  ylab("Proteins") +
  ggtitle(title) +
  theme(title = element_text(size=10), text = element_text(size = 20), axis.text.x = element_text(angle = 90, size=8), axis.text.y = element_text(size=2))
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# PLOT HEATMAP BASIC (not GGPLOT)
# better for very long protein lists

# create title
title <- paste("Qvalue < ", Qvalue.cutoff, sep="") %>%
  paste(" \n |Log2.Ratio| > ", Log.Ratio.cutoff, sep="")

# function for calculating font size for heatmap row labels
heatmap.font.size.function <- function(x){
  # input x is the number of proteins to be displayed in the heatmap
  if(all(x>0, x<=135)){
    # use this linear function to calculate font size
    font.size <- -0.007*x + 1.07
  }
  if(x>135){
    # default to a very small font size for large lists of proteins
    font.size <- 0.1 
  }
  if(x<0){
    # x must be positive
    font.size <- NA
  }
  return(font.size)
}

# wide format heatmap for non-ggplot heatmap
df.basic <- df.heatmap %>%
  dcast(formula = Uni_Gene + UniProtIds + Genes ~ Comparison..group1.group2., value.var=c("AVG.Log2.Ratio"))

# keep only data 
df.heatmap.input.temp <- df.basic %>%
  select(-Uni_Gene, -UniProtIds, -Genes)

# make data frame numeric
df.heatmap.input <- apply(df.heatmap.input.temp, 2, as.numeric)

# assign rownames for heatmap 
rownames(df.heatmap.input) <- df.basic$Uni_Gene

# heatmap layout options to customize the size of the dendrograms, key, etc...
lmat = rbind(c(4,3),c(2,1))
lwid = c(1,4)
lhei = c(1,4)
layout(mat = lmat, widths = lwid, heights = lhei)


graphics.off()
pdf(file="Heatmap_test.pdf", 8.25, 11)

# Heatmap
heatmap.2(df.heatmap.input, 
          main=paste("TITLE \n", title, sep=" "), # heatmap title
          scale="none", # scale to row or column or none
          col=bluered, # color key palette
          trace="none",
          breaks = seq(-4,4,2*0.001), # specifies when to saturate with color, and how many steps
          key=TRUE, # includes a color key
          symkey=FALSE,
          density.info="none", 
          margins = c(10, 10), 
          Rowv=TRUE, # row clustering   
          Colv=FALSE, # column clustering   
          na.color = "gray", # color of NA values
          dendrogram="none", 
          cexCol=0.75, # column text size
          cexRow=heatmap.font.size.function(nrow(df.heatmap.basic)), # row text size
          key.title="", # title of the key
          key.xlab="", # key x axis labels
          lmat=lmat, # requires heatmap layout above to be specified
          lwid=lwid, # requires heatmap layout above to be specified
          lhei=lhei, # requires heatmap layout above to be specified 
          key.par = list(cex=1.05, mar=c(5,3,3.5,0)), # mar is A numerical vector of the form c(bottom, left, top, right) 
          #which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
          keysize = 1.5
)
# #add timestamp
# mtext(date(), side=1, line=4, adj=0) # side (1=bottom, 2=left, 3=top, 4=right)
graphics.off()

# # optional
# # same heatmap, with column clustering
# graphics.off()
# pdf(file="Heatmap_colclustered_test.pdf", 8.25, 11)
# 
# #Heatmap
# heatmap.2(df.heatmap.input, 
#           main=paste("TITLE \n", title, sep=" "), # heatmap title
#           scale="none", # scale to row or column or none
#           col=bluered, # color key palette
#           trace="none",
#           breaks = seq(-4,4,2*0.001), # specifies when to saturate with color, and how many steps
#           key=TRUE, # includes a color key
#           symkey=FALSE,
#           density.info="none", 
#           margins = c(10, 10), 
#           Rowv=TRUE, # row clustering   
#           Colv=TRUE, # column clustering   
#           na.color = "gray", # color of NA values
#           dendrogram="none", 
#           cexCol=0.75, # column text size
#           cexRow=heatmap.font.size.function(nrow(df.heatmap.basic)), # row text size
#           key.title="", # title of the key
#           key.xlab="", # key x axis labels
#           lmat=lmat, # requires heatmap layout above to be specified
#           lwid=lwid, # requires heatmap layout above to be specified
#           lhei=lhei, # requires heatmap layout above to be specified 
#           key.par = list(cex=1.05, mar=c(5,3,3.5,0)), # mar is A numerical vector of the form c(bottom, left, top, right) 
#           #which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
#           keysize = 1.5
# )
# # #add timestamp
# # mtext(date(), side=1, line=4, adj=0) # side (1=bottom, 2=left, 3=top, 4=right)
# graphics.off()
#-----------------------------------------------------------------------------------------------------


# END
