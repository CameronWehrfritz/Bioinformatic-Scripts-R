#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#October 16, 2020


# BARECIA BC2 Analysis

######################
#### Begin Program ###
######################

#------------------------------------------------------------------------------------
#set working directory
setwd("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/R_workspace") # VPN mac
# setwd("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/R_workspace") # VPN windows
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("gplots", "dplyr", "svMisc", "tidyr", "stringr", "reshape", "reshape2", "data.table", "hablar")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# # nate uses this for qvalues
# BiocManager::install("qvalue")
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# read in data for filtering and plotting:

# # data with THREE replicates:
# # data for filtering
# df.stats <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1013_2020_3replicates/Table_Barecia_BC2_all-proteins_1013_2020_v1.csv", stringsAsFactors = FALSE)
# df.stats <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1013_2020_3replicates/Table_Barecia_BC2_all-proteins_1013_2020_v1.csv", stringsAsFactors = FALSE)
# # # data for plotting
# # df.areas <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1013_2020_3replicates/Table_Barecia_BC2_areas_1013_2020_v1.csv", stringsAsFactors = FALSE)
# df.areas <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1013_2020_3replicates/Table_Barecia_BC2_areas_1013_2020_v1.csv", stringsAsFactors = FALSE)
# df.areas <- df.areas %>%
#   select(-X) # remove unnecessary column


# # data with FOUR replicates:
# df.stats <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1015_2020_4replicates/Table_Barecia_4replicates_all-proteins_1015_2020_v1.csv", stringsAsFactors = FALSE)
df.stats <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1015_2020_4replicates/Table_Barecia_4replicates_all-proteins_1015_2020_v1.csv", stringsAsFactors = FALSE)
# data for plotting
# df.areas <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1015_2020_4replicates/Table_Barecia_4replicates_areas_1015_2020_v1.csv", stringsAsFactors = FALSE)
df.areas <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Results/Tables/Area-method/1015_2020_4replicates/Table_Barecia_4replicates_areas_1015_2020_v1.csv", stringsAsFactors = FALSE)
df.areas <- df.areas %>%
  select(-X) # remove unnecessary column
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Filters
# Significant P<0.05 and Ratio (Treatment/Control)>1 (simultaneously)

df.stats.sig <- df.stats %>%
  filter(Pvalue<0.05 & Ratio>1) %>%
  arrange(Pvalue)

# write out summary report
write.csv(df.stats.sig, file="Table_summary_report_date.csv", row.names = FALSE)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Plots

# filter for significant proteins
df.areas.sig <- df.areas %>% 
  filter(Genes %in% df.stats.sig$Genes)

# calculate the mean of FG.Quantity across protein/gene, treatment and temperature
# this averages across the replicates
df.areas.sig.avg <- df.areas.sig %>%
  group_by(Genes, PG.ProteinGroups, Treatment, Temperature) %>%
  summarise(FG.Quantity=mean(FG.Quantity))

# relevel factor in descending order of significance
df.areas.sig.avg$Genes <- factor(df.areas.sig.avg$Genes, levels=df.stats.sig$Genes)

# get vector of genes from df.stats.sig in ascending order of Pvalue
sig.genes <- df.stats.sig$Genes

#inititate PDF
pdf(file="Plots_CETSA_test.pdf")
par(mfrow=c(2,3))

# gene - for looping through proteins
gene <- unique(df.areas.sig.avg$Genes)

# treatment conditions
treatments <- unique(df.areas.sig.avg$Treatment)

# temperatures, for x-axis
x <- unique(df.areas.sig.avg$Temperature)

for(i in 1:length(gene)){
  # subset data for one protein
  df.loop <- df.areas.sig.avg %>%
    filter(Genes==gene[i])
  
  # subset for the treatment conditions
  # control:
  df.ctrl <- df.loop %>%
    filter(Treatment==treatments[1]) %>%
    mutate(Norm.FG.Quantity=FG.Quantity/FG.Quantity[1]) # normalize to first data point, at lowest temp.
  # treatment:
  df.treat <- df.loop %>%
    filter(Treatment==treatments[2]) %>%
    mutate(Norm.FG.Quantity=FG.Quantity/FG.Quantity[1]) # normalize to first data point, at lowest temp.

  # plot
  plot(0,
       xlim=c(min(x), max(x)),
       ylim=c(0,max(c(df.treat$Norm.FG.Quantity, df.ctrl$Norm.FG.Quantity))*1.4), 
       xlab="Temperature (°C)", 
       ylab="Relative Abundance", 
       main = unique(df.loop$Genes)
       )
  # plot data points
  with(df.ctrl, points(Temperature, Norm.FG.Quantity, col="blue")) # control data
  with(df.treat, points(Temperature, Norm.FG.Quantity, col="red")) # treatment data
  
  # plot line connecting the data points
  with(df.ctrl, lines(Temperature, Norm.FG.Quantity, col="blue")) # control data
  with(df.treat, lines(Temperature, Norm.FG.Quantity, col="red")) # plot treatment data
  
  # legend
  legend("topleft", legend=c("Treatment", "Control"), col=c("Red", "Blue"), lty = 1:1, cex=0.8)
} #end for loop
# mtext(date(), side=1, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right) # timestamp
graphics.off()
#------------------------------------------------------------------------------------



# END
