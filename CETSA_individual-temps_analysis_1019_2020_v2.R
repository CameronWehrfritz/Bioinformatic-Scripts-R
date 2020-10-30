#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#October 15, 2020


# BARECIA BC2 Analysis
# cetsa analysis 

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
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# LOAD DATA #

# # BGS Report with THREE replicates:
# df <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/20201013_111205_200909_BRC23_1009_2020_v1_BGS_Report.csv", sep=",", stringsAsFactors = F) #VPN mac
# df <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/20201013_111205_200909_BRC23_1009_2020_v1_BGS_Report.csv", sep=",", stringsAsFactors = F) #VPN windows
# # fix first column name, which is messed up for some reason
# names(df)[1] <- "R.Condition"


# BGS Report with FOUR replicates:
# df <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/20201015_121004_200909_BRC123_1012_2020_v3_BGS_Report.csv", sep=",", stringsAsFactors = F) #VPN mac
df <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/20201015_121004_200909_BRC123_1012_2020_v3_BGS_Report.csv", sep=",", stringsAsFactors = F) #VPN windows
# fix first column name, which is messed up for some reason
names(df)[1] <- "R.Condition"


# # candidates file:
# df.barecia <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/Candidates_BRC23_1009_2020_v1.tsv", sep="\t", stringsAsFactors = F) #VPN mac
df.barecia <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/Candidates_BRC23_1009_2020_v1.tsv", sep="\t", stringsAsFactors = F) #VPN mac


# Human Proteome UniProtIds and Gene Names
# df.human.proteome <- read.csv("/Volumes/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/uniprot_human_table.csv", stringsAsFactors = F)
df.human.proteome <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_1012_Barecia_BC2/Input/uniprot_human_table.csv", stringsAsFactors = F)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Keep only necessary columns
dat <- df %>%
  select(PG.ProteinGroups, R.Condition, R.Replicate, FG.Quantity) %>%
  separate(col=R.Condition, into=c("Temperature", "Treatment"), sep="_", remove=FALSE) # make two new columns
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# REDUCE DATA TABLE

# Calculate  Total FG.Quantity across same variables: Protein, Treatment, Temperature and Replicate

dat.total <- dat %>%
  group_by(PG.ProteinGroups, Treatment, Temperature, R.Replicate) %>%
  summarise(FG.Quantity=sum(FG.Quantity)) 

dat.total <- dat.total %>%
  group_by(PG.ProteinGroups, Treatment, R.Replicate) %>%
  mutate(Normalized.FG.Quantity=FG.Quantity/FG.Quantity[1])
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Add Human Proteome information by matching with UniProtIds
dat.total <- dat.total %>%
  left_join(df.human.proteome, by=c("PG.ProteinGroups"="Entry"))

# Make new Gene column
dat.total <- dat.total %>%
  mutate(Gene=unlist(str_split(Gene.names, " "))[1]) %>%
  select(PG.ProteinGroups, Gene, everything()) # move Gene column left
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Calculate Statistical Differences between Treatment and Control at each Temperature Point

# this also works, and keeps all of the data
# P.values and T.values are hardcoded to p=1 and t=0 at 37 degrees because the values are all normalized to 1 there
df.stats <- dat.total %>%
  group_by(PG.ProteinGroups, Gene, Treatment, Temperature) %>%
  summarise(values=list(Normalized.FG.Quantity)) %>% # create list-column
  spread(Treatment, values) %>% # spread values for each Treatment group to the same row
  group_by(PG.ProteinGroups, Gene, Temperature) %>%
  mutate(P.value=ifelse(Temperature>37, t.test(unlist(contr), unlist(NMN))$p.value, 1), 
         T.value=ifelse(Temperature>37, t.test(unlist(contr), unlist(NMN))$statistic, 0))

# Make Average value columns for plotting
df.stats <- df.stats %>%
  mutate(AVG.contr=mean(unlist(contr)), AVG.NMN=mean(unlist(NMN))) %>% # calculate averages
  mutate(Ratio=AVG.NMN/AVG.contr) # calculate ratio of averages
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Calculate the number of significant data points for each protein
P.threshold <- 0.05

df.sig.number <- df.stats %>%
  group_by(PG.ProteinGroups) %>%
  summarise(Number.Significant.Points=sum(P.value<P.threshold & Ratio>1)) # significant Pvalue and Ratio>1
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Which Proteins have at least 1 significant data point?

# vector for plotting statistically significant proteins (Pvalue filter)
sig.proteins <- df.sig.number %>%
  filter(Number.Significant.Points>0) %>% # must have at least 1 significant data point: significant Pvalue and Ratio>1
  pull(PG.ProteinGroups) %>%
  unique()


# subset df.stats for these significant proteins
df.stats.sig <- df.stats %>%
  filter(PG.ProteinGroups %in% sig.proteins) # statistically significant filter - Pvalue


# remove list-columns and convert to data frame - this makes writing out and plotting easier
df.stats.sig <- df.stats.sig %>%
  select(-contr, -NMN) %>% 
  as.data.frame() %>%
  left_join(df.human.proteome, by=c("PG.ProteinGroups"="Entry"))


# apply additional filter: at the last temperature point the Ratio (Treatment/Control) must be greater than 1
last.temp.sig.proteins <- df.stats.sig %>%
  filter(Temperature==70 & Ratio>1) %>%
  pull(PG.ProteinGroups) %>%
  unique()

# filter for these last-temperature significant proteins
df.stats.sig <- df.stats.sig %>%
  filter(PG.ProteinGroups %in% last.temp.sig.proteins) # last temperature point Ratio filter


# write out data frame of statistically significant proteins, which we will plot
write.csv(df.stats.sig, file="Table_significant-proteins_test.csv", row.names = FALSE)

# temperatures - for the x-axis
temps <- df.stats.sig$Temperature %>%
  unique() %>%
  as.numeric() # convert to numeric

# Plotting:
#inititate PDF
pdf(file="Plots_test.pdf")
par(mfrow=c(2,3))

for(i in seq_along(unique(df.stats.sig$PG.ProteinGroups))){
  # subset data for specific gene 
  df.loop <- df.stats.sig %>%
    filter(PG.ProteinGroups==last.temp.sig.proteins[i])
  
  # plot grid
  plot(0, xlim=c(min(temps), max(temps)),
       ylim=c(0, max(c(df.loop$AVG.contr, df.loop$AVG.NMN))*1.4), 
       xlab="Temperature (Â°C)", 
       ylab="Relative Abundance", 
       main=unique(df.loop$Gene)
  )

  # plot data points
  points(x=df.loop$Temperature, y=df.loop$AVG.contr, col="blue") # control - blue
  points(x=df.loop$Temperature, y=df.loop$AVG.NMN, col="red") # treatment - red
  
  # flag the significant data points
  significant.points.indx <- which(df.loop$P.value<P.threshold) # indeces for where to flag data with asterisk
  points(x=temps[significant.points.indx], 
         y=rep(max(df.loop$AVG.contr[significant.points.indx], df.loop$AVG.NMN[significant.points.indx])*1.15, length=length(significant.points.indx)),
         pch=8) # add significance flags (asterisks)
  
  #add a line connecting the data points
  lines(x=temps, y=df.loop$AVG.contr, col="blue") # control
  lines(x=temps, y=df.loop$AVG.NMN, col="red") # treatment
  
  #Legend
  legend("topleft", legend=c("Treatment", "Control"), col=c("Red", "Blue"), lty = 1:1, cex=0.8)
} #end for loop
# mtext(date(), side=2, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right) # timestamp
graphics.off()
#------------------------------------------------------------------------------------


# #------------------------------------------------------------------------------------
# # How many of our significant results are membrane or plasma membrane associated?
# 
# # Four Reps
# df.membrane <- df.barecia %>%
#   filter(ProteinGroups %in% last.temp.sig.proteins) %>% # significant proteins
#   filter(grepl("membrane", GO.Cellular.Component)) %>% # membrane associated
#   filter(Comparison..group1.group2.=="70_NMN / 70_contr") # take only this row, shortens the data frame
# 
# write.csv(df.membrane, file="Table_membrane-associated_significant-proteins_test.csv", row.names = FALSE)
# 
# # Three Reps
# df.membrane <- df.barecia %>%
#   filter(ProteinGroups %in% unique(df.three.reps$PG.ProteinGroups)) %>% # significant proteins
#   filter(grepl("membrane", GO.Cellular.Component)) %>% # membrane associated
#   filter(Comparison..group1.group2.=="70_NMN / 70_contr") # take only this row, shortens the data frame
# 
# 
# # combined four and three reps
# df.membrane <- df.barecia %>%
#   filter(ProteinGroups %in% unique(df.three.reps$PG.ProteinGroups) | ProteinGroups %in% last.temp.sig.proteins) %>% # significant proteins
#   filter(grepl("membrane", GO.Cellular.Component)) %>% # membrane associated
#   filter(Comparison..group1.group2.=="70_NMN / 70_contr") %>% # take only this row, shortens the data frame
#   select(ProteinGroups, Genes, GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component, everything())
# 
# write.csv(df.membrane, file="Table_membrane-associated_significant-proteins_test.csv", row.names = FALSE)
# #------------------------------------------------------------------------------------


# #------------------------------------------------------------------------------------
# # combine 3replicate and 4replicate data sets
# df.combined <- df.stats.sig %>%
#   filter(PG.ProteinGroups %in% overlap.proteins) %>%
#   cbind(Number.Replicates=4) %>%
#   rbind(df.three.reps %>% filter(PG.ProteinGroups %in% overlap.proteins) %>% cbind(Number.Replicates=3)) %>%
#   arrange(PG.ProteinGroups)
# 
# # write table
# write.csv(df.combined, file="Table_significant-proteins_test.csv", row.names = FALSE)
# 
# 
# # Plot
# 
# pdf(file="Plots_test.pdf")
# par(mfrow=c(2,3))
# 
# reps <- unique(df.combined$Number.Replicates)
# 
# for(i in seq_along(unique(df.combined$PG.ProteinGroups))){
#   for(j in seq_along(unique(df.combined$Number.Replicates))){
#     # subset data for specific gene 
#     df.loop <- df.combined %>%
#       filter(PG.ProteinGroups==overlap.proteins[i]) %>%
#       filter(Number.Replicates==reps[j])
#     
#     # plot grid
#     plot(0, xlim=c(min(temps), max(temps)),
#          ylim=c(0, max(c(df.loop$AVG.contr, df.loop$AVG.NMN))*1.4), 
#          xlab="Temperature (°C)", 
#          ylab="Relative Abundance", 
#          main=paste(unique(df.loop$Gene), "# Replicates=",unique(df.loop$Number.Replicates))
#     )
#     
#     # plot data points
#     points(x=df.loop$Temperature, y=df.loop$AVG.contr, col="blue") # control - blue
#     points(x=df.loop$Temperature, y=df.loop$AVG.NMN, col="red") # treatment - red
#     
#     # flag the significant data points
#     significant.points.indx <- which(df.loop$P.value<P.threshold) # indeces for where to flag data with asterisk
#     points(x=temps[significant.points.indx], 
#            y=rep(max(df.loop$AVG.contr[significant.points.indx], df.loop$AVG.NMN[significant.points.indx])*1.15, length=length(significant.points.indx)),
#            pch=8) # add significance flags (asterisks)
#     
#     #add a line connecting the data points
#     lines(x=temps, y=df.loop$AVG.contr, col="blue") # control
#     lines(x=temps, y=df.loop$AVG.NMN, col="red") # treatment
#     
#     #Legend
#     legend("topleft", legend=c("Treatment", "Control"), col=c("Red", "Blue"), lty = 1:2, cex=0.8)
#   } #end for reps loop
# } #end for protein loop
# # mtext(date(), side=2, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right) # timestamp
# graphics.off()
# #------------------------------------------------------------------------------------


# END #