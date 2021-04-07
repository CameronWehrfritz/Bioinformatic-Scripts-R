#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#April 5, 20201

# Performs statistical analysis on CETSA mass spec data via Individual Temperature Method
# This script calculates:
# i. Statistics via Welch Two Sample t-test for each temperature point
# ii. qvalues using bioconductor's qvalue package

# OUTPUT: 
# i. processed excel workbook with multiple sheets
# ii. plots of melting curves on the protein level - averaged across replicates

### Begin Script ###

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2021_0330_BRC4/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0330_BRC4/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("stringr", "hablar",
             "gdata",
             "readxl", "writexl", "openxlsx",
             "dplyr", "tidyr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# # to generate qvalues
# BiocManager::install("qvalue")
# suppressPackageStartupMessages(library(qvalue))
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# load data

# BGS Report of CETSA data on the peptide level (.tsv file)
df.input <- read.csv("//bigrock/GibsonLab/users/birgit/Barecia/BRC4_CDK/21_0325_BRC4_V1/Spectronaut/20210325_174612_210325_BRC4_all files_V1/20210401_164439_210325_BRC4_all_files_V1_BGS_V2_Report.tsv",
                     sep="\t", stringsAsFactors = FALSE)

# check out the unique files names
df.input$R.FileName %>% unique()

# how many files names?
df.input$R.FileName %>% unique() %>% length()

# candidates file
df.candidates <- read.csv("//bigrock/GibsonLab/users/birgit/Barecia/BRC4_CDK/21_0325_BRC4_V1/Spectronaut/20210325_174612_210325_BRC4_all files_V1/Candidates.tsv", sep="\t", stringsAsFactors = FALSE)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# remove one peptide wonders

# one peptide wonders
wonder.prots <- df.candidates %>%
  filter(X..Unique.Total.Peptides == 1) %>%
  pull(UniProtIds) %>%
  unique()

# remove one peptide wonders from input data
df.input <- df.input %>%
  filter(! PG.ProteinAccessions %in% wonder.prots) # exclude one peptide wonders

# remove one peptide wonders from candidates data
df.candidates <- df.candidates %>%
  filter(X..Unique.Total.Peptides > 1) # exclude one peptide wonders
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# clean data

# candidates file
df.candidates.clean <- df.candidates %>%
  select(-contains("DEPRECATED")) %>% # drop DEPRECATED variables
  select(-Valid) %>% # drop Valid variable
  rename(Number.of.Ratios=X..of.Ratios, # rename these variables
         Percent.Change=X..Change,
         Number.of.Unique.Total.Peptides=X..Unique.Total.Peptides,
         Number.of.Unique.Total.EG.Id=X..Unique.Total.EG.Id) %>%
  separate(col = Condition.Numerator, into = c("Treatment.Numerator", "Temperature.Numerator"), sep="-") %>% # NUMERATOR treatment and temperature separated by dash
  separate(col = Condition.Denominator, into = c("Treatment.Denominator", "Temperature.Denominator"), sep="-") %>% # DENOMINATOR treatment and temperature separated by dash
  filter(Temperature.Numerator==Temperature.Denominator) %>% # keep comparisons at the same temperature
  mutate(Temperature = Temperature.Numerator %>% as.numeric()) %>% # make temperature variable (should be same in numerator and denominator)
  select(UniProtIds, Genes, # move these variables to the beginning 
         ProteinGroups, ProteinNames, ProteinDescriptions,
         Comparison..group1.group2., Temperature,
         AVG.Log2.Ratio, Qvalue, Absolute.AVG.Log2.Ratio,
         GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component,
         everything()) 
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# prepare data

# BGS Report
# split treatment group and temperature from R.Condition variable
df <- df.input %>%
  separate(col = R.Condition, into = c("Treatment", "Temperature"), sep="-") %>% # treatment and temperature separated by dash
  mutate(Temperature = Temperature %>% as.numeric()) %>% # convert temperature variable to numeric
  group_by(PG.ProteinAccessions, PG.Genes, R.Replicate, Treatment, Temperature) %>% # sum FG.Quantity across peptides
  summarise(Total.FG.Quantity = sum(FG.Quantity)) %>% 
  ungroup() 
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# check out data a bit before proceeding with analysis

# how many temperature points and proteins were measured for each condition and replicate?
df %>%
  group_by(R.Replicate, Treatment) %>%
  mutate(Number.of.Temperature.Points = unique(Temperature) %>% length()) %>%
  mutate(Number.of.Proteins = unique(PG.ProteinAccessions) %>% length()) %>%
  select(R.Replicate, Treatment, Number.of.Temperature.Points, Number.of.Proteins) %>%
  unique()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# normalize to first temperature point

df <- df %>%
  group_by(PG.ProteinAccessions, PG.Genes, R.Replicate, Treatment) %>%
  mutate(Normalized.FG.Quantity=Total.FG.Quantity/Total.FG.Quantity[1]) %>%
  ungroup()
#------------------------------------------------------------------------------------


# WORKING ON THIS CHUNK 4/5/2021
#------------------------------------------------------------------------------------
# Calculate Statistical Differences between Treatment and Control at each Temperature Point

# this also works, and keeps all of the data
# Pvalues and T.values are hardcoded to p=1 and t=0 at 37 degrees (the first temp.) because those values are all normalized to 1
df.stats <- df %>%
  group_by(PG.ProteinAccessions, PG.Genes, Treatment, Temperature) %>%
  mutate(values=list(Normalized.FG.Quantity)) %>% # create list-column of quantities to perform statistics on
  spread(Treatment, values) %>% # spread values for each treatment group to the same row
  ungroup() %>%
  group_by(PG.ProteinAccessions, PG.Genes, Temperature) %>% 
  # perform statistical test by protein+gene+temperature
  mutate(Pvalue=ifelse(Temperature>37, t.test(unlist(CTL), unlist(NMN))$p.value, 1), 
         T.value=ifelse(Temperature>37, t.test(unlist(CTL), unlist(NMN))$statistic, 0)) %>%
  ungroup()

# # using pivot_wider ... not working yet - update later at some point
# df.test <- df %>%
#   group_by(PG.ProteinAccessions, PG.Genes, Treatment, Temperature) %>%
#   mutate(values=list(Normalized.FG.Quantity)) %>% # create list-column of quantities to perform statistics on
#   pivot_wider(names_from = Treatment, values_from = Normalized.FG.Quantity) %>% # spread values for each treatment group to the same row
#   ungroup() %>%
#   group_by(PG.ProteinAccessions, PG.Genes, Temperature) %>% 
#   # perform statistical test by protein+gene+temperature
#   mutate(Pvalue=ifelse(Temperature>37, t.test(unlist(CTL), unlist(NMN))$p.value, 1), 
#          T.value=ifelse(Temperature>37, t.test(unlist(CTL), unlist(NMN))$statistic, 0)) %>%
#   ungroup()

# Make Average value columns for plotting
df.stats <- df.stats %>%
  group_by(PG.ProteinAccessions, PG.Genes, Temperature) %>% 
  mutate(AVG.CONTROL=mean(unlist(CTL)), AVG.NMN=mean(unlist(NMN))) %>% # calculate averages
  mutate(Ratio=AVG.NMN/AVG.CONTROL) %>% # calculate ratio of averages
  ungroup()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Calculate the number of significant data points for each protein
# set threshold
Pvalue.threshold <- 0.05

df.sig.number <- df.stats %>%
  group_by(PG.ProteinAccessions) %>%
  summarise(Number.Significant.Points=sum(Pvalue < Pvalue.threshold & Ratio > 1)) # significant Pvalue and Ratio>1

# add this calculation on to df.stats
df.stats <- df.stats %>%
  group_by(PG.ProteinAccessions) %>%
  mutate(Number.Significant.Points=sum(Pvalue < Pvalue.threshold & Ratio > 1)) %>% # significant Pvalue and Ratio>1
  ungroup()
#------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------
# Which Proteins have at least 1 significant data point?

# statistically significant proteins
sig.proteins <- df.sig.number %>%
  filter(Number.Significant.Points>0) %>% # must have at least 1 significant data point: significant Pvalue and Ratio>1
  arrange(desc(Number.Significant.Points)) %>% # arrange most significant first
  pull(PG.ProteinAccessions) %>%
  unique()

# get statistically significant proteins
df.stats.sig <- df.stats %>%
  filter(Number.Significant.Points>0) %>%
  arrange(desc(Number.Significant.Points)) # arrange most significant first
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# prepare for writing out and plotting

# remove list-columns and convert to data frame - this makes writing out and plotting easier
df.stats.sig <- df.stats.sig %>%
  select(-CTL, -NMN) %>% # remove list-column
  as.data.frame()

df.stats.out <- df.stats %>%
  select(-CTL, -NMN) %>% # remove list-column
  as.data.frame()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# write out

# full dataset
# statistically significant proteins
write_xlsx(list("All" = df.stats.out, "Significant Hits" = df.stats.sig), path = "CETSA_Individual_Temperature_significant_proteins.xlsx")
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# plot

# define temperatures - for the x-axis
temps <- df.stats.sig %>%
  pull(Temperature) %>%
  unique() %>%
  as.numeric() # convert to numeric

#inititate PDF
pdf(file="CETSA_Individual_Temperature_Significant_test.pdf")
par(mfrow=c(2,3))

for(i in seq_along(unique(df.stats.sig$PG.ProteinAccessions))){
  # subset data for specific gene 
  df.loop <- df.stats.sig %>%
    filter(PG.ProteinAccessions==sig.proteins[i])
  
  # plot grid
  plot(0, xlim=c(min(temps), max(temps)),
       ylim=c(0, max(c(df.loop$AVG.CONTROL, df.loop$AVG.NMN))*1.4), 
       xlab=expression("Temperature " (degree*C)), 
       ylab="Relative Abundance", 
       main=unique(df.loop$PG.Genes)
  )
  
  # plot data points
  points(x=df.loop$Temperature %>% unique(), y=df.loop$AVG.CONTROL %>% unique(), col="blue") # control - blue
  points(x=df.loop$Temperature %>% unique(), y=df.loop$AVG.NMN %>% unique(), col="red") # treatment - red
  
  # add lines connecting the data points
  lines(x=df.loop$Temperature %>% unique(), y=df.loop$AVG.CONTROL %>% unique(), col="blue") # control - blue
  lines(x=df.loop$Temperature %>% unique(), y=df.loop$AVG.NMN %>% unique(), col="red") # treatment - red
  
  # flag the significant data points
  significant.points.indx <- which(df.loop$Pvalue < Pvalue.threshold) # grab indeces for where to flag data with asterisk
  # add flags
  points(x=temps[significant.points.indx], 
         y=rep(max(df.loop$AVG.CONTROL[significant.points.indx], df.loop$AVG.NMN[significant.points.indx])*1.15, length=length(significant.points.indx)),
         pch=8) # add significance flags (asterisks)
  
  # legend
  legend("topleft", legend=c("Treatment", "Control"), col=c("Red", "Blue"), lty = 1:1, cex=0.8)
} #end for loop
# mtext(date(), side=2, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right) # timestamp
graphics.off()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# summarize key points for each protein
# need to improve this...

df.test <- df.stats.sig %>%
  select(PG.ProteinAccessions, PG.Genes, Temperature, Pvalue, T.value, Number.Significant.Points) %>%
  unique()

df.test.2 <- df.stats.sig %>%
  select(PG.ProteinAccessions, PG.Genes, Number.Significant.Points) %>%
  unique()
#------------------------------------------------------------------------------------

# #------------------------------------------------------------------------------------
# # Qvalue 
# 
# qobj <- qvalue(p = df.stats$Pvalue)
# pi0 <- qobj$pi0
# 
# #adjust pvalues
# padj <- p.adjust(p = df.statistics$Pvalue, method=c("BH"), n = length(df.statistics$Pvalue))
# hist(padj, breaks=100)
# unique(padj)
# #------------------------------------------------------------------------------------


# END
