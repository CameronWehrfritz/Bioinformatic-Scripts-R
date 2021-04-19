#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#April 5, 2021

# Performs statistical analysis on CETSA mass spec data via Individual-Temperature-Method
# This script calculates:
# i. Statistics via Welch Two Sample t-test at each temperature, comparing relative abundances of each treatment group
# ii. Rank of each protein as a potential target (higher the better)

# OUTPUT: 
# i. processed excel workbook with multiple sheets
# ii. plots of melting curves in ranked order


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


#------------------------------------------------------------------------------------
# Calculate Statistical Differences between Treatment and Control at each Temperature Point

df.stats <- df %>%
  group_by(PG.ProteinAccessions, PG.Genes, Treatment, Temperature) %>%
  mutate(values=list(Normalized.FG.Quantity)) %>% # create list-column of quantities to perform statistics on
  spread(Treatment, values) %>% # spread values for each treatment group to the same row
  ungroup() %>%
  group_by(PG.ProteinAccessions, PG.Genes, Temperature) %>% 
  # perform statistical test by protein+gene+temperature
  # Pvalues and T.values are hardcoded to Pvalue=1 and T.value=0 at 37 degrees (the first temperature) because those values are all normalized to 1
  mutate(Pvalue=ifelse(Temperature>37, t.test(unlist(CTL), unlist(NMN))$p.value, 1), 
         T.value=ifelse(Temperature>37, t.test(unlist(CTL), unlist(NMN))$statistic, 0)) %>%
  ungroup()

# # using pivot_wider ... not working yet - update later at some point?
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

# make average value columns for plotting
df.stats <- df.stats %>%
  group_by(PG.ProteinAccessions, PG.Genes, Temperature) %>% 
  mutate(AVG.CONTROL=mean(unlist(CTL)), AVG.NMN=mean(unlist(NMN))) %>% # calculate averages
  mutate(Ratio=AVG.NMN/AVG.CONTROL) %>% # calculate ratio of averages - ensure control is in the denominator
  ungroup()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Calculate the number of significant data points for each protein

# set pvalue threshold
Pvalue.threshold <- 0.05

# reduce observations and calculate number of significant data points
df.stats.reduced <- df.stats %>%
  # boil down observations
  select(PG.ProteinAccessions, PG.Genes, Temperature, AVG.CONTROL, AVG.NMN, Ratio, T.value, Pvalue) %>%
  unique() %>%
  # determine significant yes/no for each observation
  mutate(Significant = ifelse(Pvalue < Pvalue.threshold & Ratio > 1, "yes", "no")) %>%
  # calculate number of statistically significant data points for each protein
  group_by(PG.ProteinAccessions, PG.Genes) %>%
  mutate(Number.Significant.Points = sum(Pvalue < Pvalue.threshold & Ratio > 1)) %>%
  ungroup() %>%
  arrange(desc(Number.Significant.Points)) # arrange most significant first

# add this calculation on to the full data set for completeness
df.stats <- df.stats %>%
  left_join(df.stats.reduced %>% select(PG.ProteinAccessions, PG.Genes, Number.Significant.Points), by=c("PG.ProteinAccessions", "PG.Genes"))
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# which proteins have at least 1 significant data point?

# statistically significant proteins
sig.proteins <- df.stats.reduced %>%
  filter(Number.Significant.Points>0) %>% # must have at least 1 significant data point
  arrange(desc(Number.Significant.Points)) %>% # arrange most significant first
  pull(PG.ProteinAccessions) %>%
  unique()

# subset statistically significant protein data
df.significant <- df.stats.reduced %>%
  filter(Number.Significant.Points>0) %>% # must have at least 1 significant temperature point
  arrange(desc(Number.Significant.Points)) # arrange most significant first
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# identify two types of melting curves: sigmoidal and non-sigmoidal
# based on the control group
# sigmoidal = last two temperature points both have a relative abundance less than 1
# non-sigmoidal = last two temperature points have a relative abundance greater than or equal to 1
# these are of course only approximating the sigmoidal curve, since it doesn't explicitly check the earlier temperature points

# yes or no sigmoidal curve?
df.significant <- df.significant %>%
  group_by(PG.ProteinAccessions) %>%
  # is the relative abundance of the control group less than 1 at the last two temperature points?
  mutate(Sigmoidal = ifelse(Temperature==65 | Temperature==70,  #
                            ifelse(AVG.CONTROL<1, "yes", "no"), "no check")) %>%
  # are either of the last points non-sigmoidal? if so, the protein is not sigmoidal
  mutate(Sigmoidal = ifelse("no" %in% Sigmoidal, "no", "yes")) %>% # yes or no sigmoidal curve?
  ungroup()

# sigmoidal proteins
sigmoidal.proteins <- df.significant %>%
  filter(Sigmoidal == "yes") %>%
  pull(PG.ProteinAccessions) %>%
  unique()

# non-sigmoidal proteins
nonsigmoidal.proteins <- df.significant %>%
  filter(Sigmoidal == "no") %>%
  pull(PG.ProteinAccessions) %>%
  unique()
#------------------------------------------------------------------------------------                            


#------------------------------------------------------------------------------------
# rank the significant hits

# score according to the following criteria, in this order below: 
# 1. number of significant points (which must be in the right direction, ie stabilized compared to control) - 1 point each
# 2. number of hits at middle three temperatures (50, 55, 60) - 1 point each

# first calculate score
df.ranked <- df.significant %>%
  group_by(PG.ProteinAccessions, PG.Genes) %>%
  filter(Temperature %in% c(50, 55, 60)) %>% # keep data at middle three temperatures (50, 55, 60)
  mutate(N = sum(Significant=="yes")) %>% # count number of significant hits at middle three temperatures
  mutate(Score = Number.Significant.Points + N) # calculate score

# calculate and order by rank
df.ranked <- df.ranked %>%
  mutate(Rank = max(df.ranked$Score)-Score+1) %>% # calculate rank (so lower is better)
  select(PG.ProteinAccessions, PG.Genes, Number.Significant.Points, Score, Rank) %>%
  unique() %>%
  arrange(Rank) # ascending by rank - best at the top

# add rank to the full data set 
df.stats <- df.stats %>%
  left_join(df.ranked %>% select(PG.ProteinAccessions, PG.Genes, Rank), by=c("PG.ProteinAccessions", "PG.Genes")) %>%
  mutate(Rank = ifelse(is.na(Rank), max(df.ranked$Rank)+1, Rank)) %>% # fill in missing Rank from NA to one above the max
  arrange(Rank) # descending by rank - best at the top

# add rank to reduced data set
df.stats.reduced <- df.stats.reduced %>%
  left_join(df.ranked %>% select(PG.ProteinAccessions, PG.Genes, Rank), by=c("PG.ProteinAccessions", "PG.Genes")) %>%
  mutate(Rank = ifelse(is.na(Rank), max(df.ranked$Rank)+1, Rank)) %>% # fill in missing Rank from NA to one above the max
  arrange(Rank) # descending by rank - best at the top

# add rank to statistically significant protein data
df.significant <- df.significant %>%
  left_join(df.ranked %>% select(PG.ProteinAccessions, PG.Genes, Rank), by=c("PG.ProteinAccessions", "PG.Genes")) %>%
  arrange(Rank) # ascending by rank - best at the top
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# add description attributes to output files

# add attributes from candidates file
df.stats.reduced.out <- df.stats.reduced %>%
  left_join(df.candidates.clean %>%
              select(UniProtIds, ProteinDescriptions, Organisms, GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component), by=c("PG.ProteinAccessions"="UniProtIds")) %>%
  unique() %>%
  select(Rank, PG.ProteinAccessions, PG.Genes, ProteinDescriptions, everything())

# add attributes from candidates file
df.significant.out <- df.significant %>%
  left_join(df.candidates.clean %>%
              select(UniProtIds, ProteinDescriptions, Organisms, GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component), by=c("PG.ProteinAccessions"="UniProtIds")) %>%
  unique() %>%
  select(Rank, Sigmoidal, PG.ProteinAccessions, PG.Genes, ProteinDescriptions, everything())
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# create protein summary

df.protein.summary <- df.significant.out %>%
  select(Rank, Sigmoidal, PG.ProteinAccessions, PG.Genes, ProteinDescriptions, Number.Significant.Points, GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component, Organisms) %>%
  unique() %>%
  arrange(Rank) %>% # best rank at top
  arrange(desc(Sigmoidal)) # sigmoidal at top
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# create overall summary

df.summary <- df.stats.reduced %>%
  select(PG.ProteinAccessions, PG.Genes, Rank) %>%
  unique() %>%
  group_by(Rank) %>%
  summarise(Number.of.Proteins = n()) %>% # calculate count for each rank
  arrange(Rank) %>% # best at top
  ungroup() %>%
  select(Rank, Number.of.Proteins)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# write out results

# full dataset
# statistically significant proteins
write_xlsx(list("Summary" = df.summary, # overall summary
                "Protein Summary" = df.protein.summary, # protein level summary
                "Sigmoidal Protein Summary" = df.protein.summary %>% filter(PG.ProteinAccessions %in% sigmoidal.proteins), # protein level summary - sigmoidal curves
                "Non-Sigmoidal Protein Summary" = df.protein.summary %>% filter(PG.ProteinAccessions %in% nonsigmoidal.proteins), # protein level summary - non-sigmoidal curves
                "Significant Proteins" = df.significant.out, # significant proteins
                "All Proteins" = df.stats.reduced.out), # all proteins
           path = "CETSA_Individual_Temperature_significant_proteins.xlsx")
#------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------
# plot significant sigmoidal proteins

# define temperatures - for the x-axis
temps <- df.significant %>%
  pull(Temperature) %>%
  unique() %>%
  as.numeric() # convert to numeric

# sigmoidal proteins
df.significant.sigmoidal <- df.significant %>%
  filter(Sigmoidal == "yes") 

# inititate PDF
pdf(file="CETSA_Individual_Temperature_Sigmoidal_Proteins_test.pdf")
par(mfrow=c(2,3))

for(i in seq_along(unique(df.significant.sigmoidal$PG.ProteinAccessions))){
  # grab data for specific protein 
  df.loop <- df.significant.sigmoidal %>%
    filter(PG.ProteinAccessions==sigmoidal.proteins[i])
  
  # initiate plot
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
  significant.points.indeces <- which(df.loop$Pvalue < Pvalue.threshold & df.loop$Ratio > 1) # grab indeces for where to flag data with asterisk
  # add flags
  points(x=temps[significant.points.indeces], 
         y=rep(max(df.loop$AVG.CONTROL[significant.points.indeces], df.loop$AVG.NMN[significant.points.indeces])*1.15, length=length(significant.points.indeces)),
         pch=8) # add significance flags (asterisks)
  
  # legend
  legend("topleft", legend=c("CTL", "NMN"), col=c("Blue", "Red"), lty = 1:1, cex=0.8)
} #end for loop
# mtext(date(), side=2, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right) # timestamp
graphics.off()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# plot significant non-sigmoidal proteins

# define temperatures - for the x-axis
temps <- df.significant %>%
  pull(Temperature) %>%
  unique() %>%
  as.numeric() # convert to numeric

# sigmoidal proteins
df.significant.nonsigmoidal <- df.significant %>%
  filter(Sigmoidal == "no") 

# inititate PDF
pdf(file="CETSA_Individual_Temperature_NonSigmoidal_Proteins_test.pdf")
par(mfrow=c(2,3))

for(i in seq_along(unique(df.significant.nonsigmoidal$PG.ProteinAccessions))){
  # grab data for specific protein 
  df.loop <- df.significant.nonsigmoidal %>%
    filter(PG.ProteinAccessions==nonsigmoidal.proteins[i])
  
  # initiate plot
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
  significant.points.indeces <- which(df.loop$Pvalue < Pvalue.threshold & df.loop$Ratio > 1) # grab indeces for where to flag data with asterisk
  # add flags
  points(x=temps[significant.points.indeces], 
         y=rep(max(df.loop$AVG.CONTROL[significant.points.indeces], df.loop$AVG.NMN[significant.points.indeces])*1.15, length=length(significant.points.indeces)),
         pch=8) # add significance flags (asterisks)
  
  # legend
  legend("topleft", legend=c("CTL", "NMN"), col=c("Blue", "Red"), lty = 1:1, cex=0.8)
} #end for loop
# mtext(date(), side=2, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right) # timestamp
graphics.off()
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
