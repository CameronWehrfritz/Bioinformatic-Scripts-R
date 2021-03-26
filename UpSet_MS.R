#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#March 26, 2021


# Comparing several data sets
# Using {UpSetR} to visualize the sets and their interactions

# OUTPUT: UpSet plots

### Begin Script ###

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/0303_2021_AC1/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #
packages = c("writexl", "VennDiagram", "stringr", "purrr", "dplyr", "ggplot2", "tidyr", "UpSetR")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# LOAD DATA #

# macrophage library
df.macro <- read.csv("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/Input/AC1_0228_2021_v2_XX_MacrophLIB.tsv", sep="\t", stringsAsFactors = F) #PC
# assign name for this dataset

# mouse library
df.mouse <- read.csv("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/Input/AC1_0228_2021_v1_largelibr_MouseREF.tsv", sep="\t", stringsAsFactors = F) #PC

# AC1 library
df.ac1 <- read.csv("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/Input/AC1_0302_2021_v3_DDA_AC1_libr.tsv", sep="\t", stringsAsFactors = F) #PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Make list of data frames
df.list <- list(df.macro, df.mouse, df.ac1)
names(df.list) <- c("df.macro", "df.mouse", "df.ac1") # name each dataframe # do this manually
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Cleaning Function
# cleans the data

cleaning.function <- function(x){
  x <- x %>%
    select(-contains("DEPRECATED")) %>% # drop DEPRECATED variables
    select(-Valid) %>% # drop Valid variable
    rename(Number.of.Ratios=X..of.Ratios, # rename these variables
           Percent.Change=X..Change,
           Number.of.Unique.Total.Peptides=X..Unique.Total.Peptides,
           Number.of.Unique.Total.EG.Id=X..Unique.Total.EG.Id) %>%
    select(UniProtIds, Genes, # move these variables to the beginning 
           ProteinGroups, ProteinNames, ProteinDescriptions, 
           Comparison..group1.group2., 
           AVG.Log2.Ratio, Qvalue, Number.of.Unique.Total.Peptides, Number.of.Unique.Total.EG.Id,
           Absolute.AVG.Log2.Ratio,
           GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component, 
           Group, Organisms,
           everything()) %>% # all other variables follow
    relocate(Degrees.Of.Freedom, .after = last_col()) %>% # move to end
    relocate(Standard.Error, .after = last_col()) %>% # move to end
    filter(Number.of.Unique.Total.Peptides > 1) %>% # filter for at least 2 unique peptides per protein - exclude one peptide wonders
    arrange(Qvalue) # arrange by Qvalue, ascending
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Clean Data

# Apply Cleaning function over df.list
df.list.cleaned <- lapply(df.list, cleaning.function) 

# use list2env to assign each object in the list to an object in the global environment
list2env(df.list.cleaned, envir = .GlobalEnv)

# # make list of sets of UniProtIds
# sets.list <- list(MacrophageLib = df.macro$UniProtIds,
#                   MouseLib = df.mouse$UniProtIds,
#                   AC1Lib = df.ac1$UniProtIds
# )
# #print names
# names(sets.list)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# Which comparisons are present in the data?

# inner join unique comparisons from each dataset
comparisons <- df.ac1 %>% 
  select(Comparison..group1.group2.) %>% 
  unique() %>%
  inner_join(df.macro %>% select(Comparison..group1.group2.) %>% unique(), by="Comparison..group1.group2.") %>% 
  inner_join(df.mouse %>% select(Comparison..group1.group2.) %>% unique(), by="Comparison..group1.group2.") %>%
  pull()
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# make an UpSet Plot for each comparison - showing all proteins ID'd 
for(i in comparisons){
  print(i)
  # filter for one comparison
  df.list.singlecomparison <- lapply(df.list.cleaned, function(x){x <- x %>% filter(Comparison..group1.group2.==i)})
  
  # use list2env to assign each object in the list to an object into the environment
  list2env(df.list.singlecomparison, envir = environment())
  
  # make list of sets of UniProtIds
  sets.list <- list(MacrophageLib = df.macro$UniProtIds,
                    MouseLib = df.mouse$UniProtIds,
                    AC1Lib = df.ac1$UniProtIds
  )
  print(lengths(sets.list)) # make sure these numbers for each dataset looks right
  
  # make list of sets of UniProtIds
  sets.list <- list(MacrophageLib = df.macro$UniProtIds,
                    MouseLib = df.mouse$UniProtIds,
                    AC1Lib = df.ac1$UniProtIds
  )
  
  # UpSet Plot
  # ordered by frequency
  # x and y bar labels
  # custom text scaling
  filename <- paste("UpSet_comparison_", gsub(" / ", "_vs_", i), sep="") # first create file name
  # open png file
  png(filename = paste(filename, ".png", sep=""), # add .png to the filename
      width = 800, height = 800, units = "px", pointsize = 20,
      bg = "white")
  p <- upset(data = fromList(sets.list),
             nsets = length(sets.list), # plot all datasets in sets.list
             order.by = "freq", 
             mainbar.y.label = paste("Number of Proteins in", i, sep= " "), # mention the comparison here
             sets.x.label = "All Proteins per Library", # mention 'all proteins' here since upset doesn't allow for a title
             text.scale = c(2.0, # intersection size title
                            1.5, # intersection size tick labels
                            1.4, # set size title
                            1.5, # set size tick labels1.5, # set names
                            1) # numbers above bars
  )
  # inside of a loop one must print the graphical figure for it to save
  print(p)
  # print(p + theme_bw(base_size=16*(800/800)))
  # close png file
  dev.off()
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# make an UpSet Plot for each comparison - showing upregulated significantly changed proteins
for(i in comparisons){
  print(i)
  # filter for one comparison and for significance (upregulated)
  df.list.singlecomparison <- lapply(df.list.cleaned, 
                                     function(x){x <- x %>% 
                                       filter(Comparison..group1.group2.==i) %>%
                                       filter(Qvalue<0.01 & AVG.Log2.Ratio > 0.58)}) # upregulated
  
  # use list2env to assign each object in the list to an object into the environment
  list2env(df.list.singlecomparison, envir = environment())
  
  # make list of sets of UniProtIds
  sets.list <- list(MacrophageLib = df.macro$UniProtIds,
                    MouseLib = df.mouse$UniProtIds,
                    AC1Lib = df.ac1$UniProtIds
  )
  print(lengths(sets.list)) # make sure these numbers for each dataset looks right
  
  # make list of sets of UniProtIds
  sets.list <- list(MacrophageLib = df.macro$UniProtIds,
                    MouseLib = df.mouse$UniProtIds,
                    AC1Lib = df.ac1$UniProtIds
  )
  
  # UpSet Plot
  # ordered by frequency
  # x and y bar labels
  # custom text scaling
  filename <- paste("UpSet_up_significant_comparison_", gsub(" / ", "_vs_", i), sep="") # first create file name
  # open png file
  png(filename = paste(filename, ".png", sep=""), # add .png to the filename
      width = 800, height = 800, units = "px", pointsize = 20,
      bg = "white")
  p <- upset(data = fromList(sets.list),
             nsets = length(sets.list), # plot all datasets in sets.list
             order.by = "freq", 
             mainbar.y.label = paste("Number of Significantly Upregulated Proteins in", i, sep= " "), # mention the comparison here
             sets.x.label = "Sig. Proteins per Library", # mention 'significant proteins' here since upset doesn't allow for a title
             text.scale = c(2.0, # intersection size title
                            1.5, # intersection size tick labels
                            1.4, # set size title
                            1.5, # set size tick labels1.5, # set names
                            1) # numbers above bars
  )
  # inside of a loop one must print the graphical figure for it to save
  print(p)
  # print(p + theme_bw(base_size=16*(800/800)))
  # close png file
  dev.off()
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# make an UpSet Plot for each comparison - showing downregulated significantly changed proteins
for(i in comparisons){
  print(i)
  # filter for one comparison and for significance (downregulated)
  df.list.singlecomparison <- lapply(df.list.cleaned, 
                                     function(x){x <- x %>% 
                                       filter(Comparison..group1.group2.==i) %>%
                                       filter(Qvalue<0.01 & AVG.Log2.Ratio < -0.58)}) # downregulated
  
  # use list2env to assign each object in the list to an object into the environment
  list2env(df.list.singlecomparison, envir = environment())
  
  # make list of sets of UniProtIds
  sets.list <- list(MacrophageLib = df.macro$UniProtIds,
                    MouseLib = df.mouse$UniProtIds,
                    AC1Lib = df.ac1$UniProtIds
  )
  print(lengths(sets.list)) # make sure these numbers for each dataset looks right
  
  # make list of sets of UniProtIds
  sets.list <- list(MacrophageLib = df.macro$UniProtIds,
                    MouseLib = df.mouse$UniProtIds,
                    AC1Lib = df.ac1$UniProtIds
  )
  
  # UpSet Plot
  # ordered by frequency
  # x and y bar labels
  # custom text scaling
  filename <- paste("UpSet_significant_down_comparison_", gsub(" / ", "_vs_", i), sep="") # first create file name
  # open png file
  png(filename = paste(filename, ".png", sep=""), # add .png to the filename
      width = 800, height = 800, units = "px", pointsize = 20,
      bg = "white")
  p <- upset(data = fromList(sets.list),
             nsets = length(sets.list), # plot all datasets in sets.list
             order.by = "freq", 
             mainbar.y.label = paste("Number of Significantly Downregulated Proteins in", i, sep= " "), # mention the comparison here
             sets.x.label = "Sig. Proteins per Library", # mention 'significant proteins' here since upset doesn't allow for a title
             text.scale = c(2.0, # intersection size title
                            1.5, # intersection size tick labels
                            1.4, # set size title
                            1.5, # set size tick labels1.5, # set names
                            1) # numbers above bars
  )
  # inside of a loop one must print the graphical figure for it to save
  print(p)
  # print(p + theme_bw(base_size=16*(800/800)))
  # close png file
  dev.off()
}
#------------------------------------------------------------------------------------


# END
