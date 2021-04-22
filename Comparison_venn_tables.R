#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#April 21, 2021

# Comparison between two protein lists
# Generates two-way venn diagram and corresponding tables

### Begin Script ###

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/0409_2021_BARECIA/2021_0421_Barecia_Comparisons/BRC4_vs_BRC5/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/0409_2021_BARECIA/2021_0421_Barecia_Comparisons/BRC4_vs_BRC5/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("stringr", "hablar", "VennDiagram",
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

# BRC4 CETSA
df.one <- read_xlsx("//bigrock/GibsonLab/users/Cameron/0409_2021_BARECIA/2021_0330_Barecia_BRC4/Results/Individual_Temperature/Tables/CETSA_BRC4_Individual_Temperature_Results_0419_2021_v7.xlsx",
                     sheet=3, skip=5)

# BRC5 CETSA
df.two <- read_xlsx("//bigrock/GibsonLab/users/Cameron/0409_2021_BARECIA/2021_0419_Barecia_BRC5/Results/Individual_Temperature/Tables/CETSA_BRC5_Individual_Temperature_Results_0420_2021_v1.xlsx",
                     sheet=3, skip=5)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# venn of significantly altered proteins with sigmoidal melting curves

# Venn
venn.diagram(
  list("BRC4" = df.one$PG.ProteinAccessions, "BRC5" = df.two$PG.ProteinAccessions), # manually name the datasets here
  filename="Venn_test.png",
  main = "Significantly Altered Proteins \n P<0.05", # customize title
  fill = c("red", "blue"),
  main.cex = 1,
  ext.text = TRUE,
  #cat.pos = c(3, 8),
  cex = 1.5
  #inverted = TRUE
)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# generate tables

# unique to first dataset df.one
df.one.unique <- df.one %>%
  filter(! PG.ProteinAccessions %in% df.two$PG.ProteinAccessions)

# unique to second dataset df.two
df.two.unique <- df.two %>%
  filter(! PG.ProteinAccessions %in% df.one$PG.ProteinAccessions)

# overlap between both datasets
df.overlap <- df.one %>%
  rename("Rank.BRC4"=Rank, "Temperature.Points.BRC4"=`Temperature Points`) %>%
  inner_join(df.two %>% select("Rank.BRC5"=Rank, "Temperature.Points.BRC5"=Temperature.Points, PG.ProteinAccessions), by="PG.ProteinAccessions") 

# move some variables around - customized
df.overlap <- df.overlap %>%
  select(Rank.BRC4, Rank.BRC5, everything()) %>% # move rank to front
  relocate(Temperature.Points.BRC5, .after=Temperature.Points.BRC4) # move temperature points together
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# write tables as one excel workbook

write_xlsx(list(df.one.unique, df.two.unique, df.overlap), "Tables_Comparisons_test.xlsx")
#------------------------------------------------------------------------------------


# END
