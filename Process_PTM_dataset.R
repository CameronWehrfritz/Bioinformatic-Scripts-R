#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#May 6, 2021

# DESCRIPTION: This script processes Mass Spectrometry PTM datasets
# i. cleans data
# ii. corrects certain comparisons (if they are in the inverted order - this must be specified by hand, see 'CORRECT SOME COMPARISONS' chunk below)
# iii. joins spectronaut and protein pilot datasets in order to retrieve peptide cleavage information from protein pilot dataset
# iv. write processed data to excel workbook with multiple sheets

# OUTPUT: PROCESSED EXCEL FILE

### Begin Script ###


#-----------------------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/R_workspace") # PC
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# PACKAGES #
packages = c("purrr", "stringr", "tidyr", "dplyr", "writexl", "openxlsx")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# LOAD DATA #

# protein pilot - contains the cleavage locations - used to determine if peptide is semi-tryptic 
df.pilot <- read.delim("//bigrock/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/Input/201203_IK1_DDA_ALL_PeptideSummary.txt", stringsAsFactors = F) # PC


# NMN datasets:
# # load candidates file
# df.input <- read.csv("//bigrock/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/Input/Candidates_IK1_pept_semitrypt_NMN_0504_2021_v5.tsv", sep="\t", stringsAsFactors = F) # PC
# 
# # define file name for output file
# input.filename <- "Candidates_IK1_pept_semitrypt_NMN_0504_2021_v5.tsv" # manually edit this
# 
# # load peptide data
# df.peptide <- read.xlsx("//bigrock/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/Input/20210504_170403_201203_IK1_pept_semitrypt_NMN_0504_2021_v5_Report.xlsx")


# NAD BAD datasets:
# load candidates file
df.input <- read.csv("//bigrock/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/Input/Candidates_IK1_pept_semitrypt_NAD_BAD_0504_2021_v4.tsv", sep="\t", stringsAsFactors = F) # PC

# define file name for output file
input.filename <- "Candidates_IK1_pept_semitrypt_NAD_BAD_0504_2021_v4.tsv" # manually edit this

# load peptide data
df.peptide <- read.xlsx("//bigrock/GibsonLab/users/Cameron/2021_0503_LimitedProteolysis/Input/report_NAD_BAD__0504_2021_v4.xlsx")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# define numerical values

# significance criteria
# Qvalue.threshold <- 0.01 
Pvalue.threshold <- 0.01 # using pvalue threshold for IK1 (since Qvalues were not good)
log2.ratio.threshold <- 0.58
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Functions

# function to extract entry, splitting on semi-colon
extract.uniprotid.fun <- function(x){
  strsplit(x, split="\\|") %>% # split on vertical bar
    unlist() %>% # unlist
    nth(n=2) # take the second element of the vector
}

# function to extract last entry, splitting on semi-colon
extract.semicolon.last.fun <- function(x){
  strsplit(x, split=";") %>% # split on semicolon
    unlist() %>% # unlist
    last() # take the last element of the vector
}

# function to extract first entry, splitting on dash
remove.isoform.info.fun <- function(x){
  strsplit(x, split="-") %>% # split on dash
    unlist() %>% # unlist
    first() # take the first element of the vector
}
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# CLEAN THE DATA

df <- df.input %>%
  select(-contains("DEPRECATED")) %>% # drop DEPRECATED variables
  select(-Valid) %>% # drop Valid variable
  rename(Number.of.Ratios=X..of.Ratios, # rename these variables
         Percent.Change=X..Change,
         Number.of.Unique.Total.Peptides=X..Unique.Total.Peptides,
         Number.of.Unique.Total.EG.Id=X..Unique.Total.EG.Id) %>%
  select(UniProtIds, Genes, # move these variables to the beginning (we will widen these columns in the excel workbook)
         ProteinGroups, ProteinNames, ProteinDescriptions, 
         Comparison..group1.group2., 
         AVG.Log2.Ratio, Absolute.AVG.Log2.Ratio, Qvalue,
         GO.Biological.Process, GO.Molecular.Function, GO.Cellular.Component, 
         everything()) %>% # all other variables follow
  rename(UniProtIds.OG = UniProtIds) %>% # rename original UniProtIds, before creating the new UniProtIds column which we will use to join datasets
  mutate(UniProtIds = map_chr(UniProtIds.OG, extract.semicolon.last.fun)) %>% # extract last uniprotid 
  mutate(UniProtIds = map_chr(UniProtIds, remove.isoform.info.fun)) %>% # remove isoform suffix if present - retainingly only uniprotid
  select(UniProtIds, everything()) %>% # move UniProtIds to front
  rename(Group.OG = Group) %>% # rename original Group, before creating a new Group column which we will use to join datasets
  mutate(Group = str_replace_all(Group.OG, regex("(\\[.+\\])"), replacement = "")) %>% # remove peptide modifications (which are contained inside brackets "[]"), so that we can use this to join with protein pilot Group attribute
  relocate(Group, .after=Group.OG) %>%
  arrange(Qvalue) # arrange by Qvalue, ascending
#------------------------------------------------------------------------------------


# #-----------------------------------------------------------------------------------------------------
# # check which comparisons are present in the data
# unique(df$Comparison..group1.group2.) # run this to print to screen comparisons which are present in the data
# # optional
# # CORRECT SOME COMPARISONS
# # which comparisons do we want? 
# # which comparisons need to be flipped?
# comparisons.to.flip <- c("E_21_C / E_21_B") # add the comparisons you wish to flip and correct by hand
# # check that these comparisons are in the dataset (confirm no typos)
# comparisons.to.flip %in% unique(df$Comparison..group1.group2.) # should all be true
# 
# # correction function
# # some comparisons are inverted, and need to be flipped and corrected 
# flip.and.correct.fun <- function(df.x, comparisons){
#   
#   # split data that needs correction
#   df.y <- df.x %>%
#     filter(Comparison..group1.group2. %in% comparisons) %>% # filter for comparisons which need correction
#     mutate(Condition.Numerator.TEMP=Condition.Denominator) %>% # swap numerator & denominator conditions # temporary name
#     mutate(Condition.Denominator.TEMP=Condition.Numerator) %>% # swap numerator & denominator conditions # temporary name
#     mutate(AVG.Group.Quantity.Numerator.TEMP=AVG.Group.Quantity.Denominator) %>% # swap numerator & denominator values # temporary name
#     mutate(AVG.Group.Quantity.Denominator.TEMP=AVG.Group.Quantity.Numerator) %>% # swap numerator & denominator values # temporary name
#     select(-c(Condition.Numerator, Condition.Denominator, AVG.Group.Quantity.Numerator, AVG.Group.Quantity.Denominator)) %>% # drop these original variables, before renaming TEMP variables
#     rename(Condition.Numerator=Condition.Numerator.TEMP,  # rename (dropping TEMP suffix)
#            Condition.Denominator=Condition.Denominator.TEMP,
#            AVG.Group.Quantity.Numerator=AVG.Group.Quantity.Numerator.TEMP,
#            AVG.Group.Quantity.Denominator=AVG.Group.Quantity.Denominator.TEMP) %>%
#     mutate(AVG.Log2.Ratio=-AVG.Log2.Ratio) %>% # correct LogRatio (simply with a minus sign)
#     mutate(Ratio=1/Ratio) %>% # correct Ratio (simply reciprocating)
#     mutate(Percent.Change=(Ratio-1)*100) %>% # calculate new percent.change
#     mutate(Comparison..group1.group2.=paste(Condition.Numerator, Condition.Denominator, sep=" / "))
#   
#   # split data that does not need correction
#   df.x <- df.x %>%
#     filter(!Comparison..group1.group2. %in% comparisons) # filter for comparisons which are not being corrected
#   
#   # join newly corrected data back with the rest of the data
#   df.out <- df.x %>%
#     rbind(df.y) %>%
#     select(UniProtIds, Genes, ProteinGroups, ProteinNames, ProteinDescriptions, # rearrange variables
#            everything()) %>% # all other variables follow
#     arrange(Qvalue) # arrange by Qvalue, ascending
#   
#   # output
#   return(df.out)
# }
# 
# # apply correction function to dataframe
# df.corrected <- flip.and.correct.fun(df, comparisons.to.flip) # correct these comparisons
# 
# # print comparison ratios to screen - are these the ones we want?
# unique(df.corrected$Comparison..group1.group2.)
# 
# # join the old and the new and inspect the values - are they correct?
# df.inspect <- df.corrected %>%
#   select(UniProtIds, 
#          Genes, 
#          Comparison..group1.group2., 
#          AVG.Group.Quantity.Numerator, 
#          AVG.Group.Quantity.Denominator, 
#          AVG.Log2.Ratio, 
#          Absolute.AVG.Log2.Ratio, 
#          Ratio) %>%
#   full_join(df %>% select(UniProtIds, 
#                           Genes, 
#                           Comparison..group1.group2., 
#                           AVG.Group.Quantity.Numerator, 
#                           AVG.Group.Quantity.Denominator, 
#                           AVG.Log2.Ratio, 
#                           Absolute.AVG.Log2.Ratio, 
#                           Ratio), 
#             by=c("UniProtIds", "Genes", "Absolute.AVG.Log2.Ratio")) # also joining by Absolute.Log2.Ratio in order to get the correct match
# 
# # rename corrected dataframe for simplicity
# df <- df.corrected
# #-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# add peptide position to candidates data from df.peptide 

# simplify peptide information for joining
df.peptide.unique <- df.peptide %>%
  select(PG.ProteinAccessions, PEP.GroupingKey, PEP.PeptidePosition) %>% # want to join PeptidePosition to candidates file
  unique()

# join
df <- df %>%
  left_join(df.peptide.unique, by=c("UniProtIds"="PG.ProteinAccessions", "Group"="PEP.GroupingKey")) %>%
  relocate(Group, PEP.PeptidePosition, .after=Comparison..group1.group2.)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# add cleavage data from peptide data from protein pilot

# clean up this dataset 
df.pilot <- df.pilot %>%
  mutate(UniProtIds = map_chr(Accessions, extract.uniprotid.fun)) %>% # extract uniprotid 
  select(UniProtIds, everything())

# simplify cleavage information for joining
df.pilot.unique <- df.pilot %>%
  select(UniProtIds, Sequence, Cleavages) %>%
  mutate(Group = paste0("_", Sequence, "_")) %>% # add underscores to both ends of peptide sequence to match with Group from other dataset
  unique()

# check overlap before joining
sum(df$Group %in% df.pilot.unique$Group)
sum(df$UniProtIds %in% df.pilot.unique$UniProtIds)

# join - on the peptide sequence level (Birgit's suggestion ... seems to work better than joining on both protein and peptide)
df <- df %>%
  left_join(df.pilot.unique %>% select(-UniProtIds), by="Group") %>% # drop UniProtIds since it already exists in df and we aren't joining by it here
  relocate(Cleavages, .after=PEP.PeptidePosition)

# # join - by protein and peptide
# df <- df %>%
#   left_join(df.pilot.unique, by=c("UniProtIds", "Group")) %>%
#   relocate(Cleavages, .after=PEP.PeptidePosition)

# how many observations have cleavage information?
sum(!is.na(df$Cleavages))

# classify peptides as either semi-tryptic or tryptic based on cleave information
# tryptic if both ends were cleaved after R or K (have to investigate the N-terminus cleavage to confirm this)
# otherwise semi-tryptic (since Birgit did not allow for non-tryptic peptides in her workflow)
# any cases which did not join between datasets (spectronaut and protein pilot) will remain as NA (and will not be classified)
df <- df %>%
  mutate(Classification = ifelse(!startsWith(Cleavages, "cleaved"), "tryptic", # if no cleavage information then default to tryptic (tryptic peptides is the default for the search engine Birgit uses)
                                 ifelse(startsWith(Cleavages, "cleaved K-") & endsWith(Sequence, "K"), "tryptic", # cleaved at Lysine K on both sides
                                        ifelse(startsWith(Cleavages, "cleaved K-") & endsWith(Sequence, "R"), "tryptic", # cleaved K@N-term and R@C-term
                                               ifelse(startsWith(Cleavages, "cleaved R-") & endsWith(Sequence, "R"), "tryptic", # cleaved at Arginine R on both sides
                                                      ifelse(startsWith(Cleavages, "cleaved R-") & endsWith(Sequence, "K"), "tryptic", # cleaved R@N-term and K@C-term
                                                             "semi-tryptic"))))) # otherwise, classify as semi-tryptic since we could not confirm it was fully tryptic from the above cases
  ) %>%
  relocate(Classification, .after=Cleavages)

# check classification types
unique(df$Classification)

# check how many of each classification
df.class <- df %>%
  group_by(Classification) %>%
  summarise(COUNT = n()) %>%
  arrange(desc(COUNT))

# check out all the unique trypsin cleavages
df %>%
  filter(Classification=="tryptic") %>%
  pull(Cleavages) %>%
  unique()
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# for IK1, arrange by Pvalue (since Qvalue were not very good)
df <- df %>%
  relocate(Pvalue, .before=Qvalue) %>%
  arrange(Pvalue)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Make First Sheet
# This sheet contains all of the data (with the exception of always filtering out the one peptide wonders)

sheet.name <- "All proteins"

# ROUND
# Round the values for the following variables:
# columns_to_round <- c("AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio", "Ratio")
# df <- df %>% 
#   mutate_at(columns_to_round, round, digits = 2)


# HIDE
# before hiding columns (which needs to be done in the workbook itself, we must create the workbook using openxlsx):
# Create workbook:
wb <- createWorkbook("ECM")
addWorksheet(wb, sheet.name, gridLines = TRUE)
writeData(wb, sheet.name, df, rowNames = FALSE)
headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD", border="TopBottomLeftRight", borderColour = "black", wrapText = TRUE)
addStyle(wb, sheet.name, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
setColWidths(wb, sheet.name, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later


# Hide the following columns:
columns_to_hide <- c("AVG.Group.Quantity.Denominator", "AVG.Group.Quantity.Numerator", 
                     "Number.of.Ratios", "Percent.Change", "Standard.Deviation", 
                     "Condition.Numerator", "Condition.Denominator", 
                     "Degrees.Of.Freedom", "Standard.Error")
column.hide.indeces <- which(colnames(df) %in% columns_to_hide)

# hide columns
setColWidths(wb, sheet.name, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)


# WIDEN
# Widen the following columns for better legibility:
# "UniProtIds" # set to 15
# "Genes" # set to 10
# "ProteinGroups" # set to width of 18
# "ProteinNames" # set to width of 20
# "ProteinDescriptions" # set to width of 58
# "Comparison..group1.group2." # set to width of 20
# "AVG.Log2.Ratio" # set to width of 19
# "Qvalue" # set to width of 12
# "Absolute.AVG.Log2.Ratio" # set to width of 19
# "GO.Biological.Process" # set to width of 55
# "GO.Molecular.Function" # set to width of 55
# "GO.Cellular.Component" # set to width of 55
columns_to_widen <- c("UniProtIds", "Genes", 
                      "ProteinGroups", "ProteinNames", "ProteinDescriptions", 
                      "Comparison..group1.group2.", 
                      "AVG.Log2.Ratio", "Qvalue", "Absolute.AVG.Log2.Ratio",
                      "GO.Biological.Process", "GO.Molecular.Function", "GO.Cellular.Component")

column.widen.indeces <- which(colnames(df) %in% columns_to_widen)

# set widths for each column above (the order must match)
column.widths <- c(15, 10, 
                   18, 20, 58, 
                   20,
                   19, 12, 19, 
                   55, 55, 55)

# assign names to column.widths vector - to keep things straight
names(column.widths) <- columns_to_widen

# print them vector out to make sure widths are correct
print(column.widths)

# set widths in workbook
setColWidths(wb, sheet = sheet.name, cols = column.widen.indeces, widths = column.widths) # column widths


# CENTER
# center the following columns:
columns_to_center <- c("AVG.Log2.Ratio", "Qvalue", "Absolute.AVG.Log2.Ratio", "Ratio")
column.center.indeces <- which(colnames(df) %in% columns_to_center)

# align center:
# create style
align_center_style <- createStyle(fontName = "calibri", fontSize = 14, fontColour = NULL, numFmt = "GENERAL", border = NULL,
                                  borderColour = getOption("openxlsx.borderColour", "black"),
                                  borderStyle = getOption("openxlsx.borderStyle", "thin"),
                                  bgFill = NULL, fgFill = NULL, 
                                  halign = "center",  ## align text in center of the cell
                                  valign = NULL, textDecoration = NULL,
                                  wrapText = FALSE, textRotation = NULL, indent = NULL, locked = NULL, hidden = NULL)

# add style
addStyle(wb, sheet.name, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
# END SHEET 
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Make one separate sheet for each comparison
# which contains proteins with P < Pvalue.threshold

for(i in unique(df$Comparison..group1.group2.)){
  print(i)
  # sheet name is the comparison + descriptor
  # change comparison name slightly for writing to excel file
  sheet.name <- paste(gsub("\\/", "vs", i), "P-filter", sep=" ") # replace the slash division with 'vs'
  
  # significantly altered proteins for this specific comparison
  df.sig <- df %>%
    filter(Comparison..group1.group2.==i) %>% # filter for specific comparison for this sheet
    filter(Pvalue < Pvalue.threshold)
  
  addWorksheet(wb, sheet.name, gridLines = TRUE)
  writeData(wb, sheet.name, df.sig, rowNames = FALSE)
  
  addStyle(wb, sheet.name, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
  setColWidths(wb, sheet.name, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later
  
  # HIDE
  setColWidths(wb, sheet.name, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)
  
  # WIDEN
  setColWidths(wb, sheet = sheet.name, cols = column.widen.indeces, widths = column.widths) # column widths
  
  # CENTER
  addStyle(wb, sheet = sheet.name, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
  # END SHEET
}
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Make one separate sheet for each comparison
# which contains significantly altered proteins 

for(i in unique(df$Comparison..group1.group2.)){
  print(i)
  # sheet name is the comparison + descriptor
  # change comparison name slightly for writing to excel file
  sheet.name <- paste(gsub("\\/", "vs", i), "Significant", sep=" ") # replace the slash division with 'vs'
  
  # significantly altered proteins for this specific comparison
  df.sig <- df %>%
    filter(Comparison..group1.group2.==i) %>% # filter for specific comparison for this sheet
    filter(Pvalue < Pvalue.threshold) %>%
    filter(Absolute.AVG.Log2.Ratio > log2.ratio.threshold)
  
  addWorksheet(wb, sheet.name, gridLines = TRUE)
  writeData(wb, sheet.name, df.sig, rowNames = FALSE)
  
  addStyle(wb, sheet.name, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
  setColWidths(wb, sheet.name, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later
  
  # HIDE
  setColWidths(wb, sheet.name, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)
  
  # WIDEN
  setColWidths(wb, sheet = sheet.name, cols = column.widen.indeces, widths = column.widths) # column widths
  
  # CENTER
  addStyle(wb, sheet = sheet.name, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
  # END SHEET
}
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Make one separate sheet for each comparison
# which contains all of the proteins

for(i in unique(df$Comparison..group1.group2.)){
  print(i)
  # sheet name is the comparison + descriptor
  # change comparison name slightly for writing to excel file
  sheet.name <- paste(gsub("\\/", "vs", i), "All", sep=" ") # replace the slash division with 'vs'
  
  # all proteins for this specific comparison
  df.all <- df %>%
    filter(Comparison..group1.group2.==i) # filter for specific comparison for this sheet

  addWorksheet(wb, sheet.name, gridLines = TRUE)
  writeData(wb, sheet.name, df.all, rowNames = FALSE)
  
  addStyle(wb, sheet.name, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
  setColWidths(wb, sheet.name, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later
  
  # HIDE
  setColWidths(wb, sheet.name, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)
  
  # WIDEN
  setColWidths(wb, sheet = sheet.name, cols = column.widen.indeces, widths = column.widths) # column widths
  
  # CENTER
  addStyle(wb, sheet = sheet.name, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
  # END SHEET
}
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# RENAME THE SHEETS
# optional

# print names
names(wb)

# names(wb) <- c("All proteins", "Sheet 2", "Sheet 3", "etc...")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# create name for processed output file
output.filename <- input.filename %>%
  strsplit(split=".tsv") %>%
  unlist() %>%
  paste("_processed.xlsx", sep="")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# SAVE WORKBOOK
saveWorkbook(wb, file = output.filename, overwrite = TRUE)
#-----------------------------------------------------------------------------------------------------


# #-----------------------------------------------------------------------------------------------------
# # calculate number of proteins by condition and by direction
# df.grouped.all <- df %>%
#   # filter(Number.of.Unique.Total.Peptides > 1) %>% # filter for at least 2 unique peptides per protein - exclude one peptide wonders
#   select(Comparison..group1.group2., AVG.Log2.Ratio) %>%
#   dplyr::group_by(Comparison..group1.group2.) %>%
#   mutate(COUNT_ALL = n()) %>%
#   mutate(DIRECTION = ifelse(AVG.Log2.Ratio>0, "UP", "DOWN")) %>%
#   dplyr::group_by(Comparison..group1.group2., DIRECTION) %>%
#   mutate(COUNT_DIR = n()) %>%
#   select(-AVG.Log2.Ratio) %>%
#   unique()
# #-----------------------------------------------------------------------------------------------------
# 
# 
# #-----------------------------------------------------------------------------------------------------
# # calculate number of significantly changed proteins by condition and by direction
# df.grouped.sig <- df %>%
#   # filter(Number.of.Unique.Total.Peptides > 1) %>% # filter for at least 2 unique peptides per protein - exclude one peptide wonders
#   filter(Qvalue< Qvalue.threshold) %>%
#   filter(Absolute.AVG.Log2.Ratio> log2.ratio.threshold) %>%
#   select(Comparison..group1.group2., AVG.Log2.Ratio) %>%
#   mutate(DIRECTION = ifelse(AVG.Log2.Ratio>0, "UP", "DOWN")) %>%
#   dplyr::group_by(Comparison..group1.group2., DIRECTION) %>%
#   mutate(COUNT_SIG = n()) %>%
#   select(-AVG.Log2.Ratio) %>%
#   unique()
# #-----------------------------------------------------------------------------------------------------


# END
