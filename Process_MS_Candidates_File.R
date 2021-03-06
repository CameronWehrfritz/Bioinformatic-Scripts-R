#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute for Research on Aging
#Novato, California, USA
#March 11, 2021

# DESCRIPTION: This script processes Mass Spectrometry Candidates Files
# i. cleans data
# ii. corrects certain comparisons (if they are in the inverted order - this must be specified by hand, see 'CORRECT SOME COMPARISONS' chunk below)
# iii. write processed data to excel workbook with multiple sheets

# OUTPUT: PROCESSED EXCEL FILE

### Begin Script ###

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2021_0310_E20_CRUK/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0310_E20_CRUK/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #
packages = c("tidyr", "dplyr", "writexl", "openxlsx")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# LOAD DATA #

df.input <- read.csv("//bigrock/GibsonLab/users/Cameron/2021_0310_E20_CRUK/Input/Candidates_PG20_E20_0305_2021_v2.tsv", sep="\t", stringsAsFactors = F) #PC

# define file name for output file
input.filename <- "Candidates_PG20_E20_0305_2021_v2.tsv"
#------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# CLEAN THE DATA

df <- df.input %>%
  select(-contains("DEPRECATED")) %>% # drop DEPRECATED variables
  select(-Valid) %>% # drop Valid variable
  rename(Number.of.Ratios=X..of.Ratios, # rename these variables
         Percent.Change=X..Change,
         Number.of.Unique.Total.Peptides=X..Unique.Total.Peptides,
         Number.of.Unique.Total.EG.Id=X..Unique.Total.EG.Id) %>%
  select(UniProtIds, Genes, ProteinGroups, ProteinNames, ProteinDescriptions, # move these variables to the beginning
         everything()) %>% # all other variables follow
  arrange(Qvalue) # arrange by Qvalue, ascending
#------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# CORRECT SOME COMPARISONS
# which comparisons do we want? 
unique(df$Comparison..group1.group2.) # run this to print to screen comparisons which are present in the data
# which comparisons need to be flipped?
comparisons.to.flip <- c() # currently empty # add the comparisons you wish to flip and correct by hand
                         

# correction function
# some comparisons are inverted, and need to be flipped and corrected 
flip.and.correct.fun <- function(df.x, comparisons){
  
  # split data that needs correction
  df.y <- df.x %>%
    filter(Comparison..group1.group2. %in% comparisons) %>% # filter for comparisons which need correction
    mutate(Condition.Numerator.TEMP=Condition.Denominator) %>% # swap numerator & denominator conditions # temporary name
    mutate(Condition.Denominator.TEMP=Condition.Numerator) %>% # swap numerator & denominator conditions # temporary name
    mutate(AVG.Group.Quantity.Numerator.TEMP=AVG.Group.Quantity.Denominator) %>% # swap numerator & denominator values # temporary name
    mutate(AVG.Group.Quantity.Denominator.TEMP=AVG.Group.Quantity.Numerator) %>% # swap numerator & denominator values # temporary name
    select(-c(Condition.Numerator, Condition.Denominator, AVG.Group.Quantity.Numerator, AVG.Group.Quantity.Denominator)) %>% # drop these original variables, before renaming TEMP variables
    rename(Condition.Numerator=Condition.Numerator.TEMP,  # rename (dropping TEMP suffix)
           Condition.Denominator=Condition.Denominator.TEMP,
           AVG.Group.Quantity.Numerator=AVG.Group.Quantity.Numerator.TEMP,
           AVG.Group.Quantity.Denominator=AVG.Group.Quantity.Denominator.TEMP) %>%
    mutate(AVG.Log2.Ratio=-AVG.Log2.Ratio) %>% # correct LogRatio (simply with a minus sign)
    mutate(Ratio=1/Ratio) %>% # correct Ratio (simply reciprocating)
    mutate(Percent.Change=(Ratio-1)*100) %>% # calculate new percent.change
    mutate(Comparison..group1.group2.=paste(Condition.Numerator, Condition.Denominator, sep=" / "))
  
  # split data that does not need correction
  df.x <- df.x %>%
    filter(!Comparison..group1.group2. %in% comparisons) # filter for comparisons which are not being corrected
  
  # join newly corrected data back with the rest of the data
  df.out <- df.x %>%
    rbind(df.y) %>%
    select(UniProtIds, Genes, ProteinGroups, ProteinNames, ProteinDescriptions, # rearrange variables
           everything()) %>% # all other variables follow
    arrange(Qvalue) # arrange by Qvalue, ascending
  
  # output
  return(df.out)
}


# apply correction function to dataframe
df.corrected <- flip.and.correct.fun(df, comparisons.to.flip) # correct these comparisons - vector created above by hand, at the beginning of this chunk

# print comparison ratios to screen - are these the ones we want?
unique(df.corrected$Comparison..group1.group2.)

# join the old and the new and inspect the values - are they correct?
df.inspect <- df.corrected %>%
  select(UniProtIds, 
         Genes, 
         Comparison..group1.group2., 
         AVG.Group.Quantity.Numerator, 
         AVG.Group.Quantity.Denominator, 
         AVG.Log2.Ratio, 
         Absolute.AVG.Log2.Ratio, 
         Ratio) %>%
  full_join(df %>% select(UniProtIds, 
                          Genes, 
                          Comparison..group1.group2., 
                          AVG.Group.Quantity.Numerator, 
                          AVG.Group.Quantity.Denominator, 
                          AVG.Log2.Ratio, 
                          Absolute.AVG.Log2.Ratio, 
                          Ratio), 
            by=c("UniProtIds", "Genes", "Absolute.AVG.Log2.Ratio"))

# rename corrected dataframe for simplicity
df <- df.corrected
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Make First Sheet
# This sheet contains all of the data (no filters)

# # first set sheet name
# x <- 0 # set x equal to 0 for the very first sheet
# sheet.name <- paste("sheet", x+1)
# x <- x+1 # increment x

sheet.name <- "All proteins"

# ROUND
# Round the values in the following columns to 2 decimal places:
columns_to_round <- c("AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio", "Ratio")
df <- df %>% 
  mutate_at(columns_to_round, round, digits = 2)


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
# "AVG.Group.Quantity.Denominator"
# "AVG.Group.Quantity.Numerator"
# "X..Change" 
# "Standard.Deviation"
columns_to_hide <- c("AVG.Group.Quantity.Denominator", "AVG.Group.Quantity.Numerator", "Number.of.Ratios", "Percent.Change", "Standard.Deviation")
column.hide.indeces <- which(colnames(df) %in% columns_to_hide)

# hide columns:
setColWidths(wb, sheet.name, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)


# WIDEN
# Widen the following columns for better legibility:
# "UniProtIds" # set to 15
# "ProteinGroups" # set to width of 18
# "ProteinNames" # set to width of 20
# "ProteinDescriptions" # set to width of 40
# "Comparison..group1.group2." # set to width of 20
# "Condition.Numerator" # set to width of 15
# "Condition.Denominator" # set to width of 15
# "AVG.Log2.Ratio" # set to width of 19
# "Absolute.AVG.Log2.Ratio" # set to width of 19
# "GO.Biological.Process" # set to width of 55
# "GO.Molecular.Function" # set to width of 55
# "GO.Cellular.Component" # set to width of 55
columns_to_widen <- c("UniProtIds", "ProteinGroups", "ProteinNames", "ProteinDescriptions", 
                      "Comparison..group1.group2.", "Condition.Numerator", "Condition.Denominator", 
                      "AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio",
                      "GO.Biological.Process", "GO.Molecular.Function", "GO.Cellular.Component")

column.widen.indeces <- which(colnames(df) %in% columns_to_widen)

# set widths for each column above (the order must match)
column.widths <- c(15, 18, 20, 40, 
                   20, 15, 15, 
                   19, 19, 
                   55, 55, 55)

# assign names to column.widths vector - to keep things straight
names(column.widths) <- columns_to_widen

setColWidths(wb, sheet = 1, cols = column.widen.indeces, widths = column.widths) # column widths



# CENTER
# Center text:
# columns to center:
# "AVG.Log2.Ratio"
# "Absolute.AVG.Log2.Ratio"
# "Ratio"
columns_to_center <- c("AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio", "Ratio")
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
# Make Second Sheet
# This sheet contains all proteins with at least 2 unique peptides

# # first set sheet name
# sheet.name <- paste("sheet", x+1)
# x <- x+1 # increment x

sheet.name <- "All proteins 2+ peptides"


# all proteins with at least 2 unique peptides
df.sig <- df %>%
  filter(Number.of.Unique.Total.Peptides>1) # must have at least 2 unique peptides per protein

addWorksheet(wb, sheet.name, gridLines = TRUE)
writeData(wb, sheet.name, df.sig, rowNames = FALSE) # significantly altered proteins

addStyle(wb, sheet.name, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
setColWidths(wb, sheet.name, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later

# HIDE
setColWidths(wb, sheet.name, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)

# WIDEN
setColWidths(wb, sheet = sheet.name, cols = column.widen.indeces, widths = column.widths) # column widths

# CENTER
# add style
addStyle(wb, sheet = sheet.name, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)

# END SHEET
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# Make one separate sheet for each comparison
for(i in unique(df$Comparison..group1.group2.)){
  print(i)
  # sheet name is the comparison 
  # tweaked comparison name slightly for writing to excel file
  sheet.name <- gsub("\\/", "vs", i)
    
  # # first set sheet name
  # sheet.name <- paste("sheet", x+1)
  # x <- x+1 # increment x

  # significantly altered proteins with at least 2 unique peptides for this specific comparison
  df.sig <- df %>%
    filter(Qvalue<=0.05) %>%
    filter(Absolute.AVG.Log2.Ratio>0.58) %>%
    filter(Number.of.Unique.Total.Peptides>1) %>% # must have at least 2 unique peptides per protein
    filter(Comparison..group1.group2.==i) # i is the specific comparison for this sheet
  
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
# RENAME THE SHEETS
# optional

# names(wb) <- c("All Proteins", "Significant Proteins")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# create name for processed file
output.filename <- input.filename %>%
  strsplit(split=".tsv") %>%
  unlist() %>%
  paste("_processed.xlsx", sep="")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# SAVE WORKBOOK
saveWorkbook(wb, file = output.filename, overwrite = TRUE)
#-----------------------------------------------------------------------------------------------------
