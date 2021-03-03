#Written by Cameron Wehrfritz
#March 3, 2021

# Processing Mass Spec Candidates File

# OUTPUT: PROCESSED DATA TABLES

####################
### Begin Script ###
####################


#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/0303_2021_AC1/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/R_workspace") # PC
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

df.input <- read.csv("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/Input/AC1_0302_2021_v3_DDA_AC1_libr.tsv", sep="\t", stringsAsFactors = F) #PC
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


#------------------------------------------------------------------------------------
# Make Sheet #1

# first set sheet name
sheet.number <- "sheet 1"

# ROUND
# Round the values in the following columns to 2 decimal places:
columns_to_round <- c("AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio", "Ratio")
df <- df %>% 
  mutate_at(columns_to_round, round, digits = 2)


# HIDE
# before hiding columns (which needs to be done in the workbook itself, we must create the workbook using openxlsx):
# Create workbook:
wb <- createWorkbook("ECM")
addWorksheet(wb, sheet.number, gridLines = TRUE)
writeData(wb, sheet.number, df, rowNames = FALSE)
headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD", border="TopBottomLeftRight", borderColour = "black") #"#4F81BD")
addStyle(wb, sheet.number, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
setColWidths(wb, sheet.number, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later


# Hide the following columns:
# "AVG.Group.Quantity.Denominator"
# "AVG.Group.Quantity.Numerator"
# "X..Change" 
# "Standard.Deviation"
columns_to_hide <- c("AVG.Group.Quantity.Denominator", "AVG.Group.Quantity.Numerator", "Number.of.Ratios", "Percent.Change", "Standard.Deviation")
column.hide.indeces <- which(colnames(df) %in% columns_to_hide)

# hide columns:
setColWidths(wb, sheet.number, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)


# WIDEN
# Widen the following columns for better legibility:
# "Comparison..group1.group2." # set to width of 34
# "Condition.Numerator" # set to width of 24
# "Condition.Denominator" # set to width of 28
# "AVG.Log2.Ratio" # set to width of 18
# "Absolute.AVG.Log2.Ratio" # set to width of 24
# "ProteinGroups" # set to width of 16
# "ProteinNames" # set to width of 24
# "ProteinDescriptions" # set to width of 32
# "UniProtIds" # set to 14 
# "GO.Biological.Process" # set to width of 50
# "GO.Molecular.Function" # set to width of 50
# "GO.Cellular.Component" # set to width of 50
columns_to_widen <- c("Comparison..group1.group2.", "Condition.Numerator", "Condition.Denominator", "AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio", 
                      "ProteinGroups", "ProteinNames", "ProteinDescriptions", "UniProtIds", "GO.Biological.Process", "GO.Molecular.Function", "GO.Cellular.Component")
column.widen.indeces <- which(colnames(df) %in% columns_to_widen)
setColWidths(wb, sheet = 1, cols = column.widen.indeces, widths = c(34, 24, 28, 18, 28, 16, 24, 32, 14, 50, 50, 50)) # column widths



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
addStyle(wb, sheet.number, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
# END SHEET 1 #
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
# Make Sheet #2

# significantly altered proteins with at least 2 unique peptides
df.sig <- df %>%
  filter(Qvalue<=0.01) %>%
  filter(Absolute.AVG.Log2.Ratio>0.58) %>%
  filter(Number.of.Unique.Total.Peptides>1) # must have at least 2 unique peptides per protein

# first set sheet name
sheet.number <- "sheet 2"

addWorksheet(wb, sheet.number, gridLines = TRUE)
writeData(wb, sheet.number, df.sig, rowNames = FALSE) # significantly altered proteins

addStyle(wb, sheet.number, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
setColWidths(wb, sheet.number, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later

# HIDE
setColWidths(wb, sheet.number, cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)

# WIDEN
setColWidths(wb, sheet = sheet.number, cols = column.widen.indeces, widths = c(34, 24, 28, 18, 28, 16, 24, 32, 14, 50, 50, 50))

# CENTER
# add style
addStyle(wb, sheet = sheet.number, style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)

# END SHEET 2 #
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# RENAME THE SHEETS #
names(wb) <- c("All Proteins", "Significant Proteins")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
## SAVE WORKBOOK ##
saveWorkbook(wb, file = "Table_processed.xlsx", overwrite = TRUE)
#-----------------------------------------------------------------------------------------------------
