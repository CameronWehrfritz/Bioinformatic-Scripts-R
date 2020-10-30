#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute
#April 22, 2020

# VERSION 1
# PROCESSING of sample L08 from Tlsty CRUK ECM project
# Conditions: A_contr, B_tumor
# OUTPUT: EXCEL File (with multiple sheets)


######################
#### Begin Program ###
######################

#set working directory
#setwd("/Volumes/GibsonLab/users/Cameron/2020_0422_L08_CRUK/R_workspace")
setwd("//bigrock/GibsonLab/users/Cameron/2020_0422_L08_CRUK/R_workspace") # VPN


#-----------------------------------------------------------------------------------------------------
# PACKAGES #
packages = c("plyr", "dplyr", "openxlsx", "readxl", "writexl", "reshape2", "gplots")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# LOAD DATA #
#############
df <- read.csv("//bigrock/GibsonLab/users/Cameron/2020_0422_L08_CRUK/Tables/Input/Candidates_PG17_L08_0328_2020_v1.tsv", sep="\t", stringsAsFactors = FALSE) # VPN
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
# NOTES FOR DATA:
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
# CLEAN THE DATA #
# For all sheets:

# keep only ratios with Normal condition ('A_N') in the denominator
df <- subset(df, grepl("*A_contr", Condition.Denominator) )


# Delete the following columns:
# anything with "DEPRECATED"
df <- df[ , !grepl("*DEPRECATED*", colnames(df))]

# Order ascending based on Qvalue
df <- arrange(df, Qvalue)


#########
# ROUND
# Round the values in the following columns to 2 decimal places:
columns_to_round <- c("AVG.Log2.Ratio", "Absolute.AVG.Log2.Ratio", "Ratio")
df <- df %>% 
  mutate_at(columns_to_round, round, digits = 2)
#########


#########
# HIDE
# before hiding columns (which needs to be done in the workbook itself, we must create the workbook using openxlsx):
# Create workbook:
wb <- createWorkbook("ECM")
addWorksheet(wb, "Sheet 1", gridLines = TRUE)
writeData(wb, sheet = 1, df, rowNames = FALSE)
headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD", border="TopBottomLeftRight", borderColour = "black") #"#4F81BD")
addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
setColWidths(wb, sheet = 1, cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later


# Hide the following columns:
# "AVG.Group.Quantity.Denominator"
# "AVG.Group.Quantity.Numerator"
# "X..Change" 
# "Standard.Deviation"
columns_to_hide <- c("AVG.Group.Quantity.Denominator", "AVG.Group.Quantity.Numerator", "X..of.Ratios", "X..Change", "Standard.Deviation")
column.hide.indeces <- which(colnames(df) %in% columns_to_hide)

# hide columns:
setColWidths(wb, "Sheet 1", cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)
#########


#########
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
#########


#########
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
addStyle(wb, sheet = "Sheet 1", style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
#########
# END SHEET 1 #
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
# BEGIN SHEET 2 #
# B_tumor/A_contr
# Q < 0.01

addWorksheet(wb, "Sheet 2", gridLines = TRUE)
writeData(wb, sheet = "Sheet 2", df[ df$Qvalue<0.01 & grepl("*B_tumor", df$Condition.Numerator), ], rowNames = FALSE) # take only B_cancer data

addStyle(wb, sheet = "Sheet 2", headerStyle, rows = 1, cols = 1:length(colnames(df)), gridExpand = TRUE)
setColWidths(wb, sheet = "Sheet 2", cols=1:length(colnames(df)), widths = 10) # set all column widths to 10; more tailoring to come later

########
# HIDE
setColWidths(wb, "Sheet 2", cols = column.hide.indeces, widths = 8.43, hidden = TRUE, ignoreMergedCells = FALSE)
#########


#########
# WIDEN
setColWidths(wb, sheet = "Sheet 2", cols = column.widen.indeces, widths = c(34, 24, 28, 18, 28, 16, 24, 32, 14, 50, 50, 50))
#########


##########
# CENTER
# add style
addStyle(wb, sheet = "Sheet 2", style = align_center_style, rows = 1:dim(df)[1], cols = column.center.indeces, gridExpand = TRUE, stack = TRUE)
##########
# END SHEET 2 #
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# RENAME THE SHEETS #
names(wb) <- c("All_vs_Acontr", "Btumor_vs_Acontr Q_0.01")
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
## SAVE WORKBOOK ##
saveWorkbook(wb, file = "test_wb.xlsx", overwrite = TRUE)
#-----------------------------------------------------------------------------------------------------


# DONE CLEANING #
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------


#######
# END #
#######