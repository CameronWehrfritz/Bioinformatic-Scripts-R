#Written by Cameron Wehrfritz
#Schilling Lab, Buck Institute
#March 29, 2021

# Process peptide and protein level mass spectrometry DDA data

# Output: 
# i. peptide level summary
# ii. protein level summary

#### Begin Program ###

#-----------------------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/0303_2021_AC1/R_workspace") # mac
setwd("//bigrock/GibsonLab/users/Cameron/0303_2021_AC1/R_workspace") # pc
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# PACKAGES #
packages = c("dplyr", "openxlsx", "readxl", "writexl", "reshape2", "ggplot2", "VennDiagram",
             "gplots", "scales", "forcats", "tidyr")

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# LOAD DATA

df.input <- read.delim("//bigrock/GibsonLab/users/birgit/Verdin_Covarrubias_2020/AC1_2021/Cameron_Birgit/210224_0006_AC1_MASTER_cm_ALL_PeptideSummary.txt")

# check confidence level range
range(df.input$Conf)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# filter by confidence level
df <- df.input %>%
  filter(Conf>98.9) 
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# calculate number of peptides per protein
df <- df %>%
  group_by(Accessions) %>%
  mutate(Number.of.Peptides=Sequence %>% unique() %>% length()) %>%
  ungroup()
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# make protein level table - one entry per protein
# df.protein.level <- df %>%
#   select(N, Unused, Total, `%Cov(95)`, Accession=Accessions, Name=Names, Number.of.Peptides) %>%
#   unique()

# for .txt files
df.protein.level <- df %>%
  select(N, Unused, Total, X.Cov.95., Accession=Accessions, Name=Names, Number.of.Peptides) %>%
  unique()

# check for duplicates
sum(duplicated(df.protein.level$Accession))
sum(duplicated(df.protein.level$N))

# check range of number.of.proteins
range(df.protein.level$Number.of.Peptides)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# make table with 2+ peptides per protein
df.protein.level.multi <- df.protein.level %>%
  filter(Number.of.Peptides>1)

# check range of number.of.proteins
range(df.protein.level.multi$Number.of.Peptides)
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# histogram number.of.peptide distribution

# all proteins
df.protein.level %>%
  ggplot(aes(x=Number.of.Peptides)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.5) +
  xlab("Number of Peptides") +
  ylab("Number of Proteins") +
  ggtitle("Secretome Macrophage Senescence Project AC1 DDA Results") +
  # stat_bin(binwidth = 1, geom="text", aes(label=..count..), position = position_stack(vjust = 0.5)) + # add count as annotation
  theme_minimal()

# exclude one peptide wonders
df.protein.level %>%
  filter(Number.of.Peptides>1) %>% # exclude one peptide wonders
  ggplot(aes(x=Number.of.Peptides)) +
  geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.5) +
  xlab("Number of Peptides") +
  ylab("Number of Proteins") +
  ggtitle("Secretome Macrophage Senescence Project AC1 DDA Results") +
  # stat_bin(binwidth = 1, geom="text", aes(label=..count..), position = position_stack(vjust = 0.5)) + # add count as annotation
  theme_minimal()

# df.protein.level %>%
#   ggplot(aes(x=Number.of.Peptides)) +
#   geom_histogram(binwidth = 1, fill="#69b3a2", color="#e9ecef", alpha=0.5) +
#   xlab("Number of Peptides") +
#   ylab("Number of Proteins") +
#   ggtitle("Secretome Macrophage Senescence Project AC1 DDA Results") +
#   xlim(0, 50) +
#   # stat_bin(binwidth = 1, geom="text", aes(label=..count..), position = position_stack(vjust = 0.5)) + # add count as annotation
#   theme_minimal()
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# write out excel files
write_xlsx(df, "Table_peptide_level.xlsx")
write_xlsx(df.protein.level, "Table_protein_level_all.xlsx")
write_xlsx(df.protein.level.multi, "Table_protein_level_atleast2peptides.xlsx")
#-----------------------------------------------------------------------------------------------------


# END
