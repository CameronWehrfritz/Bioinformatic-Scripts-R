#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#May 10, 2021

# Generates Dotplot of Gene Set Enrichment Analysis (GSEA) Gene Ontology (GO) terms
# Input data is from a GSEA from Consesus Pathway Database

### Begin Script ###

#------------------------------------------------------------------------------------
# set working directory

# setwd("/Volumes/GibsonLab/users/Cameron/2021_0507_EG_PTM/2021_0507_EG7B/Pathway_analysis/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0507_EG_PTM/2021_0507_EG7B/Pathway_analysis/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# packages

packages = c("purrr", "ggplot2", "dplyr", "openxlsx")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# load data

# ID vs CD - Upregulated
df.input <- read.xlsx("//bigrock/GibsonLab/users/Cameron/2021_0507_EG_PTM/2021_0507_EG7B/Pathway_analysis/CPDB_output/EG7B_2021_0422_Lysate_directDIA_v3_MinusEG7B_12_Processed_2021_0423_v1_ID_vs_CD_Upregulated_Results.xlsx")
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# functions

# function to count elements separated by semicolon
length.semicolon.fun <- function(x){
  strsplit(x, split=";") %>% # split on semicolon
    unlist() %>% # unlist
    length() # calculate length of vector
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# clean data

df <- df.input %>%
  dplyr::rename(q.value=`q-value`) %>% # rename
  dplyr::rename(pathway=term_name) # rename to pathway
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# calculate data attributes
# count, gene.ratio (not corrected) and gene.ratio (corrected)

df <- df %>%
  mutate(count = map_dbl(members_input_overlap_geneids, length.semicolon.fun)) %>% # count number of geneids per GO term (each row)
  mutate(Gene.Ratio.noncorrected = count/size) %>% # gene ratio noncorrected is count divided by size
  mutate(Gene.Ratio = count/effective_size) %>% # gene ratio is count divided by effective size
  arrange(desc(Gene.Ratio)) # arrange descending by Gene.Ratio
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# filter data for plotting

df.sig <- df %>%
  filter(q.value<0.05) %>%
  filter(term_category=="b") %>%
  arrange(desc(Gene.Ratio)) # arrange descending by Gene.Ratio
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# dotplot of Gene Ontology pathways

S1 <- df.sig %>% 
  ggplot(aes(x=Gene.Ratio , y=factor(pathway, levels=rev(pathway)), size=count, color=q.value)) + geom_point(alpha = 0.8) + 
  theme_classic()

S1 # print plot to screen
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# dotplot of Gene Ontology pathways
# better visualization

# change color scheme
S2 <- S1 + 
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(df$q.value), max(df$q.value))) +
  # clean up plot labels
  scale_y_discrete(name="") +
  scale_size(range = c(2, 8)) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12), 
        axis.text.y = element_text(size=12))

S2 # print plot to screen
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# top N GO pathways
# N=40

S.top <- df.sig %>% 
  top_n(n=40) %>% # top 40 pathways is a good number to focus on, if there are more than 40 pathways the plot is too cluttered
  ggplot(aes(x=Gene.Ratio , y=factor(pathway, levels=rev(pathway)), size=count, color=q.value)) + geom_point(alpha = 0.8) + 
  theme_classic() +
  # change color scheme
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(df$q.value), max(df$q.value))) +
  # clean up plot labels
  scale_y_discrete(name="") +
  scale_size(range = c(2, 8)) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12), 
        axis.text.y = element_text(size=12))

S.top # print plot to screen
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# write out results in excel table
write.xlsx(list(df), "Table_Pathway_Analysis_output.xlsx")
# write.xlsx(list(df.id.vs.ic, df.cd.vs.cc, df.id.vs.cd, df.ic.vs.cc), "Table_output.xlsx")
#------------------------------------------------------------------------------------


# END
