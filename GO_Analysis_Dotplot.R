#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#April 16, 2021

# Generates Dotplot of Gene Set Enrichment Analysis (GSEA) Gene Ontology (GO) terms
# Input data is from a GSEA from Consesus Pathway Database

### Begin Script ###

#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2021_0222_EG5/Pathway_analysis/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0222_EG5/Pathway_analysis/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("clusterProfiler", "AnnotationHub", "org.Hs.eg.db", "dplyr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# load data
data.input <- read.xlsx("//bigrock/GibsonLab/users/Cameron/2021_0222_EG5/Pathway_analysis/CPDB_output/EG5_Enriched_gene_ontology_based_sets_Downregulated_T5KO_DC_vs_T5K0_ND_ORA_results.xlsx", sheet=2)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# clean data

data <- data.input %>%
  filter(`q-value`<0.05) %>%
  arrange(desc(Gene.Ratio)) %>%
  filter(term_category=="b") %>%
  dplyr::rename(q.value=`q-value`) %>% # rename
  dplyr::rename(pathway=term_name) # rename to pathway
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# dotplot of Gene Ontology pathways

S1 <- data %>% 
  ggplot(aes(x=Gene.Ratio , y=factor(pathway, levels=rev(pathway)), size=count, color=q.value)) + geom_point(alpha = 0.8) + 
  theme_classic()

S1 
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# dotplot of Gene Ontology pathways
# better visualization

# change color scheme
S2 <- S1 + scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(min(data$q.value), 
                                                                                          max(data$q.value)))
# clean up plot labels

S2 + scale_y_discrete(name="") +
  scale_size(range = c(2, 8)) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=12), 
        axis.text.y = element_text(size=12))
#------------------------------------------------------------------------------------


# END
