#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#April 15, 20201

# This script produces volcano plots


### Begin Script ###


#------------------------------------------------------------------------------------
#set working directory
# setwd("/Volumes/GibsonLab/users/Cameron/2021_0415_EG7/R_workspace") # MAC
setwd("//bigrock/GibsonLab/users/Cameron/2021_0415_EG7/R_workspace") # PC
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# PACKAGES #

packages = c("calibrate", "stringr", "VennDiagram", "hablar", "readxl", "writexl", "openxlsx",
             "purrr", "dplyr", "tidyr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# load data

df.input <- read_xlsx("//bigrock/GibsonLab/users/Cameron/2021_0415_EG7/Input/Candidates_EG7_2021_0412_directDIA_2MC_QSp_v2_Processed_2021_0415_v1_CW_Input.xlsx", sheet=2)
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# set significance thresholds
Qvalue.threshold <- 0.05
log2.ratio.threshold <- 0.58
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# volcano plot

# grab data to be plotted
data <- df.input %>% 
  filter(Comparison=="CD / CC")
  
# plot
with(data, plot(AVG.Log2.Ratio, -log10(Qvalue), pch=20, main=paste0("Comparison: ", unique(data$Comparison)), xlab = "Log2 Fold Change", col="gray"))

# add blue colored points for significantly downregulated proteins 
with(subset(data, Qvalue<Qvalue.threshold & AVG.Log2.Ratio < -log2.ratio.threshold ), points(AVG.Log2.Ratio, -log10(Qvalue), pch=20, col="blue"))

# add red colored points for significantly upregulated proteins 
with(subset(data, Qvalue<Qvalue.threshold & AVG.Log2.Ratio > log2.ratio.threshold ), points(AVG.Log2.Ratio, -log10(Qvalue), pch=20, col="red"))

# label significantly dysregulated proteins with textxy function from {calibrate}
# nB m=(x,y) is the origin, and should be chosen specifically to be the location of a good location/a certain protein's coordinates
# so that the labels are legible (as they fan out away from the origin)
# with(subset(data, Qvalue<Qvalue.threshold & Absolute.AVG.Log2.Ratio > log2.ratio.threshold), textxy(AVG.Log2.Ratio, -log10(Qvalue), labs=Genes, cex=0.75))

# add significance threshold lines
abline(h=-log10(Qvalue.threshold), lty=2)
abline(v=log2.ratio.threshold, lty=2)
abline(v=-log2.ratio.threshold, lty=2)

# how many observations in each group? SIGNIFICANT UP, SIGNIFICANT DOWN or NOT SIGNIFICANT
data.summary <- data %>%
  mutate(GROUP = ifelse( Qvalue<Qvalue.threshold & AVG.Log2.Ratio < -log2.ratio.threshold, "SIGNIFICANT DOWN",
                         ifelse(Qvalue<Qvalue.threshold & AVG.Log2.Ratio > log2.ratio.threshold, "SIGNIFICANT UP",
                                "NOT SIGNIFICANT"))) %>%
  group_by(GROUP) %>%
  summarise(COUNT = n()) %>% # calculate number in each group
  mutate(COUNT_PERCENTAGE = COUNT/sum(COUNT)*100 %>% round(digits=1)) # calculate percentage in each group
#------------------------------------------------------------------------------------


# END
