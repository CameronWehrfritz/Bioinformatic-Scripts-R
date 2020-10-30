#written by Cameron Wehrfritz 
#Schilling Lab, Buck Institute
#August 1, 2019

# VERSION 1
# VOLCANO PLOTTER

#### Begin Program ###

#set working directory
setwd("/Volumes/GibsonLab/users/Cameron/2019_0715_TT2_CETSA")


############
# PACKAGES #
############
install.packages('calibrate')
suppressPackageStartupMessages(library(calibrate))


#############
# LOAD DATA #
#############
data <- read.csv("Table_rapamycin_df_ouput_0802_2019_v1.csv", stringsAsFactors = F)


########
# PLOT #
########
# plotting log2(ratio) and log10(Pvalue) from df "data"

graphics.off()

# initiate pdf
pdf(file="volcano_plot_test.pdf")

# basic plot
with(data, plot(log2(Ratio), -log10(Pvalue), pch=20, main="Volcano Plot"))

# add red colored points for significant proteins (P<0.05)
with(subset(data, Pvalue<0.05), points(log2(Ratio), -log10(Pvalue), pch=20, col="red"))

# label points with textxy function from calibrate package
# nB m=(x,y) is the origin, and should be chosen specifically to be the location of a good location/a certain protein's coordinates
#so that the labels are legible (as they fan out away from the origin)
with(subset(data, Pvalue<0.05), textxy(log2(Ratio), -log10(Pvalue), labs=Genes, m= c(log2(2.0841232),-log10(0.033687459)), cex=0.75)) 
# NAPRT m= c(log2(2.0841232),-log10(0.033687459)) # fairly legible labels

# add horizontal line for Pvalue=0.05
abline(h=-log10(0.05))

#add timestamp
mtext(date(), side=3, line=0, adj=0) # side (1=bottom, 2=left, 3=top, 4=right)

graphics.off()

#######
# END #
#######
