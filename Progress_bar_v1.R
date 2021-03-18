# Cameron Wehrfritz
# March 18, 2021

# Progress Bar
# Progress Bars under MS Windows

#------------------------------------------------------------------------------------
# example with some toy data
mydata=matrix(rnorm(6000*300),ncol = 300)
result=as.data.frame(matrix(nrow = 6000,ncol = 2))
progression <- winProgressBar(title = "Progress bar", min = 0, max = 6000 , width = 300)
for (i in 1:6000) {
  result[i,1]=mean(mydata[i,]) # calculate mean
  result[i,2]=median(mydata[i,]) # calculate median
  setWinProgressBar(progression, i, title=paste(round(i/6000)*100,"% done"))
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------
# example of a progress bar with multiple for loops (incrementing through indeces i,j,k)
# use in TurnoveR 

# set progress bar
progression <- winProgressBar(title = "Precursor Pool Progress bar", min = 0, max = nrow(df.areas.charge), width = 300)

# update progress bar
setWinProgressBar(progression, i*j*k, title=paste(round(i*j*k/nrow(df.areas.charge), digits=4)*100,"% done - Precursor Pool Progress bar"))

# close progress bar
close(progression)
#------------------------------------------------------------------------------------
