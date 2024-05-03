#sources files and finds best-fit effective resistance for example site
#Right-click on the file name and used "copy path" to obtain path string
pathstr <- "~/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Research/Ellen/Model Files/" #need to input machine-specific path string, or update

#source necessary files
source(paste(pathstr,"recordedFieldData.R",sep=""))
source(paste(pathstr,"createriverRandomSpacing.R",sep=""))
source(paste(pathstr,"createriverUniformSpacing.R",sep=""))
source(paste(pathstr,"makelookuparrayv6.R",sep=""))
source(paste(pathstr,"findViterate.R",sep=""))
source(paste(pathstr,"findBWvolv2.R",sep=""))
source(paste(pathstr,"findoutflowv7.R",sep=""))
source(paste(pathstr,"getinflow.R",sep=""))
source(paste(pathstr,"odes3.R",sep=""))
source(paste(pathstr,"findcoeffminerror.R",sep=""))

#required libraries
library(pracma)  
library(deSolve)
library(Matrix)
library(Rfast)
library(TruncExpFam)

#run example for Ouzel 2, uniformly distributed jams
#site and jam arrangement can be changed in findcoeffminerror.R
findcoeffminerror()
