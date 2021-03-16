## ---------------------------
##
## Script name: IBM_analysis.R
##
## Purpose of script: Analyze the new format of simulation outputs
##
## Author: George A. Maynard
##
## Date Created: 2021-03-03
##
## Copyright (c) George Alphonse Maynard, 2020
## Email: galphonsemaynard@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory
## Set the working directory to the local SAMBA server (different path depending
##    on machine)
if(Sys.info()[[4]]=="Acanthias"){
  setwd("/run/user/1000/gvfs/smb-share:server=salmo.local,share=sambashare/george/NewData")
}
## ---------------------------

## Set options
options(scipen = 6, digits = 4) # eliminates scientific notation
## ---------------------------

## load necessary packages:
library(qpcR)
## ---------------------------

## load necessary functions:

## ---------------------------
## Read in the data
historic=read.csv("historic.csv")
modern=read.csv("modern.csv")
## Add a column to differentiate them
historic$start="historic"
modern$start="modern"
## Combine the frames
data=rbind(historic,modern)
## subset out only those simulations that ran and exported data correctly
data=subset(
  data,
  data$nDams%in%seq(0,5,1)
)
data=subset(
  data,
  data$Year%in%seq(1,100,1)
)
## Format the data to numbers
data$SimNum=as.numeric(as.character(data$SimNum))
data$nDams=as.numeric(as.character(data$nDams))
data$gHerit=as.numeric(as.character(data$gHerit))
data$sHerit=as.numeric(as.character(data$sHerit))
data$M.F=as.numeric(as.character(data$M.F))
data$Year=as.numeric(as.character(data$Year))
data$n1SW=as.numeric(as.character(data$n1SW))
data$n2SW=as.numeric(as.character(data$n2SW))
data$nMSW=as.numeric(as.character(data$nMSW))
data$FL1SW=as.numeric(as.character(data$FL1SW))
data$FL2SW=as.numeric(as.character(data$FL2SW))
data$FL3SW=as.numeric(as.character(data$FL3SW))
data$sd1SW=as.numeric(as.character(data$sd1SW))
data$sd2SW=as.numeric(as.character(data$sd2SW))
data$sd3SW=as.numeric(as.character(data$sd3SW))
data$mGrowth=as.numeric(as.character(data$mGrowth))
data$mGrilse=as.numeric(as.character(data$mGrilse))
data$m2SW=as.numeric(as.character(data$m2SW))
data$S1SW=as.numeric(as.character(data$S1SW))
data$S2SW=as.numeric(as.character(data$S2SW))
data$SMSW=as.numeric(as.character(data$SMSW))
##### ANALYSIS OF SPAWNING POPULATION #####
data$Spawners=data$n1SW+data$n2SW+data$nMSW

##### ANALYSIS OF END OF SIMULATION RESULTS #####
results=subset(data,data$Year==100)
results=subset(results,results$Spawners>0)
H=subset(results,results$start=="historic")
M=subset(results,results$start=="modern")
## Historic 1SW
mod1=glm(
  H$FL1SW~H$nDams+H$gHerit+H$sHerit
)
mod2=glm(
  H$FL1SW~H$nDams+H$gHerit
)
mod3=glm(
  H$FL1SW~H$nDams+H$sHerit
)
mod4=glm(
  H$FL1SW~H$nDams
)
mod5=glm(
  H$FL1SW~H$sHerit+H$gHerit
)
mod6=glm(
  H$FL1SW~H$sHerit
)
mod7=glm(
  H$FL1SW~H$gHerit
)
modList=list(
  mod1,
  mod2,
  mod3,
  mod4,
  mod5,
  mod6,
  mod7
)
scores=vector(length=length(modList))
for(i in 1:length(modList)){
  scores[i]=AICc(modList[[i]])
}
akaike.weights(scores)

## Historic 2SW
mod1=glm(
  H$FL2SW~H$nDams+H$gHerit+H$sHerit
)
mod2=glm(
  H$FL2SW~H$nDams+H$gHerit
)
mod3=glm(
  H$FL2SW~H$nDams+H$sHerit
)
mod4=glm(
  H$FL2SW~H$nDams
)
mod5=glm(
  H$FL2SW~H$sHerit+H$gHerit
)
mod6=glm(
  H$FL2SW~H$sHerit
)
mod7=glm(
  H$FL2SW~H$gHerit
)
modList=list(
  mod1,
  mod2,
  mod3,
  mod4,
  mod5,
  mod6,
  mod7
)
scores=vector(length=length(modList))
for(i in 1:length(modList)){
  scores[i]=AICc(modList[[i]])
}
akaike.weights(scores)

## Historic MSW
mod1=glm(
  H$FL3SW~H$nDams+H$gHerit+H$sHerit
)
mod2=glm(
  H$FL3SW~H$nDams+H$gHerit
)
mod3=glm(
  H$FL3SW~H$nDams+H$sHerit
)
mod4=glm(
  H$FL3SW~H$nDams
)
mod5=glm(
  H$FL3SW~H$sHerit+H$gHerit
)
mod6=glm(
  H$FL3SW~H$sHerit
)
mod7=glm(
  H$FL3SW~H$gHerit
)
modList=list(
  mod1,
  mod2,
  mod3,
  mod4,
  mod5,
  mod6,
  mod7
)
scores=vector(length=length(modList))
for(i in 1:length(modList)){
  scores[i]=AICc(modList[[i]])
}
akaike.weights(scores)

## Modern 1SW
mod1=glm(
  M$FL1SW~M$nDams+M$gHerit+M$sHerit
)
mod2=glm(
  M$FL1SW~M$nDams+M$gHerit
)
mod3=glm(
  M$FL1SW~M$nDams+M$sHerit
)
mod4=glm(
  M$FL1SW~M$nDams
)
mod5=glm(
  M$FL1SW~M$sHerit+M$gHerit
)
mod6=glm(
  M$FL1SW~M$sHerit
)
mod7=glm(
  M$FL1SW~M$gHerit
)
modList=list(
  mod1,
  mod2,
  mod3,
  mod4,
  mod5,
  mod6,
  mod7
)
scores=vector(length=length(modList))
for(i in 1:length(modList)){
  scores[i]=AICc(modList[[i]])
}
akaike.weights(scores)

## Modern 2SW
mod1=glm(
  M$FL2SW~M$nDams+M$gHerit+M$sHerit
)
mod2=glm(
  M$FL2SW~M$nDams+M$gHerit
)
mod3=glm(
  M$FL2SW~M$nDams+M$sHerit
)
mod4=glm(
  M$FL2SW~M$nDams
)
mod5=glm(
  M$FL2SW~M$sHerit+M$gHerit
)
mod6=glm(
  M$FL2SW~M$sHerit
)
mod7=glm(
  M$FL2SW~M$gHerit
)
modList=list(
  mod1,
  mod2,
  mod3,
  mod4,
  mod5,
  mod6,
  mod7
)
scores=vector(length=length(modList))
for(i in 1:length(modList)){
  scores[i]=AICc(modList[[i]])
}
akaike.weights(scores)

## Modern MSW
mod1=glm(
  M$FL3SW~M$nDams+M$gHerit+M$sHerit
)
mod2=glm(
  M$FL3SW~M$nDams+M$gHerit
)
mod3=glm(
  M$FL3SW~M$nDams+M$sHerit
)
mod4=glm(
  M$FL3SW~M$nDams
)
mod5=glm(
  M$FL3SW~M$sHerit+M$gHerit
)
mod6=glm(
  M$FL3SW~M$sHerit
)
mod7=glm(
  M$FL3SW~M$gHerit
)
modList=list(
  mod1,
  mod2,
  mod3,
  mod4,
  mod5,
  mod6,
  mod7
)
scores=vector(length=length(modList))
for(i in 1:length(modList)){
  scores[i]=AICc(modList[[i]])
}
akaike.weights(scores)

##### Plotting FL over time #####
plot(data$FL2SW~data$Year,type="n")
points(subset(data,data$start=="historic")$FL2SW~subset(data,data$start=="historic")$Year,col='blue')
points(subset(data,data$start=="modern")$FL2SW~subset(data,data$start=="modern")$Year,col='red')
