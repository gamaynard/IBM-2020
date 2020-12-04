## -----------------------------------
##
## Script Name: IBM_burnin.R
##
## Purpose of Script: This script runs as an undammed system for a period
##    of years specified by the user before handing off the virtual 
##    salmon to the simulation script. It ensures that the simulations
##    are all stable before starting. 
## 
## Author: George A. Maynard
## 
## Date Created: 2020-12-03
##
## Copyright (c) George Alphonse Maynard, 2020
## galphonsemaynard@gmail.com
##
## -----------------------------------
##
## Notes: This is the cleaned up version of my dissertation code
##
## -----------------------------------
##
## Set working directory
##
## -----------------------------------
##
## Set options
## Turn off scientific notation
options(
  scipen=6,
  digits=4
)
## -----------------------------------
##
## Load necessary packages
##
## -----------------------------------
##
## Load necessary functions
##
## -----------------------------------
## Clear the workspace of all extraneous objects (i.e., those not
##    specified in the list).
rm(
  list=setdiff(
    x=ls(),
    y=c(
      "betas",
      "Results",
      "sheritability",
      "gheritability",
      "nDams",
      "simNum",
      "mating",
      "nYears",
      "maxBurn",
      "rkm",
      "W",
      "A",
      "B",
      "startSize",
      "freshwaterSurvival",
      "damSurvive",
      "marineSurvival",
      "keltSurvival",
      "maxGrowth",
      "s",
      "i",
      "cc"
    )
  )
)

## Create a progress bar
bpb=txtProgressBar(
  min=0,
  max=1,
  initial=0,
  char=paste0(
    ".",
    s,
    "."
  ),
  style=3
)

## Set starting size distributions based on "startSize"
if(startSize=="modern"){
  m1=55.5
  sd1=4
  m2=76.6
  sd2=5.7
  m3=91.9
  sd3=7.6
} else {
  if(startSize=="historic"){
    m1=75
    sd1=6.59
    m2=83
    sd2=4.51
    m3=99
    sd3=2.07
  }
}

## Set a starting number of eggs based on carrying capacity of eggs per
##    m2 of spawning habitat
S=cc*W

## Initialize number of parr present based on the number of eggs and
##    freshwater survival
## New parr
Parr0=round(
  S*exp(A-B*S)*freshwaterSurvival,
  0
)
## 1 year old parr
Parr1=Parr0
## 2 year old parr
Parr2=Parr1

## Initialize the number of juveniles based on parr presence, smolt
##    survival and at sea survival
## 1SW juveniles
Juveniles1=round(
  Parr2*0.99^rkm*marineSurvival^2,
  0
)
## 2SW juveniles
Juveniles2=round(
  Parr2*0.99^rkm*marineSurvival^3,
  0
)
## 3SW juveniles
Juveniles3=round(
  Parr2*0.99^rkm*marineSurvival^4,
  0
)

## The number of spawning adults is set at the upper limit of the
##    estimated capacity of the Narraguagus River (USASAC 2015)
Adults=800

## The number of reconditioning adults is based on kelt survival in
##    an undammed system and marine survival.
Recon=round(
  Adults*0.63*marineSurvival,
  0
)

## The number of repeat spawn adults is based on the number of 
##    reconditioning adults and marine survival
Repeats=round(
  Recon*marineSurvival,
  0
)

## Create matrices to store each life history stage. These matrices
##    will each have eight columns as follows:
##    1) Sex (0=Male, 1=Female)
##    2) Size (Fork Length measured in cm)
##    3) Age (Years)
##    4) Growth coefficient (z distribution, lower=slower growing)
##    5) Individual threshold to become a grilse (z distribution)
##    6) Individual threshold to become a 2SW (z distribution)
##    7) Seawinters at maturity
##    8) Prior spawns

## Reconditioned spawners --------------------
Recon=matrix(
  ncol=8,
  nrow=Recon
)
colnames(Recon)=c(
  "Sex",
  "FL",
  "Age",
  "Growth",
  "grilseThresh",
  "twoSWThresh",
  "SWatMat",
  "Spawns"
)
## The growth coefficient is a random draw from a z distribution
Recon[,4]=rnorm(
  n=nrow(Recon),
  mean=0,
  sd=1
)
## Fish in the upper 10% of the growth distribution are classified as
##    grilse to match sea age proportions in Izzo and Zydlewski 2017
Recon[,7]=ifelse(
  Recon[,4]>=1.28,
  1,
  2
)
## Fish in the lower 3% of the growth distribution are classified as
##    MSW to match sea age proportions in Izzo and Zydlewski 2017
Recon[,7]=ifelse(
  Recon[,4]<(-1.87),
  3,
  Recon[,7]
)
## Set the sizes of fish based on the starting size distribution
##    selected and at sea growth
Recon[,2]=ifelse(
  Recon[,7]==1,
  rnorm(
    n=nrow(Recon),
    mean=m1,
    sd=sd1
  ),
  ifelse(
    Recon[,7]==2,
    rnorm(
      n=nrow(Recon),
      mean=m2,
      sd=sd2
    ),
    rnorm(
      n=nrow(Recon),
      mean=m3,
      sd=sd3
    )
  )
)
Recon[,2]=Recon[,2]+sample(
  x=maxGrowth,
  size=nrow(Recon),
  replace=TRUE
)
## Set sex of fish based on life history (i.e., most grilse are males
##    and most 2SW and MSW are females)
Recon[,1]=ifelse(
  Recon[,7]==1,
  rbinom(
    n=nrow(Recon),
    size=1,
    prob=0.015
  ),
  rbinom(
    n=nrow(Recon),
    size=1,
    prob=0.55
  )
)
## Set maturity thresholds based on life history
Recon[,c(5,6)]=0
Recon[,5]=ifelse(
  Recon[,7]>=2,
  rtruncnorm(
    n=nrow(Recon),
    a=Recon[,4],
    b=10000,
    mean=0,
    sd=1
  ),
  rtruncnorm(
    n=nrow(Recon),
    a=-10000,
    b=Recon[,4],
    mean=0,
    sd=1
  )
)
Recon[,6]=ifelse(
  Recon[,7]==3,
  rtruncnorm(
    n=nrow(Recon),
    a=-10000,
    b=Recon[,5],
    mean=0,
    sd=1
  ),
  Recon[,6]
)
Recon[,6]=ifelse(
  Recon[,7]==2,
  rtruncnorm(
    n=nrow(Recon),
    a=-10000,
    b=Recon[,4],
    mean=0,
    sd=1
  ),
  Recon[,6]
)
Recon[,6]=ifelse(
  Recon[,7]==1,
  rtruncnorm(
    n=nrow(Recon),
    a=Recon[,4],
    b=10000,
    mean=0,
    sd=1
  ),
  Recon[,6]
)