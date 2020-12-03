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
## Set a starting number of eggs based on carrying capacity of eggs per
##    m2 of spawning habitat
S=cc*W

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
## Calculate number of parr present based on the number of eggs and
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
## Calculate the number of juveniles based on parr presence, smolt
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
