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
  Parr2*0.99^rkm*marineSurvival^1,
  0
)
## 2SW juveniles
Juveniles2=round(
  Parr2*0.99^rkm*marineSurvival^2,
  0
)
## 3SW juveniles
Juveniles3=round(
  Parr2*0.99^rkm*marineSurvival^3,
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
## For the purposes of initializing the model, all reconditioned 
##    spawners are assumed to have one prior spawning event
Recon[,8]=1
## Repeat Spawners  --------------------
Repeats=matrix(
  ncol=8,
  nrow=Repeats
)
## The growth coefficient is a random draw from a z distribution
Repeats[,4]=rnorm(
  n=nrow(Repeats),
  mean=0,
  sd=1
)
## Fish in the upper 10% of the distribution are grilse
Repeats[,7]=2
Repeats[,7]=ifelse(
  Repeats[,4]>=1.28,
  1,
  Repeats[,7]
)
## Fish in the lower 3% of the distribution are MSW
Repeats[,7]=ifelse(
  Repeats[,4]<(-1.87),
  3,
  Repeats[,7]
)
## Set the sizes of fish based on the starting size distribution
##    selected and at sea growth
Repeats[,2]=ifelse(
  Repeats[,7]==1,
  rnorm(
    n=nrow(Repeats),
    mean=m1,
    sd=sd1
  ),
  ifelse(
    Repeats[,7]==2,
    rnorm(
      n=nrow(Repeats),
      mean=m2,
      sd=sd2
    ),
    rnorm(
      n=nrow(Repeats),
      mean=m3,
      sd=sd3
    )
  )
)
Repeats[,2]=Repeats[,2]+sample(
  x=maxGrowth,
  size=nrow(Repeats),
  replace=TRUE
)
## Set sex of fish based on life history (i.e., most grilse are males
##    and most 2SW and MSW are females)
Repeats[,1]=ifelse(
  Repeats[,7]==1,
  rbinom(
    n=nrow(Repeats),
    size=1,
    prob=0.015
  ),
  rbinom(
    n=nrow(Repeats),
    size=1,
    prob=0.55
  )
)
## Set maturity thresholds based on life history
Repeats[,c(5,6)]=0
Repeats[,5]=ifelse(
  Repeats[,7]>=2,
  rtruncnorm(
    n=nrow(Repeats),
    a=Repeats[,4],
    b=10000,
    mean=0,
    sd=1
  ),
  rtruncnorm(
    n=nrow(Repeats),
    a=-10000,
    b=Repeats[,4],
    mean=0,
    sd=1
  )
)
Repeats[,6]=ifelse(
  Repeats[,7]==3,
  rtruncnorm(
    n=nrow(Repeats),
    a=-10000,
    b=Repeats[,5],
    mean=0,
    sd=1
  ),
  Repeats[,6]
)
Repeats[,6]=ifelse(
  Repeats[,7]==2,
  rtruncnorm(
    n=nrow(Repeats),
    a=-10000,
    b=Repeats[,4],
    mean=0,
    sd=1
  ),
  Repeats[,6]
)
Repeats[,6]=ifelse(
  Repeats[,7]==1,
  rtruncnorm(
    n=nrow(Repeats),
    a=Repeats[,4],
    b=10000,
    mean=0,
    sd=1
  ),
  Repeats[,6]
)
## For the purposes of initializing the model, all repeat 
##    spawners are assumed to have one prior spawning event
Repeats[,8]=1

## Virgin Adults --------------------
## The number of virgin spawners is the number of adults minus the 
##    number of repeat spawners present
Virgins=Adults-nrow(Repeats)
Virgins=matrix(
  ncol=8,
  nrow=Virgins
)
## The growth coefficient is a random draw from a z distribution
Virgins[,4]=rnorm(
  n=nrow(Virgins),
  mean=0,
  sd=1
)
## Fish in the upper 10% of the distribution are grilse
Virgins[,7]=2
Virgins[,7]=ifelse(
  Virgins[,4]>=1.28,
  1,
  Virgins[,7]
)
## Fish in the lower 3% of the distribution are MSW
Virgins[,7]=ifelse(
  Virgins[,4]<(-1.87),
  3,
  Virgins[,7]
)
## Set the sizes of fish based on the starting size distribution
##    selected and at sea growth
Virgins[,2]=ifelse(
  Virgins[,7]==1,
  rnorm(
    n=nrow(Virgins),
    mean=m1,
    sd=sd1
  ),
  ifelse(
    Virgins[,7]==2,
    rnorm(
      n=nrow(Virgins),
      mean=m2,
      sd=sd2
    ),
    rnorm(
      n=nrow(Virgins),
      mean=m3,
      sd=sd3
    )
  )
)
Virgins[,2]=Virgins[,2]+sample(
  x=maxGrowth,
  size=nrow(Virgins),
  replace=TRUE
)
## Set sex of fish based on life history (i.e., most grilse are males
##    and most 2SW and MSW are females)
Virgins[,1]=ifelse(
  Virgins[,7]==1,
  rbinom(
    n=nrow(Virgins),
    size=1,
    prob=0.015
  ),
  rbinom(
    n=nrow(Virgins),
    size=1,
    prob=0.55
  )
)
## Set maturity thresholds based on life history
Virgins[,c(5,6)]=0
Virgins[,5]=ifelse(
  Virgins[,7]>=2,
  rtruncnorm(
    n=nrow(Virgins),
    a=Virgins[,4],
    b=10000,
    mean=0,
    sd=1
  ),
  rtruncnorm(
    n=nrow(Virgins),
    a=-10000,
    b=Virgins[,4],
    mean=0,
    sd=1
  )
)
Virgins[,6]=ifelse(
  Virgins[,7]==3,
  rtruncnorm(
    n=nrow(Virgins),
    a=-10000,
    b=Virgins[,5],
    mean=0,
    sd=1
  ),
  Virgins[,6]
)
Virgins[,6]=ifelse(
  Virgins[,7]==2,
  rtruncnorm(
    n=nrow(Virgins),
    a=-10000,
    b=Virgins[,4],
    mean=0,
    sd=1
  ),
  Virgins[,6]
)
Virgins[,6]=ifelse(
  Virgins[,7]==1,
  rtruncnorm(
    n=nrow(Virgins),
    a=Virgins[,4],
    b=10000,
    mean=0,
    sd=1
  ),
  Virgins[,6]
)
## By definition, virgin spawners have no previous spawning events
Virgins[,8]=0
## 3SW Juveniles -------------------
Juveniles3=matrix(
  ncol=8,
  nrow=Juveniles3
)
## Growth coefficient should be in the lowest 3%
Juveniles3[,4]=rtruncnorm(
  n=nrow(Juveniles3),
  a=-10000,
  b=-1.87,
  mean=0,
  sd=1
)
## The grilse threshold must be higher than the growth coefficient
Juveniles3[,5]=rtruncnorm(
  n=nrow(Juveniles3),
  a=Juveniles3[,4],
  b=Inf,
  mean=0,
  sd=1
)
## The 2SW threshold should be lower than the grilse threshold, but
##    higher than the growth value
Juveniles3[,6]=rtruncnorm(
  n=nrow(Juveniles3),
  a=Juveniles3[,4],
  b=Juveniles3[,5],
  mean=0,
  sd=1
)
## 3SW juveniles mature at 3SW
Juveniles3[,7]=3
## Which means they are 5 years old
Juveniles3[,3]=5
## Fork length is drawn from the distributions specified at the start
Juveniles3[,2]=rnorm(
  n=nrow(Juveniles3),
  mean=m3,
  sd=sd3
)
## The sex distribution of 3SW fish skews towards females (55:45, F:M)
Juveniles3[,1]=rbinom(
  n=nrow(Juveniles3),
  size=1,
  prob=0.55
)
## Juveniles have never spawned previously
Juveniles3[,8]=0

## 2SW Juveniles -------------------
Juveniles2=matrix(
  ncol=8,
  nrow=Juveniles2
)
## Growth coefficient should be set to mature as a 3SW or 2SW
Juveniles2[,4]=rtruncnorm(
  n=nrow(Juveniles2),
  a=-10000,
  b=1.28,
  mean=0,
  sd=1
)
## The grilse threshold must be higher than the growth coefficient
Juveniles2[,5]=rtruncnorm(
  n=nrow(Juveniles2),
  a=Juveniles2[,4],
  b=10000,
  mean=0,
  sd=1
)
## The 2SW threshold should be lower than the grilse threshold
Juveniles2[,6]=rtruncnorm(
  n=nrow(Juveniles2),
  a=Juveniles2[,5],
  b=10000,
  mean=0,
  sd=1
)
## 2SW juveniles mature at 3SW or 2SW depending on their growth and 
##    maturation threshold levels
Juveniles2[,7]=ifelse(
  Juveniles2[,4]<Juveniles2[,6],
  3,
  2
)
## They are currently 4 years old 
Juveniles2[,3]=4
## Fork length is drawn from the distributions specified at the start
Juveniles2[,2]=ifelse(
  Juveniles2[,7]==2,
  rnorm(
    n=nrow(Juveniles2),
    mean=m2,
    sd=sd2
  ),
  ifelse(
    Juveniles2[,7]==3,
    rnorm(
      nrow(Juveniles2),
      m3,
      sd3
    ),
    Juveniles[,2]
  )
)
## The sex distribution of 2SW fish skews towards females (55:45, F:M)
Juveniles2[,1]=rbinom(
  n=nrow(Juveniles2),
  size=1,
  prob=0.55
)
## Juveniles have never spawned previously
Juveniles2[,8]=0

## 1SW Juveniles -------------------
Juveniles1=matrix(
  ncol=8,
  nrow=Juveniles1
)
## Growth coefficient can be anywhere in the z distribution
Juveniles1[,4]=rnorm(
  n=nrow(Juveniles1),
  mean=0,
  sd=1
)
## The fastest growing fish mature as grilse, moderate growers
##    mature as 2SW fish and slow growers mature as 3SW fish
Juveniles1[,7]=ifelse(
  Juveniles1[,4]>=1.28,
  1,
  ifelse(
    Juveniles1[,4]>-1.87&Juveniles1[,4]<1.28,
    2,
    ifelse(
      Juveniles1[,4]<(-1.87),
      3,
      Juveniles1[,4]
    )
  )
)
## These fish are three years old
Juveniles1[,3]=3
## Fork length is drawn from the distributions specified at the start
Juveniles1[,2]=ifelse(
  Juveniles1[,7]==2,
  rnorm(
    n=nrow(Juveniles1),
    mean=m2,
    sd=sd2
  ),
  ifelse(
    Juveniles1[,7]==3,
    rnorm(
      nrow(Juveniles1),
      m3,
      sd3
    ),
    ifelse(
      Juveniles1[,7]==1,
      rnorm(
        nrow(Juveniles1),
        m1,
        sd1
      ),
      Juveniles1[,7]
    )
  )
)
## The sex distribution of 2+SW fish skews towards females (55:45, 
##    F:M), but the sex distribution of 1SW fish skews towards males
##    (98.5:1.5)
Juveniles1[,1]=ifelse(
  Juveniles1[,7]>=2,
  rbinom(
    n=nrow(Juveniles1),
    size=1,
    prob=0.55
  ),
  ifelse(
    Juveniles1[,7]==1,
    rbinom(
      n=nrow(Juveniles1),
      size=1,
      prob=0.015
    ),
    Juveniles1[,7]
  )
)
## Juveniles have never spawned previously
Juveniles1[,8]=0
## Assign maturity thresholds
Juveniles1[,6]=ifelse(
  Juveniles1[,7]==3,
  rtruncnorm(
    n=nrow(Juveniles1),
    a=Juveniles1[,4],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Juveniles1[,7]==2,
    rtruncnorm(
      n=nrow(Juveniles1),
      a=-Inf,
      b=Juveniles1[,4],
      mean=0,
      sd=1
    ),
    ifelse(
      Juveniles1[,7]==1,
      rtruncnorm(
        n=nrow(Juveniles1),
        a=-Inf,
        b=Juveniles1[,4],
        mean=0,
        sd=1
      ),
      Juveniles1[,6]
    )
  )
)
Juveniles1[,5]=ifelse(
  Juveniles1[,7]>=2,
  rtruncnorm(
    n=nrow(Juveniles1),
    a=Juveniles1[,6],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Juveniles1[,7]==1,
    rtruncnorm(
      n=nrow(Juveniles1),
      a=Juveniles1[,6],
      b=Juveniles1[,4],
      mean=0,
      sd=1
    ),
    Juveniles1[,5]
  )
)
## 2 year parr -------------------
Parr2=matrix(
  ncol=8,
  nrow=Parr2
)
## Growth coefficient can be anywhere in the z distribution
Parr2[,4]=rnorm(
  n=nrow(Parr2),
  mean=0,
  sd=1
)
## The fastest growing fish mature as grilse, moderate growers
##    mature as 2SW fish and slow growers mature as 3SW fish
Parr2[,7]=ifelse(
  Parr2[,4]>=1.28,
  1,
  ifelse(
    Parr2[,4]>-1.87&Parr2[,4]<1.28,
    2,
    ifelse(
      Parr2[,4]<(-1.87),
      3,
      Parr2[,4]
    )
  )
)
## These fish are two years old
Parr2[,3]=2
## Fork length is drawn from the distributions specified at the start
Parr2[,2]=ifelse(
  Parr2[,7]==2,
  rnorm(
    n=nrow(Parr2),
    mean=m2,
    sd=sd2
  ),
  ifelse(
    Parr2[,7]==3,
    rnorm(
      nrow(Parr2),
      m3,
      sd3
    ),
    ifelse(
      Parr2[,7]==1,
      rnorm(
        nrow(Parr2),
        m1,
        sd1
      ),
      Parr2[,7]
    )
  )
)
## The sex distribution of 2+SW fish skews towards females (55:45, 
##    F:M), but the sex distribution of 1SW fish skews towards males
##    (98.5:1.5)
Parr2[,1]=ifelse(
  Parr2[,7]>=2,
  rbinom(
    n=nrow(Parr2),
    size=1,
    prob=0.55
  ),
  ifelse(
    Parr2[,7]==1,
    rbinom(
      n=nrow(Parr2),
      size=1,
      prob=0.015
    ),
    Parr2[,7]
  )
)
## Juveniles have never spawned previously
Parr2[,8]=0
## Assign maturity thresholds
Parr2[,6]=ifelse(
  Parr2[,7]==3,
  rtruncnorm(
    n=nrow(Parr2),
    a=Parr2[,4],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Parr2[,7]==2,
    rtruncnorm(
      n=nrow(Parr2),
      a=-Inf,
      b=Parr2[,4],
      mean=0,
      sd=1
    ),
    ifelse(
      Parr2[,7]==1,
      rtruncnorm(
        n=nrow(Parr2),
        a=-Inf,
        b=Parr2[,4],
        mean=0,
        sd=1
      ),
      Parr2[,6]
    )
  )
)
Parr2[,5]=ifelse(
  Parr2[,7]>=2,
  rtruncnorm(
    n=nrow(Parr2),
    a=Parr2[,6],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Parr2[,7]==1,
    rtruncnorm(
      n=nrow(Parr2),
      a=Parr2[,6],
      b=Parr2[,4],
      mean=0,
      sd=1
    ),
    Parr2[,5]
  )
)

## 1 year parr -------------------
Parr1=matrix(
  ncol=8,
  nrow=Parr1
)
## Growth coefficient can be anywhere in the z distribution
Parr1[,4]=rnorm(
  n=nrow(Parr1),
  mean=0,
  sd=1
)
## The fastest growing fish mature as grilse, moderate growers
##    mature as 2SW fish and slow growers mature as 3SW fish
Parr1[,7]=ifelse(
  Parr1[,4]>=1.28,
  1,
  ifelse(
    Parr1[,4]>-1.87&Parr1[,4]<1.28,
    2,
    ifelse(
      Parr1[,4]<(-1.87),
      3,
      Parr1[,4]
    )
  )
)
## These fish are one year old
Parr1[,3]=1
## Fork length is drawn from the distributions specified at the start
Parr1[,2]=ifelse(
  Parr1[,7]==2,
  rnorm(
    n=nrow(Parr1),
    mean=m2,
    sd=sd2
  ),
  ifelse(
    Parr1[,7]==3,
    rnorm(
      nrow(Parr1),
      m3,
      sd3
    ),
    ifelse(
      Parr1[,7]==1,
      rnorm(
        nrow(Parr1),
        m1,
        sd1
      ),
      Parr1[,7]
    )
  )
)
## The sex distribution of 2+SW fish skews towards females (55:45, 
##    F:M), but the sex distribution of 1SW fish skews towards males
##    (98.5:1.5)
Parr1[,1]=ifelse(
  Parr1[,7]>=2,
  rbinom(
    n=nrow(Parr1),
    size=1,
    prob=0.55
  ),
  ifelse(
    Parr1[,7]==1,
    rbinom(
      n=nrow(Parr1),
      size=1,
      prob=0.015
    ),
    Parr1[,7]
  )
)
## Juveniles have never spawned previously
Parr1[,8]=0
## Assign maturity thresholds
Parr1[,6]=ifelse(
  Parr1[,7]==3,
  rtruncnorm(
    n=nrow(Parr1),
    a=Parr1[,4],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Parr1[,7]==2,
    rtruncnorm(
      n=nrow(Parr1),
      a=-Inf,
      b=Parr1[,4],
      mean=0,
      sd=1
    ),
    ifelse(
      Parr1[,7]==1,
      rtruncnorm(
        n=nrow(Parr1),
        a=-Inf,
        b=Parr1[,4],
        mean=0,
        sd=1
      ),
      Parr1[,6]
    )
  )
)
Parr1[,5]=ifelse(
  Parr1[,7]>=2,
  rtruncnorm(
    n=nrow(Parr1),
    a=Parr1[,6],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Parr1[,7]==1,
    rtruncnorm(
      n=nrow(Parr1),
      a=Parr1[,6],
      b=Parr1[,4],
      mean=0,
      sd=1
    ),
    Parr1[,5]
  )
)

## Young of year parr -------------------
Parr0=matrix(
  ncol=8,
  nrow=Parr0
)
## Growth coefficient can be anywhere in the z distribution
Parr0[,4]=rnorm(
  n=nrow(Parr0),
  mean=0,
  sd=1
)
## The fastest growing fish mature as grilse, moderate growers
##    mature as 2SW fish and slow growers mature as 3SW fish
Parr0[,7]=ifelse(
  Parr0[,4]>=1.28,
  1,
  ifelse(
    Parr0[,4]>-1.87&Parr0[,4]<1.28,
    2,
    ifelse(
      Parr0[,4]<(-1.87),
      3,
      Parr0[,4]
    )
  )
)
## These fish hatched this year
Parr0[,3]=0
## Fork length is drawn from the distributions specified at the start
Parr0[,2]=ifelse(
  Parr0[,7]==2,
  rnorm(
    n=nrow(Parr0),
    mean=m2,
    sd=sd2
  ),
  ifelse(
    Parr0[,7]==3,
    rnorm(
      nrow(Parr0),
      m3,
      sd3
    ),
    ifelse(
      Parr0[,7]==1,
      rnorm(
        nrow(Parr0),
        m1,
        sd1
      ),
      Parr0[,7]
    )
  )
)
## The sex distribution of 2+SW fish skews towards females (55:45, 
##    F:M), but the sex distribution of 1SW fish skews towards males
##    (98.5:1.5)
Parr0[,1]=ifelse(
  Parr0[,7]>=2,
  rbinom(
    n=nrow(Parr0),
    size=1,
    prob=0.55
  ),
  ifelse(
    Parr0[,7]==1,
    rbinom(
      n=nrow(Parr0),
      size=1,
      prob=0.015
    ),
    Parr0[,7]
  )
)
## Juveniles have never spawned previously
Parr0[,8]=0
## Assign maturity thresholds
Parr0[,6]=ifelse(
  Parr0[,7]==3,
  rtruncnorm(
    n=nrow(Parr0),
    a=Parr0[,4],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Parr0[,7]==2,
    rtruncnorm(
      n=nrow(Parr0),
      a=-Inf,
      b=Parr0[,4],
      mean=0,
      sd=1
    ),
    ifelse(
      Parr0[,7]==1,
      rtruncnorm(
        n=nrow(Parr0),
        a=-Inf,
        b=Parr0[,4],
        mean=0,
        sd=1
      ),
      Parr0[,6]
    )
  )
)
Parr0[,5]=ifelse(
  Parr0[,7]>=2,
  rtruncnorm(
    n=nrow(Parr0),
    a=Parr0[,6],
    b=Inf,
    mean=0,
    sd=1
  ),
  ifelse(
    Parr0[,7]==1,
    rtruncnorm(
      n=nrow(Parr0),
      a=Parr0[,6],
      b=Parr0[,4],
      mean=0,
      sd=1
    ),
    Parr0[,5]
  )
)
