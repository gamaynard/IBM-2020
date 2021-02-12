## -----------------------------------
##
## Script Name: IBM_header.R
##
## Purpose of Script: Codes user inputs and assumptions for the salmon
##    evolution model, creates an empty frame to store results, and 
##    calls burn-in and simulation functions
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
## Set working directory: 
##    The working directory is automatically set to the local GitHub
##    directory
## -----------------------------------
##
## Set options
## Turn off scientific notation
options(
  scipen=6,
  digits=4
)
## Do you want diagnostic plots periodically?
db=TRUE
## -----------------------------------
##
## Load necessary packages
library(truncnorm)
## -----------------------------------
##
## Load necessary functions
##
## -----------------------------------
##
## Load necessary datasets
## Read in the beta values from Milford Dam from the paper
## Maynard G.A., M.T. Kinnison, and J.D. Zydlewski. 2017. Size selection
##    from fishways and potential evolutionary responses in a threatened
##    Atlantic Salmon population. River Research and Applications.
##    33(7) 1004-1015. DOI: 10.1002/rra.3155
betas=read.csv(
  file="/mnt/md0/Manuscripts/2020 IBM/MFBetas.csv"
  )
## Remove any NA values from the beta dataframe
betas=na.omit(betas)
## -----------------------------------
## SECTION 1: USER INPUTS
## Narrow sense heritability of size at age (must be between 0 and 1)
sheritability=0.2
## Narrow sense heritability of growth rate (must be between 0 and 1)
gheritability=0.2
## Number of dams in the system
nDams=0
## Number of simulations to run
simNum=1
## Number of years to simulate
nYears=100
## Length of burn-in (when the system will run without dams for a time
##    in order to stabilize). Time is in years. 
maxBurn=100
## Parameterize the river system length (km) and wetted habitat area
##    (m2). These values were taken from the description of the 
##    Narraguagus River in the 2015 Annual Report of the U.S. Atlantic
##    Salmon Committee, which can be found at
##    http://www.nefsc.noaa.gov/USASAC/Reports/USASAC2015-Report%2327-2014-Activities.pdf
rkm=89
W=502763
## -----------------------------------
## SECTION 2: ASSUMPTIONS
## Set the advantage for >=2SW males in the mating system. For example
##    if mating is "2", then the >=2SW males are twice as likely to be
##    drawn from the pool as grilse. If mating is "3", they are three
##    times as likely, etc. This values is not limited to whole numbers
mating=1.92
## Parameterize the Ricker Curve, values based on:
## Prevost, E. 2003. Setting biological reference points for Atlantic
##    Salmon stocks: transfer of information from data-rich to
##    sparse-data situation by Bayesian hierarchical modelling. ICES
##    Journal of Marine Science. 60(6):1177-1193.
##    DOI: 10.1016/j.icesjms.2003.08.001
A=0.420250111
B=0.0000005
## Set the starting distributions of fish sizes to either "modern" or
##    "historic". "modern" sizes are based on fish measured by the 
##    Maine Dept. of Marine Resources at the Veazie Dam fish trap on
##    the Penobscot River between 1979-2015. "historic" sizes are lengths
##    that have been back-calculated from weights measured by U.S.
##    Commission on Fish and Fisheries Agents at the Bucksport, Maine
##    fish market at the mouth of the Penobscot River in 1879. Lengths
##    were back-calculated using the average of three Atlantic Salmon
##    length/weight curves (DFO Canada, Icelandic Institute for 
##    Freshwater Science, and the River Tweed Commission of Scotland). 
startSize="historic"
## Set survival values
## Hatching to smolt survival is taken from 
## Legault, Christopher M. 2005. Population viability analysis of 
##    Atlantic Salmon in Maine, USA. Transactions of the American 
##    Fisheries Society. 134(3):549-562. 
##    DOI: 10.1577/T04-017.1
freshwaterSurvival=0.041426
## Per river kilometer survival of smolts at dams is based off the 
##    target value for downstream survival at hydroelectric dams in
##    the Penobscot River system
damSurvive=0.95
## Annual survival at sea is taken from
## Chaput, G. 2012. Overview of the status of Atlantic Salmon (Salmo
##    salar) in the North Atlantic and trends in marine mortality. ICES
##    Journal of Marine Science. 69(9):1538-1548
##    DOI: 10.1093/icesjms/fss013
marineSurvival=0.24
## Outmigration survival for kelts is based on the number of dams
##    present and is taken from 
## Maynard, George A., Lisa K. Izzo, and Joseph D. Zydlewski. 2018. 
##    Movement and mortality of Atlantic Salmon kelts (Salmo salar)
##    released into the Penobscot River, Maine. Fishery Bulletin. 
##    116(3-4):281-290. DOI: 10.755/FB.116.3-4.6
if(nDams==0){
  keltSurvival=0.63
} else {
  keltSurvival=0.429^nDams
}
## At sea growth for adults is assumed to be between 0 and 5 cm annually
##    based on
## Izzo, Lisa K. and Joseph D. Zydlewski. 2018. Retrospective analysis
##    of seasonal ocean growth rates of two sea winter Atlantic Salmon 
##    in Eastern Maine using historic scales. Marine and Coastal Fisheries.
##    9(1):357-372. 
##    DOI: 10.1080/19425120.2017.1334723
## In the simulation, at sea growth is treated as a draw from a uniform
##    distribution between 0 and maxGrowth
maxGrowth=5
## Carrying capacity of 6 eggs per square meter of spawning habitat is
##    based on
## Jonsson, N., B. Jonsson, and L.P. Hansen. 1998. The relative role
##    of density-dependent and density-independent survival in the life
##    cycle of Atlantic Salmon Salmo salar. Journal of Animal Ecology.
##    67(5):751-762. 
##    DOI: 10.1046/j.1365-2656.1998.00237.x
cc=6
## -----------------------------------
## SECTION 3: SET UP
## WARNING: THE USER SHOULD NOT CHANGE ANY OF THESE VALUES
## Set the starting simulation number to "1"
s=1
## Create a matrix to store simulation results
Results=matrix(
  nrow=simNum*nYears,
  ncol=21
)
colnames(Results)=c(
  "SimNum",
  "nDams",
  "gHerit",
  "sHerit",
  "M:F",
  "Year",
  "n1SW",
  "n2SW",
  "nMSW",
  "FL1SW",
  "FL2SW",
  "FL3SW",
  "sd1SW",
  "sd2SW",
  "sd3SW",
  "mGrowth",
  "mGrilse",
  "m2SW",
  "S1SW",
  "S2SW",
  "SMSW"
)
## For each simulation, loop over the burn-in and simulation
## scripts and iterate the simulation number
for(i in 1:simNum){
  source("IBM_burnin.R")
  source("IBM_simulation.R")
  s=s+1
}
## Write the results to a .csv file
write.csv(
  x=Results,
  file=paste0(
    "g",
    gheritability*10,
    "s",
    sheritability*10,
    "Results",
    ndams,
    "Dams.csv"
  ),
  row.names=FALSE
)