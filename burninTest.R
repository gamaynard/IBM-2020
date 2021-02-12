## ---------------------------
##
## Script name: burninTest.R
##
## Purpose of script: To test whether IBM_burnin.R can produce stable
##    simulated salmon populations.
##
## Author: George A. Maynard
##
## Date Created: 2021-01-24
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
## Tested 2021-01-24, 100 runs to 50 years and all were within the 300-800 range
##    of likely values for the Narraguagus River
## ---------------------------

## Set options
options(scipen = 6, digits = 4) # eliminates scientific notation
## ---------------------------

## load necessary packages:

## ---------------------------

## load necessary functions:

## ---------------------------
db=TRUE
if(db==TRUE){
  plot(
    x=1,
    y=1,
    type='n',
    xlim=c(0,maxBurn),
    ylim=c(0,2500)
  )
  abline(h=c(300,800),lty=2)
}
for(testBurnin in 1:5){
  source("IBM_burnin.R")
}
