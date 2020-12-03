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
  file="MFBetas.csv"
  )
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
## -----------------------------------
## SECTION 2: ASSUMPTIONS
##
