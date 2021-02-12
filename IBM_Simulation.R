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
## Simulation -------------------
## For each year of burn-in
for(year in 1:nYears){
  ## Pull together all in-migrating adults into one matrix
  Adults=rbind(
    Virgins,
    Repeats
  )
  ## If there are no dams, all adults are assumed
  ## to reach the spawning grounds
  if(nDams==0){
    Spawners=Adults
  } else {
    if(nrow(Adults)>0){
      ## Draw betas from the Maynard et al. 2017 model
      sBetas=betas[sample(nrow(betas),nrow(Adults),replace=TRUE),]
      Passage=exp((sBetas[,1]+sBetas[,2]*(Adults[,2]-65.711)/9.945)+sBetas[,
        3]*runif(nrow(Adults),0,30))
      Passage=ifelse(Passage>1,1,Passage)
      ## Binomial draw based on passage probability at each dam in the system
      Passage=rbinom(nrow(Adults),nDams,Passage)
      ## Spawners are those fish that pass all dams
      Spawners=matrix(Adults[which(Passage>=nDams),],ncol=8)
      ## Calculate the M:F ratio of spawners
      MF=nrow(subset(Spawners,Spawners[,1]==0))/nrow(subset(Spawners,
                                                            Spawners[,1]==1))
      ## Fallbacks are any adults that don't make it past all the dams
      Fallbacks=matrix(Adults[which(Passage<nDams),],ncol=8)
      rm(Adults)
    }
  }
  ## Select which juveniles will mature and which will remain in
  ## salt water
  J1=which(Juveniles1[,3]-2>=Juveniles1[,7])
  J2=which(Juveniles2[,3]-2>=Juveniles2[,7])
  ## Create a new matrix of virgin spawners for next season
  Virgins=rbind(
    Juveniles3,
    Juveniles2[J2,],
    Juveniles1[J1,]
  )
  ## Progress any 2 year juveniles that didn't mature
  Juveniles3=matrix(
    data=Juveniles2[-J2,],
    ncol=8
  )
  ## Apply marine mortality
  Juveniles3=matrix(
    Juveniles3[which(
      rbinom(
        n=nrow(Juveniles3),
        1,
        marineSurvival
      )==1
    ),],
    ncol=8
  )
  ## Progress any 1 year juveniles that didn't mature
  Juveniles2=matrix(
    data=Juveniles1[-J1,],
    ncol=8
  )
  ## Apply marine mortality
  Juveniles2=matrix(
    Juveniles2[which(
      rbinom(
        n=nrow(Juveniles2),
        1,
        marineSurvival
      )==1
    ),],
    ncol=8
  )
  ## Progress any smolts that survive the outmigration
  Smolts=matrix(
    data=Parr2[which(
      rbinom(
        n=nrow(Parr2),
        1,
        0.99^rkm*0.95^nDams
      )==1
    ),],
    ncol=8
  )
  Juveniles1=matrix(
    data=Smolts[which(
      rbinom(
        n=nrow(Smolts),
        1,
        marineSurvival
      )==1
    ),],
    ncol=8
  )
  ## Year 1 Parr and Young of Year Parr progess without any additional mortality
  ##    because all freshwater mortality is applied in the egg --> young of year
  ##    parr transition
  Parr2=Parr1
  Parr1=Parr0
  ## Clean up unused matrices
  rm(Parr0)
  rm(Smolts)
  rm(Adults)
  ## If there are enough spawners, separate out males and females to create a 
  ##    spawning event
  if(nrow(Spawners)>1){
    ## Females
    fSpawners=matrix(
      data=subset(
        Spawners,
        Spawners[,1]==1
      ),
      ncol=8
    )
    ## Males
    mSpawners=matrix(
      data=subset(
        Spawners,
        Spawners[,1]==0
      ),
      ncol=8
    )
    ## If the number of female spawners is equal to the number of male spawners,
    ##    no action is necessary
    ## If there are more female spawners than male spawners, randomly assign
    ##    spawning pairs
    if(nrow(fSpawners)>nrow(mSpawners)){
      a=sample(
        x=nrow(fSpawners),
        size=nrow(mSpawners),
        replace=FALSE,
      )
      fSpawners=matrix(
        data=fSpawners[a,],
        ncol=8
      )
    }
    ## If there are more male spawners than female spawners, assign mates based
    ##    on the advantage given to large 2+SW males in the variable "mating"
    if(nrow(fSpawners)<nrow(mSpawners)){
      mateProb=ifelse(
        mSpawners[,7]==1,
        1,
        mating
      )
      mateProb=mateProb/sum(mateProb)
      a=sample(
        x=nrow(mSpawners),
        size=nrow(fSpawners),
        replace=FALSE,
        prob=mateProb
      )
      mSpawners=matrix(
        data=mSpawners[a,],
        ncol=8
      )
    }
    ## Once you have an even number of female and male fish, assign each fish a 
    ##    z-standardized fork length based on the starting distributions of fork
    ##    lengths specified at the beginning of the simulation
    ## Females
    fFL=vector(
      length=nrow(fSpawners)
    )
    fFL=ifelse(
      fSpawners[,7]==1,
      (fSpawners[,2]-m1)/sd1,
      ifelse(
        fSpawners[,7]==2,
        (fSpawners[,2]-m2)/sd2,
        ifelse(
          fSpawners[,7]==3,
          (fSpawners[,2]-m3)/sd3,
          fFL
        )
      )
    )
    ## Males
    mFL=vector(
      length=nrow(mSpawners)
    )
    mFL=ifelse(
      mSpawners[,7]==1,
      (mSpawners[,2]-m1)/sd1,
      ifelse(
        mSpawners[,7]==2,
        (mSpawners[,2]-m2)/sd2,
        ifelse(
          mSpawners[,7]==3,
          (mSpawners[,2]-m3)/sd3,
          mFL
        )
      )
    )
    ## Calculate the number of eggs produced by each spawning pair based on female
    ##    size, per 
    ## Heinimaa, S. and P. Heinimaa. 2004. Effect of female size on egg quality
    ##    and fecundity of the wild Atlantic Salmon in the sub-Arctic River Teno. 
    ##    Boreal Environmental Research. 9:55-62.
    ##    ISSN: 1239-6095
    ## S is the total number of eggs produced in the system (single number)
    S=sum(
      round(
        exp(3.07*log(fSpawners[,2])-4.46),
        0
      ),
      na.rm=TRUE
    )
    ## spawn is the total number of eggs produced per pair (vector)
    spawn=round(
      exp(3.07*log(fSpawners[,2])-4.46),
      0
    )
    ## If there are no spawning pairs, set both S and spawn to 0
  } else {
    S=0
    spawn=0
  }
  ## As long as there is at least one female fish and one male fish
  if(
    exists("fSpawners")&&
    exists("mSpawners")&&
    nrow(fSpawners)*nrow(mSpawners)!=0
  ){
    ## Pair off the adults
    if(sheritability!=0){
      ## The original code is single commented in this section
      # Sizes=c(
      #   fFL,
      #   mFL
      # )
      ## Averaged size between parents
      #Sizes=(fFL+mFL)/2
      ## Create a vector of potential sizes using the specified size heritability
      # Sizes=Sizes*sheritability+rnorm(length(Sizes),0,1)*(1-sheritability)
      # spawnProb1=rep(
      #   spawn/sum(spawn),
      #   2
      # )
      sizes=cbind(fFL,mFL)
      sizes=cbind(sizes,(sizes[,1]+sizes[,2])/2)
      sizes=cbind(sizes,abs(sizes[,3]-sizes[,2]))
      sizes=cbind(sizes,spawn/sum(spawn))
    }
    if(gheritability!=0){
      ## The original code is single commented in this section
      # Growth=c(
      #   fSpawners[,4],
      #   mSpawners[,4]
      #   )
      Growth=cbind(fSpawners[,4],mSpawners[,4])
      Growth=cbind(Growth,(Growth[,1]+Growth[,2])/2)
      Growth=cbind(Growth,abs(Growth[,3]-Growth[,2]))
      Growth=cbind(Growth,spawn/sum(spawn))
      ## Create a vector of possible growths using the growth heritability value
      # Growth=Growth*gheritability+rnorm(length(Sizes),0,1)*(1-gheritability)
      ## Grilse and 2SW thresholds are also inherited based on growth heritability
      # tGrilse=na.omit(
      #   c(
      #     fSpawners[,5],
      #     mSpawners[,5]
      #   )
      # )
      tGrilse=cbind(fSpawners[,5],mSpawners[,5])
      tGrilse=cbind(tGrilse,(tGrilse[,1]+tGrilse[,2])/2)
      tGrilse=cbind(tGrilse,abs(tGrilse[,3]-tGrilse[,2]))
      tGrilse=cbind(tGrilse,spawn/sum(spawn))
      # tGrilse=tGrilse*gheritability+rnorm(length(Sizes),0,1)*(1-gheritability)
      # t2SW=na.omit(
      #   c(
      #     fSpawners[,6],
      #     mSpawners[,6]
      #   )
      # )
      # t2SW=t2SW*gheritability+rnorm(length(Sizes),0,1)*(1-gheritability)
      t2SW=cbind(fSpawners[,6],mSpawners[,6])
      t2SW=cbind(t2SW,(t2SW[,1]+t2SW[,2])/2)
      t2SW=cbind(t2SW,abs(t2SW[,3]-t2SW[,2]))
      t2SW=cbind(t2SW,spawn/sum(spawn))
      # spawnProb2=rep(spawn/sum(spawn),2)
    }
    Spawners=rbind(
      fSpawners,
      mSpawners
    )
    ## All spawners get one spawning season added to their total
    Spawners[,8]=Spawners[,8]+1
  }
  ## Any reconditioned fish at sea prepare to join the next spawning migration
  Repeats=matrix(
    data=Recon[,1:8],
    ncol=8
  )
  ## All reconditioned fish grow
  Repeats[,2]=Repeats[,2]+sample(
    x=maxGrowth,
    nrow(Repeats),
    replace=TRUE
  )
  ## Clear the unnecessary recon matrix
  rm(Recon)
  ## Create a new recon matrix populated by surviving kelts
  Recon=matrix(
    data=Spawners[which(
      rbinom(
        n=nrow(Spawners),
        1,
        0.63*marineSurvival
      )==1
    ),1:8],
    ncol=8
  )
  Recon=matrix(
    data=Recon[,1:8],
    ncol=8
  )
  if(db==TRUE){
    points(nrow(Spawners)~b)
  }
  ## Clear the old Spawner matrix
  rm(Spawners)
  ## Build a new Young of Year Parr matrix
  Parr0=matrix(
    ncol=8,
    nrow=round(
      S*exp(A-B*S)*freshwaterSurvival,
      0
    )
  )
  if(nrow(Parr0)>0){
    ## Assign each fish a parentage
    parents=sample(
      x=seq(1,nrow(sizes),1),
      size=nrow(Parr0),
      prob=sizes[,5],
      replace=TRUE
    )
    ## Use parentage to assign traits to each individual
    for(pa in unique(parents)){
      Parr0[,2]=ifelse(
        parents==pa,
        rnorm(1,mean=sizes[pa,3],sd=sizes[pa,4]),
        Parr0[,2]
      )
      Parr0[,4]=ifelse(
        parents==pa,
        rnorm(1,mean=Growth[pa,3],sd=Growth[pa,4]),
        Parr0[,4]
      )
      Parr0[,5]=ifelse(
        parents==pa,
        rnorm(1,mean=tGrilse[pa,3],sd=tGrilse[pa,4]),
        Parr0[,5]
      )
      Parr0[,6]=ifelse(
        parents==pa,
        rnorm(1,mean=t2SW[pa,3],sd=t2SW[pa,4]),
        Parr0[,6]
      )
    }
    ## Apply environmental variability to size
    Parr0[,2]=(Parr0[,2]*sheritability)+(runif(
      n=nrow(Parr0),
      min=-1.96,
      max=1.96
    )*(1-sheritability))
    ## Apply environmental variability to growth
    Parr0[,4]=(Parr0[,4]*gheritability)+(runif(
      n=nrow(Parr0),
      min=-1.96,
      max=1.96
    )*(1-gheritability))
    Parr0[,5]=(Parr0[,5]*gheritability)+(runif(
      n=nrow(Parr0),
      min=-1.96,
      max=1.96
    )*(1-gheritability))
    Parr0[,6]=(Parr0[,6]*gheritability)+(runif(
      n=nrow(Parr0),
      min=-1.96,
      max=1.96
    )*(1-gheritability))
    ## Compare thresholds to growth parameter and assign a value for seawinters
    ##    at maturity
    Parr0[,7]=ifelse(
      Parr0[,4]<Parr0[,6],
      3,
      ifelse(
        Parr0[,4]>Parr0[,5],
        1,
        ifelse(
          is.na(Parr0[,7]),
          2,
          Parr0[,7]
        )
      )
    )
    ## Calculate size at maturity for each fish using its z-standardized growth
    ##    and the size distribution selected at the beginning of the simulation
    Sizes=Parr0[,2]
    Parr0[,2]=ifelse(
      Parr0[,7]==1,
      Sizes*sd1+m1,
      ifelse(
        Parr0[,7]==2,
        Sizes*sd2+m2,
        ifelse(
          Parr0[,7]>=3,
          Sizes*sd3+m3,
          Parr0[,2]
        )
      )
    )
    ## Assign each fish a sex based on its assigned SW at maturity
    Parr0[,1]=ifelse(
      Parr0[,7]==1,
      rbinom(
        nrow(Parr0),
        1,
        0.015
      ),
      rbinom(
        nrow(Parr0),
        1,
        0.55
      )
    )
    ## Ensure that 2SW maturation thresholds are lower than 1SW maturation 
    ##    thresholds
    while(sum(Parr0[,6]>Parr0[,5])>0){
      Parr0[,6]=ifelse(
        Parr0[,6]>Parr0[,5],
        Parr0[,6]-0.1,
        Parr0[,6]
      )
    }
    ## All YOY have 0 prior spawns and are 0 years old
    Parr0[,8]=0
    Parr0[,3]=0
  }
  ## Increase the ages of all fish except YOY
  Repeats[,3]=Repeats[,3]+1
  Recon[,3]=Recon[,3]+1
  Juveniles1[,3]=Juveniles1[,3]+1
  Juveniles2[,3]=Juveniles2[,3]+1
  Juveniles3[,3]=Juveniles3[,3]+1
  Parr2[,3]=Parr2[,3]+1
  Parr1[,3]=Parr1[,3]+1
  ## Populate the results matrix
  newRow=vector(length=ncol(Results))
  newRow[1]=i
  newRow[2]=nDams
  newRow[3]=gheritability
  newRow[4]=sheritability
  newRow[5]=nrow(mSpawners)/nrow(fSpawners)
  newRow[6]=year
  newRow[7]=nrow(subset(Spawners,Spawners[,7]==1))
  newRow[8]=nrow(subset(Spawners,Spawners[,7]==2))
  newRow[9]=nrow(subset(Spawners,Spawners[,7]>=3))
  newRow[10]=mean(subset(Spawners,Spawners[,7]==1)[,2])
  newRow[11]=mean(subset(Spawners,Spawners[,7]==2)[,2])
  newRow[12]=mean(subset(Spawners,Spawners[,7]>=3)[,2])
  newRow[13]=sd(subset(Spawners,Spawners[,7]==1)[,2])
  newRow[14]=sd(subset(Spawners,Spawners[,7]==2)[,2])
  newRow[15]=sd(subset(Spawners,Spawners[,7]>=3)[,2])
  newRow[16]=mean(Spawners[,4])
  newRow[17]=mean(Spawners[,5])
  newRow[18]=mean(Spawners[,6])
  newRow[19]=
  ## Move the progress bar forward
  setTxtProgressBar(
    pb=bpb,
    value=b/maxBurn
  )
}