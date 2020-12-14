**DETAIL ANY MAJOR CHANGES MADE IN THE NEW SCRIPT**

**DAM SURVIVAL IN BURN-IN**

In the original burnin.R file, initialized juvenile survival was subjected to outmigrating mortality at dams as specified in the header file. This would have the effect of reducing the number of juveniles available to start the "natural" population in the burn-in period. However, the outmigration mortality at dams from the header file was not used during the burn-in simulation (i.e., it was not applied repeatedly, just once at the beginning). 

**KELT SURVIVAL IN BURN-IN**

In the original burnin.R file, initialized kelt survival was subjected to outmigrating mortality at dams. This would have the effect of reducing the number of reconditioning fish available to start the "natural" population in the burnin-in period. 

**HERITABILITY**

Changed from a uniform distribution -2.5 to 2.5 to a normal distribution (mean = 0, sd = 1)

**SIZE AND GROWTH OF OFFSPRING**

In the original model, offspring were simply drawn from a distribution of individual parent sizes that had a random, environmental component already applied. In the new version, offspring are drawn from a vector of mid-parent sizes and the random environmental component is applied to each individual offspring. Instead of two sets of identical offspring from a spawning event, there are an infinite number of types of offspring from a spawning event. 
