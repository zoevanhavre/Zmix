# Zmix
R package for overfitting Bayesian mixture models with an unknown number of components.
Documentation in progress, please contact me for questions etc 

HOW TO USE

```# install package
library(devtools)
library(roxygen2)```

`install_github('zoevanhavre/Zmix')`
`library(Zmix)`
`# NOTE: please move to a personal directory at this point as plots are created automatically with the Process_Output_Zmix function.`
`# setwd("FOLDERPATH_FOR_OUTPUT")`
 
# Galaxy
set.seed(1)
runGalaxy2 <- Zmix_univ_tempered(Galaxy , tau=0.01, iter=10000, k=10)
g2<-Process_Output_Zmix(runGalaxy2,LineUp=1, Pred_Reps=1000, Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="Galaxy, tau=0.01", SaveFileName="Galaxy_2", Burn=5000)

#Simulation n1
set.seed(1)	
dat1<-c(rnorm(50, mean=0), rnorm(100, mean=3), rnorm(20, mean=7), rnorm(10, mean=12))
run1_t01 <- Zmix_univ_tempered(dat1,  tau=0.01, iter=10000, 	k=10)
pp_run1_t01<-Process_Output_Zmix(run1_t01,  Pred_Reps=1000, Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="dat1 with Tau=0.01", SaveFileName="Run1_tau01", Burn=5000)
