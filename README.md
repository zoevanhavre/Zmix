# Zmix
##R package for overfitting Bayesian mixture models with an unknown number of components.
###Documentation in progress, please contact me for questions etc 

HOW TO USE

*Install package:*

    library(devtools)
    library(roxygen2)
    install_github('zoevanhavre/Zmix')
    library(Zmix)
    # NOTE: please move to a personal directory at this point as plots are created automatically with the Process_Output_Zmix function.
    # setwd("FOLDERPATH_FOR_OUTPUT")
    set.seed(1)	
    dat1<-c(rnorm(30, mean=3), rnorm(70, mean=6))
    set.seed(1)	
    dat1<-c(rnorm(30, mean=3), rnorm(70, mean=6))
    run1 <- Zmix_univ_tempered(dat1,tau=1,iter=10000,k=10)
    pp_run1_t01<-Process_Output_Zmix(run1,Pred_Reps=1000, Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="dat1 with Tau=0.01", SaveFileName="Zmix_Run1", Burn=5000)

*Galaxy Example*
    set.seed(1)
    runGalaxy1 <- Zmix_univ_tempered(Galaxy , tau=1, iter=10000,k=10)
    g1<-Process_Output_Zmix(runGalaxy1,LineUp=1, Pred_Reps=1000,Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="Galaxy, tau=1", SaveFileName="Zmix_Galaxy_tau1", Burn=5000)

    # increase prior variance of mean: 
    run1_t01 <- Zmix_univ_tempered(dat1,  tau=0.01, iter=10000,k=10)
    pp_run1_t01<-Process_Output_Zmix(run1_t01,  Pred_Reps=1000, Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="dat1 with Tau=0.01", SaveFileName="Run1_tau01", Burn=5000)

