# Zmix (BETA)
##R package for overfitting Bayesian mixture models with an unknown number of components.

### Before use: currently updating documentation etc, please contact me for info / questions :)  

The code below is an example of a typical analysis:

*Install package:*

    library(devtools)
    library(roxygen2)
    install_github('zoevanhavre/Zmix')
    library(Zmix)
    # NOTE: please move to a personal directory at this point as plots are created automatically with the Process_Output_Zmix function.
    # setwd("FOLDERPATH_FOR_OUTPUT")
    set.seed(1)	
    dat1<-c(rnorm(20, mean=-1), rnorm(70, mean=4), rnorm(10, mean=10))
    set.seed(1)	
    run1 <- Zmix_univ_tempered(dat1,tau=1,iter=2000,k=10) # few iterations for fast example, in practice use more 
    pp_run1<-Process_Output_Zmix(run1,Pred_Reps=1000, Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="Sim 1 Example of Zmix", SaveFileName="Zmix_Run1", Burn=1000)
    #Post Processing output is a  list of 
    #1 final parameter estimates for all models
    pp_run1[[1]]
    #2  model fit values
    pp_run1[[2]]
    # 3 Posterior allocations (list over all models found)
    head(pp_run1[[3]][[1]])
    #head(pp_run1[[3]][[2]]) # When there is a second model, etc
    # 4 Posterior Allocation probabilities for each observation  and each non-empty group (list over all models found)
    head(pp_run1[[4]][[1]])
    #head(pp_run1[[4]][[2]]) # When there is a second model, etc
    #5 Unswitched parameters  (list over all models found)
    head(pp_run1[[5]][[1]][[1]])  # Parameters
    head(pp_run1[[5]][[1]][[2]])  # Allocations (Zs)

*Galaxy Example*

    set.seed(1)
    runGalaxy1 <- Zmix_univ_tempered(Galaxy , tau=1, iter=2000,k=10)
    g1<-Process_Output_Zmix(runGalaxy1,LineUp=1, Pred_Reps=1000,Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="Galaxy, tau=1", SaveFileName="Zmix_Galaxy_tau1", Burn=1000)
    # increase prior variance of mean: 
    run1_t01 <- Zmix_univ_tempered(dat1,  tau=0.01, iter=2000,k=10)
    pp_run1_t01<-Process_Output_Zmix(run1_t01,  Pred_Reps=1000, Zswitch_Sensitivity=0.01, isSim=FALSE, Plot_Title="dat1 with Tau=0.01", SaveFileName="Run1_tau01", Burn=1000)

