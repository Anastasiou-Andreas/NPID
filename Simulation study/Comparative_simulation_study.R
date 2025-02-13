#===============================================================
# Author: Andreas Anastasiou
# Project: Non-parametric multiple change-point detection

# Description: The code below carries out the comparative simulation study as
# explained in the paper. First, the user needs to run all the functions in the
# main.R script

#===============================================================
rm(list = ls())   # Clear environment

library(changepoint)
library(changepoint.np) 
library(cpm)
library(ecp)
library(wbs)

source(main.R)

#=================================
# Some utility functions
#=================================

## Based on a given set of estimated change-points, this function calculates the
## fit in the scenario of changes in the mean.
mean.from.cpt <- function(x, cpt) {
  
  n <- length(x)
  len.cpt <- length(cpt)
  if (len.cpt) cpt <- sort(cpt)
  beg <- endd <- rep(0, len.cpt+1)
  beg[1] <- 1
  endd[len.cpt+1] <- n
  if (len.cpt) {
    beg[2:(len.cpt+1)] <- cpt+1
    endd[1:len.cpt] <- cpt
  }
  means <- rep(0, len.cpt+1)
  for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])
  rep(means, endd-beg+1)
}

## The main function that returns the results for our method, NPID. This will be used
## in the wrapper function for the simulation study
rev.sim.NPID <- function(x, thr_ic_KS_rs = 1.7) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  print("L_inf with IC")
  z <- model_selection_L_inf(x, th_const = thr_ic_KS_rs, rescale = TRUE)
  cpt.z = z$cpt_ic
  
  Linf_RS_IC <- new("cpt.est")
  if(is.na(cpt.z)){Linf_RS_IC@cpt = 0
  Linf_RS_IC@nocpt = 0}
  else{Linf_RS_IC@cpt <- cpt.z
  Linf_RS_IC@nocpt <- length(cpt.z)}
  Linf_RS_IC@time <- system.time(model_selection_L_inf(x, th_const = thr_ic_KS_rs, rescale = TRUE))[[3]]
  
  list(Linf_RS_IC = Linf_RS_IC)
}

## This is the wrapper function for the results of our method in the simulation
## study for changes in the mean
rev.sim.study.NPID.mean <- function(signal, sigma = 1, no.of.cpt, true.cpt=NULL, m = 100, seed = NULL, resc = FALSE) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  n <- length(signal)
  ns <- max(diff(c(0,true.cpt, n)))
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    est <- rev.sim.NPID(x)
    
    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt - no.of.cpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@diff <- abs(matrix(est$Linf_RS_IC@cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=F))
    Linf_RS_IC@dh[i] <- max(apply(Linf_RS_IC@diff,1,min),apply(Linf_RS_IC@diff,2,min))/ns
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time
    
    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}

## The function that returns the results for the competitors. This will be used
## in a wrapper function for the simulation study
rev.sim.comp <- function(x) {
  
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric",time="numeric"), prototype(cpt=numeric(0), nocpt=0, mean=numeric(0),time=numeric(0)))
  
  print("pelt.np")
  
  lx <- length(x)
  
  z <- cpt.np(x, method="PELT", nquantiles = 10)
  temp <- cpts(z)
  pelt.np <- new("cpt.est")
  if (length(temp) == 0){pelt.np@cpt = 0
  pelt.np@nocpt = 0}else{
    pelt.np@cpt <- as.numeric(temp)
    pelt.np@nocpt <- length(pelt.np@cpt)}
  pelt.np@time <- system.time(cpt.np(x, method="PELT",nquantiles = 10))[[3]]
  
  print("CPMKS")
  z <- processStream(x,cpmType = "Kolmogorov-Smirnov",ARL0 = 500)
  CPM_points <- new("cpt.est")
  CPM_points@cpt <- as.numeric(z$changePoints)
  if (length(CPM_points@cpt) == 0){CPM_points@cpt <- 0}
  CPM_points@nocpt <- ifelse(CPM_points@cpt==0,0,length(CPM_points@cpt))
  CPM_points@time <- system.time(processStream(x,cpmType = "Kolmogorov-Smirnov"))[[3]]
  
  print("CPMKS2")
  z <- processStream(x,cpmType = "Kolmogorov-Smirnov",ARL0 = 1000)
  CPM_points2 <- new("cpt.est")
  CPM_points2@cpt <- as.numeric(z$changePoints)
  if (length(CPM_points2@cpt) == 0){CPM_points2@cpt <- 0}
  CPM_points2@nocpt <- ifelse(CPM_points2@cpt==0,0,length(CPM_points2@cpt))
  CPM_points2@time <- system.time(processStream(x,cpmType = "Kolmogorov-Smirnov",ARL0 = 1000))[[3]]
  
  print("CPMcvm")
  z <- processStream(x,cpmType = "Cramer-von-Mises",ARL0 = 500)
  CPM_points3 <- new("cpt.est")
  CPM_points3@cpt <- as.numeric(z$changePoints)
  if (length(CPM_points3@cpt) == 0){CPM_points3@cpt <- 0}
  CPM_points3@nocpt <- ifelse(CPM_points3@cpt==0,0,length(CPM_points3@cpt))
  CPM_points3@time <- system.time(processStream(x,cpmType = "Cramer-von-Mises"))[[3]]
  
  print("CPMcvm2")
  z <- processStream(x,cpmType = "Cramer-von-Mises",ARL0 = 1000)
  CPM_points4 <- new("cpt.est")
  CPM_points4@cpt <- as.numeric(z$changePoints)
  if (length(CPM_points4@cpt) == 0){CPM_points4@cpt <- 0}
  CPM_points4@nocpt <- ifelse(CPM_points4@cpt==0,0,length(CPM_points4@cpt))
  CPM_points4@time <- system.time(processStream(x,cpmType = "Cramer-von-Mises",ARL0 = 1000))[[3]]
  
  print("ecp")
  z <- e.cp3o_delta(Z=matrix(x, length(x), 1), K=4, delta=29, alpha=1, verbose=FALSE)
  ecp <- new("cpt.est")
  ecp@cpt <- z$estimates
  ecp@nocpt <- z$number
  ecp@time <- z$time
  
  print("ecpKS")
  z <- ks.cp3o_delta(Z=matrix(x, length(x), 1), K=4, minsize=30, verbose=FALSE)
  ecpKS <- new("cpt.est")
  ecpKS@cpt <- z$estimates
  ecpKS@nocpt <- z$number
  ecpKS@time <- z$time
  
  print("NBS")
  s= 1
  e =  lx
  flag = 0 
  
  S =  NULL
  
  m = floor(lx/2)
  y_new =  cbind(x[2*(1:m)-1],x[2*(1:m)])
  
  S = NBS_full(as.matrix(y_new[,1]),as.matrix(y_new[,2]),gam=2,rep(1,m))
  S = 2*S
  
  S =  sort(unique(S))
  
  nbs <- new("cpt.est")
  if (length(S) == 0){nbs@cpt = 0
  nbs@nocpt = 0}else{
    nbs@cpt <- as.numeric(S)
    nbs@nocpt <- length(nbs@cpt)}
  nbs@time <- system.time(NBS_full(as.matrix(y_new[,1]),as.matrix(y_new[,2]),gam=2,rep(1,m)))[[3]]
  
  print("NWBS")
  M =   120
  alpha =  sample.int(size =M  , n = m,replace = TRUE)
  beta =   sample.int(size =M  , n = m,replace = TRUE)#alpha + floor((T- alpha)*runif(M))
  #
  for(j in 1:M)
  {
    aux =  alpha[j]
    aux2 =  beta[j]
    #
    alpha[j] = min(aux,aux2)
    beta[j] = max(aux,aux2)
  }
  
  S2 = NWBS_full(as.matrix(y_new[,1]),as.matrix(y_new[,2]),gam=2,rep(1,m),alpha,beta)
  S2 = 2*S2 
  S2 =  sort(unique(S2))
  
  nwbs <- new("cpt.est")
  if (length(S2) == 0){nwbs@cpt = 0
  nwbs@nocpt = 0}else{
    nwbs@cpt <- as.numeric(S2)
    nwbs@nocpt <- length(nwbs@cpt)}
  nwbs@time <- system.time(NWBS_full(as.matrix(y_new[,1]),as.matrix(y_new[,2]),gam=2,rep(1,m),alpha,beta))[[3]]
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4,ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}

## This is the wrapper function for the results of the competitors in the simulation
## study for changes in the mean
rev.sim.study.comp_mean <- function(signal, sigma = 1, no.of.cpt, true.cpt=NULL, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  n <- length(signal)
  ns <- max(diff(c(0,true.cpt, n)))
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    est <- rev.sim.comp(x)
    
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt - no.of.cpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@diff <- abs(matrix(est$CPM_points3@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=F))
    CPM_points3@dh[i] <- max(apply(CPM_points3@diff,1,min),apply(CPM_points3@diff,2,min))/ns
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt - no.of.cpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@diff <- abs(matrix(est$CPM_points4@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=F))
    CPM_points4@dh[i] <- max(apply(CPM_points4@diff,1,min),apply(CPM_points4@diff,2,min))/ns
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt - no.of.cpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@diff <- abs(matrix(est$ecp@cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=F))
    ecp@dh[i] <- max(apply(ecp@diff,1,min),apply(ecp@diff,2,min))/ns
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt - no.of.cpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@diff <- abs(matrix(est$ecpKS@cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=F))
    ecpKS@dh[i] <- max(apply(ecpKS@diff,1,min),apply(ecpKS@diff,2,min))/ns
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt - no.of.cpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@diff <- abs(matrix(est$nbs@cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=F))
    nbs@dh[i] <- max(apply(nbs@diff,1,min),apply(nbs@diff,2,min))/ns
    nbs@time[i] <- est$nbs@time
    
    nwbs@dnc[i] <- est$nwbs@nocpt - no.of.cpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@diff <- abs(matrix(est$nwbs@cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=F))
    nwbs@dh[i] <- max(apply(nwbs@diff,1,min),apply(nwbs@diff,2,min))/ns
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}



#=================================
# Results for no changes
#=================================

## The main wrapper function for our method
rev.sim.study.no_cpt <- function(lx = NULL, m = 20, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- rnorm(lx)
    
    est <- rev.sim.NPID(x)

    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time
    
    
    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}

SIM_nocpt_500 <- rev.sim.study.no_cpt(lx = 500, m = 100, seed = 12)

## The main wrapper function for the competitors

rev.sim.study.no_cpt_comp <- function(lx, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- rnorm(lx)
    
    est <- rev.sim.comp(x)
    
    CPM_points@dnc[i] <- est$CPM_points@nocpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@time[i] <- est$nbs@time
    
    nwbs@dnc[i] <- est$nwbs@nocpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}

SIM_nocpt_500_comp <- rev.sim.study.no_cpt_comp(lx = 500, m = 100, seed = 12)


#=================================
# Results for changes in the mean
#=================================

## Signals used
mean_test1 <- c(rep(0,100),rep(1,100))

multi_mean_test1 <- c(rep(0,100),rep(1,100), rep(-0.2,100), rep(-1.3,100))

multi_mean_test2 <- rep(c(rep(0,80),rep(2,80), rep(0,80), rep(2,80)), 5)

## Results for our method
SIM_mean_1 <- rev.sim.study.NPID.mean(mean_test1, sigma = 1, seed = 12, m = 100, true.cpt = c(100), no.of.cpt = 1)
SIM_multi_mean_1 <- rev.sim.study.NPID.mean(multi_mean_test1, sigma = 1, seed = 12, m = 100, true.cpt = c(100, 200, 300), no.of.cpt = 3)
SIM_multi_mean_2 <- rev.sim.study.NPID.mean(multi_mean_test2, sigma = 1, seed = 12, m = 100, true.cpt = seq(80,1520,80), no.of.cpt = 19)

## Results for the competing methods
SIM_mean_3_comp <- rev.sim.study.comp_mean(mean_test1, sigma = 1, seed = 12, m = 100, true.cpt = c(100), no.of.cpt = 1)
SIM_multi_mean_1_comp <- rev.sim.study.comp_mean(multi_mean_test1, sigma = 1, seed = 12, m = 100, true.cpt = c(100, 200, 300), no.of.cpt = 3)
SIM_multi_mean_2_comp <- rev.sim.study.comp_mean(multi_mean_test2, sigma = 1, seed = 12, m = 100, true.cpt = seq(80,1520,80), no.of.cpt = 19)


#====================================
# Results for changes in the variance
#====================================

## Results for our method
rev.sim.study.NPID.var1 <- function(lx = NULL, true.cpt = NULL, m = 20, seed = NULL, sv = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(250), rnorm(250, sd = sv))
    
    est <- rev.sim.NPID(x)
    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt - no.of.cpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@diff <- abs(matrix(est$Linf_RS_IC@cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=F))
    Linf_RS_IC@dh[i] <- max(apply(Linf_RS_IC@diff,1,min),apply(Linf_RS_IC@diff,2,min))/ns
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time
    
    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}

SIM_var_1 <- rev.sim.study.NPID.var1(lx = 500, true.cpt = c(250), m = 100, seed = 12, sv = 2)


rev.sim.study.NPID.var3 <- function(lx = NULL, true.cpt = NULL, m = 20, seed = NULL, sv = NULL)) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(true.cpt[1], sd = sv[1]), rnorm((true.cpt[2] - true.cpt[1]), sd = sv[2]), rnorm((true.cpt[3] - true.cpt[2]), sd = sv[3]), rnorm((lx - true.cpt[3]), sd = sv[4]))
    
    est <- rev.sim.NPID(x)
    
    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt - no.of.cpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@diff <- abs(matrix(est$Linf_RS_IC@cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=F))
    Linf_RS_IC@dh[i] <- max(apply(Linf_RS_IC@diff,1,min),apply(Linf_RS_IC@diff,2,min))/ns
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time
    
    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}


SIM_multi_var_1 <- rev.sim.study.NPID.var3(lx = 600, true.cpt = c(150,350,500), m = 100, seed = 12, sv = c(1,3,1.2,sqrt(0.1)))


rev.sim.study.NPID.var5 <- function(lx = NULL, true.cpt = NULL, m = 20, seed = NULL, sv = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(true.cpt[1], sd = sv[1]), rnorm((true.cpt[2] - true.cpt[1]), sd = sv[2]), rnorm((true.cpt[3] - true.cpt[2]), sd = sv[3]), rnorm((true.cpt[4] - true.cpt[3]), sd = sv[4]), rnorm((true.cpt[5] - true.cpt[4]), sd = sv[5]), rnorm((lx - true.cpt[5]), sd = sv[6]))
    
    est <- rev.sim.NPID(x)
    
    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt - no.of.cpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@diff <- abs(matrix(est$Linf_RS_IC@cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=F))
    Linf_RS_IC@dh[i] <- max(apply(Linf_RS_IC@diff,1,min),apply(Linf_RS_IC@diff,2,min))/ns
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time
    
    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}

SIM_multi_var_2 <- rev.sim.study.NPID.var5(lx = 1000, true.cpt = c(200,350,550,700,900), m = 100, seed = 12, sv = c(3,sqrt(2),sqrt(0.3),2,sqrt(10),sqrt(2)))


## Results for competitors
rev.sim.study.comp.var1 <- function(lx, sv = NULL, true.cpt=NULL, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(150), rnorm(150, sd = sv)) 
    
    est <- rev.sim.comp(x)
    
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt - no.of.cpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@diff <- abs(matrix(est$CPM_points3@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=F))
    CPM_points3@dh[i] <- max(apply(CPM_points3@diff,1,min),apply(CPM_points3@diff,2,min))/ns
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt - no.of.cpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@diff <- abs(matrix(est$CPM_points4@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=F))
    CPM_points4@dh[i] <- max(apply(CPM_points4@diff,1,min),apply(CPM_points4@diff,2,min))/ns
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt - no.of.cpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@diff <- abs(matrix(est$ecp@cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=F))
    ecp@dh[i] <- max(apply(ecp@diff,1,min),apply(ecp@diff,2,min))/ns
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt - no.of.cpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@diff <- abs(matrix(est$ecpKS@cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=F))
    ecpKS@dh[i] <- max(apply(ecpKS@diff,1,min),apply(ecpKS@diff,2,min))/ns
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt - no.of.cpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@diff <- abs(matrix(est$nbs@cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=F))
    nbs@dh[i] <- max(apply(nbs@diff,1,min),apply(nbs@diff,2,min))/ns
    nbs@time[i] <- est$nbs@time
    
    nwbs@dnc[i] <- est$nwbs@nocpt - no.of.cpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@diff <- abs(matrix(est$nwbs@cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=F))
    nwbs@dh[i] <- max(apply(nwbs@diff,1,min),apply(nwbs@diff,2,min))/ns
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}

SIM_var_3_comp <- rev.sim.study.comp.var1(lx = 300, true.cpt = c(150), m = 10, seed = 12, sv = 2)

rev.sim.study.comp.var3<- function(lx, sv = NULL, true.cpt=NULL, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(true.cpt[1], sd = sv[1]), rnorm((true.cpt[2] - true.cpt[1]), sd = sv[2]), rnorm((true.cpt[3] - true.cpt[2]), sd = sv[3]), rnorm((lx - true.cpt[3]), sd = sv[4]))
    
    est <- rev.sim.comp(x)
    
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt - no.of.cpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@diff <- abs(matrix(est$CPM_points3@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=F))
    CPM_points3@dh[i] <- max(apply(CPM_points3@diff,1,min),apply(CPM_points3@diff,2,min))/ns
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt - no.of.cpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@diff <- abs(matrix(est$CPM_points4@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=F))
    CPM_points4@dh[i] <- max(apply(CPM_points4@diff,1,min),apply(CPM_points4@diff,2,min))/ns
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt - no.of.cpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@diff <- abs(matrix(est$ecp@cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=F))
    ecp@dh[i] <- max(apply(ecp@diff,1,min),apply(ecp@diff,2,min))/ns
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt - no.of.cpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@diff <- abs(matrix(est$ecpKS@cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=F))
    ecpKS@dh[i] <- max(apply(ecpKS@diff,1,min),apply(ecpKS@diff,2,min))/ns
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt - no.of.cpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@diff <- abs(matrix(est$nbs@cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=F))
    nbs@dh[i] <- max(apply(nbs@diff,1,min),apply(nbs@diff,2,min))/ns
    nbs@time[i] <- est$nbs@time
    
    nwbs@dnc[i] <- est$nwbs@nocpt - no.of.cpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@diff <- abs(matrix(est$nwbs@cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=F))
    nwbs@dh[i] <- max(apply(nwbs@diff,1,min),apply(nwbs@diff,2,min))/ns
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}
SIM_multi_var_1_comp <- rev.sim.study.comp.var3(lx = 600, true.cpt = c(150,350,500), m = 100, seed = 12, sv = c(1,3,1.2,sqrt(0.1)))

rev.sim.study.comp.var5 <- function(lx, sv = NULL, true.cpt=NULL, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(true.cpt[1], sd = sv[1]), rnorm((true.cpt[2] - true.cpt[1]), sd = sv[2]), rnorm((true.cpt[3] - true.cpt[2]), sd = sv[3]),rnorm((true.cpt[4] - true.cpt[3]), sd = sv[4]), rnorm((true.cpt[5] - true.cpt[4]), sd = sv[5]), rnorm((lx - true.cpt[5]), sd = sv[6]))
    
    est <- rev.sim.comp(x)
    
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt - no.of.cpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@diff <- abs(matrix(est$CPM_points3@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=F))
    CPM_points3@dh[i] <- max(apply(CPM_points3@diff,1,min),apply(CPM_points3@diff,2,min))/ns
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt - no.of.cpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@diff <- abs(matrix(est$CPM_points4@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=F))
    CPM_points4@dh[i] <- max(apply(CPM_points4@diff,1,min),apply(CPM_points4@diff,2,min))/ns
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt - no.of.cpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@diff <- abs(matrix(est$ecp@cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=F))
    ecp@dh[i] <- max(apply(ecp@diff,1,min),apply(ecp@diff,2,min))/ns
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt - no.of.cpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@diff <- abs(matrix(est$ecpKS@cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=F))
    ecpKS@dh[i] <- max(apply(ecpKS@diff,1,min),apply(ecpKS@diff,2,min))/ns
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt - no.of.cpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@diff <- abs(matrix(est$nbs@cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=F))
    nbs@dh[i] <- max(apply(nbs@diff,1,min),apply(nbs@diff,2,min))/ns
    nbs@time[i] <- est$nbs@time
     
    nwbs@dnc[i] <- est$nwbs@nocpt - no.of.cpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@diff <- abs(matrix(est$nwbs@cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=F))
    nwbs@dh[i] <- max(apply(nwbs@diff,1,min),apply(nwbs@diff,2,min))/ns
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}

SIM_multi_var_2_comp <- rev.sim.study.comp.var5(lx = 1000, true.cpt = c(200,350,550,700,900), m = 100, seed = 12, sv = c(3,sqrt(2),sqrt(0.3),2,sqrt(10),sqrt(2)))


#===========================================
# Results for general distributional changes
#===========================================

## Results for our method

rev.sim.study.NPID.dist3 <- function(lx = NULL, true.cpt = NULL, m = 20, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(true.cpt[1]), rchisq(true.cpt[2] - true.cpt[1], 1), rt(true.cpt[3] - true.cpt[2], 3), rnorm(lx - true.cpt[3],1,1))
    
    est <- rev.sim.NPID(x)

    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt - no.of.cpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@diff <- abs(matrix(est$Linf_RS_IC@cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=F))
    Linf_RS_IC@dh[i] <- max(apply(Linf_RS_IC@diff,1,min),apply(Linf_RS_IC@diff,2,min))/ns
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time

    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}

SIM_distr_1 <- rev.sim.study.NPID.dist3(500, true.cpt = c(100,250,350), m= 100, seed = 12)


rev.sim.study.NPID.dist3_2 <- function(lx = NULL, true.cpt = NULL, m = 20, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  Linf_RS_IC <- new("est.eval")
  
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rgamma(true.cpt[1],1,1), rchisq(true.cpt[2] - true.cpt[1], 3), rnorm(true.cpt[3] - true.cpt[2],0.5,1), rt(lx - true.cpt[3],5,0))
    
    est <- rev.sim.NPID(x)

    Linf_RS_IC@dnc[i] <- est$Linf_RS_IC@nocpt - no.of.cpt
    Linf_RS_IC@cpt[[i]] <- est$Linf_RS_IC@cpt
    Linf_RS_IC@diff <- abs(matrix(est$Linf_RS_IC@cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$Linf_RS_IC@cpt),byr=F))
    Linf_RS_IC@dh[i] <- max(apply(Linf_RS_IC@diff,1,min),apply(Linf_RS_IC@diff,2,min))/ns
    Linf_RS_IC@time[i] <- est$Linf_RS_IC@time
    
    gc()
  }
  
  list(Linf_RS_IC = Linf_RS_IC)
}

SIM_distr_2 <- rev.sim.study.NPID.dist3_2(1000, true.cpt = c(200,500,750), m= 100, seed = 12)


## Results for the competitors

rev.sim.study.comp.dist3 <- function(lx = NULL, true.cpt = NULL, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rnorm(true.cpt[1]), rchisq(true.cpt[2] - true.cpt[1], 1), rt(true.cpt[3] - true.cpt[2], 3), rnorm(lx - true.cpt[3],1,1))
    
    est <- rev.sim.comp(x)
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt - no.of.cpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@diff <- abs(matrix(est$CPM_points3@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=F))
    CPM_points3@dh[i] <- max(apply(CPM_points3@diff,1,min),apply(CPM_points3@diff,2,min))/ns
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt - no.of.cpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@diff <- abs(matrix(est$CPM_points4@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=F))
    CPM_points4@dh[i] <- max(apply(CPM_points4@diff,1,min),apply(CPM_points4@diff,2,min))/ns
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt - no.of.cpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@diff <- abs(matrix(est$ecp@cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=F))
    ecp@dh[i] <- max(apply(ecp@diff,1,min),apply(ecp@diff,2,min))/ns
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt - no.of.cpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@diff <- abs(matrix(est$ecpKS@cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=F))
    ecpKS@dh[i] <- max(apply(ecpKS@diff,1,min),apply(ecpKS@diff,2,min))/ns
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt - no.of.cpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@diff <- abs(matrix(est$nbs@cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=F))
    nbs@dh[i] <- max(apply(nbs@diff,1,min),apply(nbs@diff,2,min))/ns
    nbs@time[i] <- est$nbs@time
    
    nwbs@dnc[i] <- est$nwbs@nocpt - no.of.cpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@diff <- abs(matrix(est$nwbs@cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=F))
    nwbs@dh[i] <- max(apply(nwbs@diff,1,min),apply(nwbs@diff,2,min))/ns
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}

SIM_distr_1_comp <- rev.sim.study.comp.dist3(500, true.cpt = c(100,250,350), m= 100, seed = 12)


rev.sim.study.comp.dist3_2 <- function(lx = NULL, true.cpt = NULL, m = 100, seed = NULL) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric", mse="numeric", time= "numeric"), prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  pelt.np <- new("est.eval")
  CPM_points <- new("est.eval")
  CPM_points2 <- new("est.eval")
  CPM_points3 <- new("est.eval")
  CPM_points4 <- new("est.eval")
  ecp <- new("est.eval")
  ecpKS <- new("est.eval")
  nbs <- new("est.eval")
  nwbs <- new("est.eval")
  
  no.of.cpt <- length(true.cpt)
  ns <- max(diff(c(0,true.cpt,lx)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- c(rgamma(true.cpt[1],1,1), rchisq(true.cpt[2] - true.cpt[1], 3), rnorm(true.cpt[3] - true.cpt[2],0.5,1), rt(lx - true.cpt[3],5,0))
    
    est <- rev.sim.comp(x)
    CPM_points@dnc[i] <- est$CPM_points@nocpt - no.of.cpt
    CPM_points@cpt[[i]] <- est$CPM_points@cpt
    CPM_points@diff <- abs(matrix(est$CPM_points@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points@cpt),byr=F))
    CPM_points@dh[i] <- max(apply(CPM_points@diff,1,min),apply(CPM_points@diff,2,min))/ns
    CPM_points@time[i] <- est$CPM_points@time
    
    CPM_points2@dnc[i] <- est$CPM_points2@nocpt - no.of.cpt
    CPM_points2@cpt[[i]] <- est$CPM_points2@cpt
    CPM_points2@diff <- abs(matrix(est$CPM_points2@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points2@cpt),byr=F))
    CPM_points2@dh[i] <- max(apply(CPM_points2@diff,1,min),apply(CPM_points2@diff,2,min))/ns
    CPM_points2@time[i] <- est$CPM_points2@time
    
    CPM_points3@dnc[i] <- est$CPM_points3@nocpt - no.of.cpt
    CPM_points3@cpt[[i]] <- est$CPM_points3@cpt
    CPM_points3@diff <- abs(matrix(est$CPM_points3@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points3@cpt),byr=F))
    CPM_points3@dh[i] <- max(apply(CPM_points3@diff,1,min),apply(CPM_points3@diff,2,min))/ns
    CPM_points3@time[i] <- est$CPM_points3@time
    
    CPM_points4@dnc[i] <- est$CPM_points4@nocpt - no.of.cpt
    CPM_points4@cpt[[i]] <- est$CPM_points4@cpt
    CPM_points4@diff <- abs(matrix(est$CPM_points4@cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$CPM_points4@cpt),byr=F))
    CPM_points4@dh[i] <- max(apply(CPM_points4@diff,1,min),apply(CPM_points4@diff,2,min))/ns
    CPM_points4@time[i] <- est$CPM_points4@time
    
    pelt.np@dnc[i] <- est$pelt.np@nocpt - no.of.cpt
    pelt.np@cpt[[i]] <- est$pelt.np@cpt
    pelt.np@diff <- abs(matrix(est$pelt.np@cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt.np@cpt),byr=F))
    pelt.np@dh[i] <- max(apply(pelt.np@diff,1,min),apply(pelt.np@diff,2,min))/ns
    pelt.np@time[i] <- est$pelt.np@time
    
    ecp@dnc[i] <- est$ecp@nocpt - no.of.cpt
    ecp@cpt[[i]] <- est$ecp@cpt
    ecp@diff <- abs(matrix(est$ecp@cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecp@cpt),byr=F))
    ecp@dh[i] <- max(apply(ecp@diff,1,min),apply(ecp@diff,2,min))/ns
    ecp@time[i] <- est$ecp@time
    
    ecpKS@dnc[i] <- est$ecpKS@nocpt - no.of.cpt
    ecpKS@cpt[[i]] <- est$ecpKS@cpt
    ecpKS@diff <- abs(matrix(est$ecpKS@cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ecpKS@cpt),byr=F))
    ecpKS@dh[i] <- max(apply(ecpKS@diff,1,min),apply(ecpKS@diff,2,min))/ns
    ecpKS@time[i] <- est$ecpKS@time
    
    nbs@dnc[i] <- est$nbs@nocpt - no.of.cpt
    nbs@cpt[[i]] <- est$nbs@cpt
    nbs@diff <- abs(matrix(est$nbs@cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nbs@cpt),byr=F))
    nbs@dh[i] <- max(apply(nbs@diff,1,min),apply(nbs@diff,2,min))/ns
    nbs@time[i] <- est$nbs@time
     
    nwbs@dnc[i] <- est$nwbs@nocpt - no.of.cpt
    nwbs@cpt[[i]] <- est$nwbs@cpt
    nwbs@diff <- abs(matrix(est$nwbs@cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$nwbs@cpt),byr=F))
    nwbs@dh[i] <- max(apply(nwbs@diff,1,min),apply(nwbs@diff,2,min))/ns
    nwbs@time[i] <- est$nwbs@time
    
    gc()
  }
  
  list(pelt.np=pelt.np,CPM_points=CPM_points,CPM_points2=CPM_points2,
       CPM_points3=CPM_points3, CPM_points4=CPM_points4, ecp = ecp, ecpKS = ecpKS,
       nbs = nbs, nwbs = nwbs)
}

SIM_distr_2_comp <- rev.sim.study.comp.dist3_2(1000, true.cpt = c(200,500,750), m= 100, seed = 12)

