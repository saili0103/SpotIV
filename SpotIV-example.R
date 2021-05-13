
##estimation and construct confidence intervals for CATE(d1,d2|z0)
set.seed(123)
library(MASS)
library(mvtnorm)
library(dr)
library(orthoDr)
library(np)
library(foreach)   #on windows machine, you may need different packages for parallel computing,
library(doParallel)  #such as "library(parallel);detectCores()"
library(caret)
registerDoParallel(4)

source("SpotIV-main.R")

######experiments correspond to settings i, ii, iii & iv
beta0 = 0.25
rho.z = 0
rho.err= 0.25
d1 = -1 #for ASF(d1,d2|z0)
d2 = 2
J=7
z0<-c(rep(0,J-1), 0.1)
normal.Z=F #normal IV measurements of uniform IV measurements
Nround <- 500 
mm=2 #1 & 2 are binary outcomes; 3& 4 are continuous outcomes
set.seed(123)
for(n in c(1000, 500)){
  for(c.gam in c(0.8, 0.6, 0.4)){
    cov=0
    cov.oracle<- 0
    re.semi.mat = matrix(NA, ncol=5, nrow = Nround)  ##### result matrix
    ########
    maj.fail.temp=0
    for(round in 1:Nround){
      dat0 <- gen.bidata(n=n, c.gam = c.gam, beta0=beta0,
                         sig.v.sq = 1, sig.u.sq = 1, rho.err = 0.25, J=J, 
                         s0=5, d1=d1, d2=d2, z0=z0, c.alp=0.8, mm=mm, normal.Z=normal.Z)
      Y<-as.numeric(dat0$Y)
      D <- dat0$D
      Z <- dat0$Z
      sim.re <- main.semi(Y, D, Z,  bootstrap = T, bs.Niter = 40,
                          d1 = d1, d2 = d2, z0 = z0, M=3, M.est=T)
      maj.fail.temp=maj.fail.temp + sim.re$Maj.fail #number of times not passing the majority rule testing
      oracle.re<-main.semi(Y, D, Z,  bootstrap = T, bs.Niter = 40, V=1:5,
                           d1 = d1, d2 = d2, z0 = z0, M=3, M.est=T)
      re.semi.mat[round,] <- as.numeric(c(dat0$cace0, sim.re$cace.hat, sim.re$cace.sd, oracle.re$cace.hat, oracle.re$cace.sd))
      cov<- cov+ (abs(re.semi.mat[round,1]-re.semi.mat[round,2])<=1.96* re.semi.mat[round,3])#oracle coverage
      cat('round=', round, re.semi.mat[round,],cov/round,'\n')
    }#round
    sum.tsth.cace <- sum.result(re.semi.mat[,1:3])
    sum.oracle.cace<- sum.result(re.semi.mat[,c(1,4,5)])
    result.semi <- c(Nround, n,  c.gam, 
                     sum.tsth.cace$mae, sum.tsth.cace$ci.cov.re, sum.tsth.cace$ci.len.re, maj.fail.temp, 
                     sum.oracle.cace$mae,sum.oracle.cace$ci.cov.re, sum.oracle.cace$ci.len.re)#inference results
    cat(result.semi,'\n')
  }##c.gam
}

####experiments for majority rule violation#############
#data generation
gen.bidata.violate <- function(n=1000, c.gam=0.8, beta0, sig.v.sq=1, sig.u.sq=1, rho.err=0.25, J=7, s0=5, 
                               d1=NULL, d2=NULL, z0=NULL,c.alp=0.4, normal.Z=T, mm=1, rn=2000){
  u1<-function(z,eta,n,rho.err, v){
    as.numeric( z%*% eta)+ v*rho.err + rnorm(n) 
  }
  if(normal.Z){
    Z <- matrix(rnorm(n * J, 0, 1) , ncol = J, nrow = n)
  }else{
    Z <- matrix(runif(n * J, -1.73, 1.73) , ncol = J, nrow = n)
  }
  colnames(Z) <- paste("Z", 1:J, sep = '')
  if(mm==1){
    gam <- c(rep(c.gam, 4), rep(-c.gam, J - 4))
    pi0<-c(0,rep(c.alp,J-1))
  }else if(mm==2){
    gam <- c(rep(c.gam, 3), rep(-c.gam, J - 3))
    pi0 <- gam*runif(J,-1,1)
  }
  
  eta<-pi0/2
  v.vec<-rnorm(n)
  D = Z %*% gam + v.vec
  
  u.vec <-u1(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
  kappa <- pi0-eta
  Y<-  sapply(plogis(D * beta0 +Z%*%kappa + u.vec), 
              function(x)  rbinom(1, 1, x))
  u1.r<-u1(z=z0, eta=eta, n=rn, rho.err=rho.err, v=rnorm(rn)) 
  cace0 <- mean(plogis(as.numeric(d1 * beta0 + z0 %*% kappa)+ u1.r )) - 
    mean( plogis(as.numeric(d2 * beta0 + z0 %*% kappa) + u1.r))
  asf0= mean(plogis(as.numeric(d2 * beta0 + z0 %*% kappa)+ u1.r ))
  
  
  list(
    Z = Z, D = D, Y = Y,
    V0 = 1:s0, 
    eta=eta,
    kappa=kappa,
    asf0=asf0,
    cace0 = cace0,
    z0 = z0
  )
}



beta0 = 0.25
rho.z = 0
rho.err= 0.25
d1 = -1
d2 = 2
J=7
z0<-c(rep(0,J-1), 0.1)
normal.Z=F 
Nround <- 500 
mm=1 #mm=1 majority of the IVs fail and can be tested; 1 iv is valid
#mm=2 majority of the IVs fail and invalid effects/gamma ~ U[-1,1]
set.seed(123)
for(n in c(1000,500)){
  for(c.gam in c(0.4, 0.6, 0.8)){
    cov=0
    cov.oracle<- 0
    re.semi.mat = matrix(NA, ncol=5, nrow = Nround)  ##### result matrix
    ########
    maj.fail.temp=0
    for(round in 1:Nround){
      #generate data which violate the majority rule
      dat0 <- gen.bidata.violate(n=n, c.gam = c.gam, beta0=beta0,
                                 sig.v.sq = 1, sig.u.sq = 1, rho.err = 0.25, J=J,
                                 d1=d1, d2=d2, z0=z0, c.alp=0.8, s0=1,mm=mm, normal.Z=normal.Z)
      Y<-as.numeric(dat0$Y)
      D <- dat0$D
      Z <- dat0$Z
      #spotIV
      sim.re <- main.semi(Y, D, Z,  bootstrap = T, bs.Niter = 40,
                          d1 = d1, d2 = d2, z0 = z0, M=2, M.est=T)
      maj.fail.temp=maj.fail.temp + sim.re$Maj.fail #number of times not passing the majority rule testing
      oracle.re<-list(cace.hat=NA, cace.sd=NA)
      
      re.semi.mat[round,] <- as.numeric(c(dat0$cace0, sim.re$cace.hat, sim.re$cace.sd, oracle.re$cace.hat, oracle.re$cace.sd))
      cov<- cov+ (abs(re.semi.mat[round,1]-re.semi.mat[round,2])<=1.96* re.semi.mat[round,3])
      cat('round=', round, re.semi.mat[round,],cov/round,'\n')
      
    }#round
    sum.tsth.cace <- sum.result(re.semi.mat[,1:3])
    sum.oracle.cace<- sum.result(re.semi.mat[,c(1,4,5)])
    result.semi <- c(Nround, n,  c.gam,  #setting
                     sum.tsth.cace$mae, sum.tsth.cace$ci.cov.re, sum.tsth.cace$ci.len.re, maj.fail.temp,  #results of spotIV
                     sum.oracle.cace$mae,sum.oracle.cace$ci.cov.re, sum.oracle.cace$ci.len.re)#result of oracle method

  }##c.gam
  
}


######Test the InSIDE assumption#######
####InSIDE assumption holds
gen.InSIDE.bidata <- function(n=1000, c.gam=0.8, beta0, sig.v.sq=1, sig.u.sq=1, rho.err=0.25, J=6, s0=5, 
                              d1=NULL, d2=NULL, z0=NULL, mm=1,c.alp=0.4,normal.Z=T, rn=2000){
  u1<-function(z,eta,n,rho.err, v){
    as.numeric( z%*% eta)+ v*rho.err + rnorm(n) 
  }
  u2<-function(z,eta,n,rho.err, v){
    exp(as.numeric(z%*% eta)+ rho.err*v) + runif(n,-1,1)
  }
  rho.z = 0
  if(normal.Z){
    Z <- matrix(rnorm(n * J, 0, 1) , ncol = J, nrow = n)
  }else{
    Z <- matrix(runif(n * J, -1.73, 1.73) , ncol = J, nrow = n)
  }
  
  
  colnames(Z) <- paste("Z", 1:J, sep = '')
  gam <- c(rep(c.gam, floor(J / 2)), rep(-c.gam, J - floor(J / 2)))
  if(s0==J){
    pi0 <- rep(0,J)
  }else if (s0>J/2){
    pi0 <- c(rep(0, s0), c.alp, -c.alp)
  }else{
    pi0 <- c(rep(0,s0), rep(c.alp, J-s0))
  }
  eta<-pi0/2
  
  v.vec<-rnorm(n)
  D = Z %*% gam + v.vec
  
  if (mm==1){ #logit-i
    u.vec <-u1(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<-  sapply(plogis(D * beta0 +Z%*%kappa + u.vec), 
                function(x)  rbinom(1, 1, x))
    u1.r<-u1(z=z0, eta=eta, n=rn, rho.err=rho.err, v=rnorm(rn)) 
    cace0 <- mean(plogis(as.numeric(d1 * beta0 + z0 %*% kappa)+ u1.r )) - 
      mean( plogis(as.numeric(d2 * beta0 + z0 %*% kappa) + u1.r))
    asf0= mean(plogis(as.numeric(d2 * beta0 + z0 %*% kappa)+ u1.r ))
  }else if (mm==2){ #logit-ii
    u.vec <- u2(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<- sapply(plogis((D * beta0 + Z%*%kappa+u.vec)+(D * beta0 + Z%*%kappa+u.vec)^2/3), function(x) rbinom(1,1,x))
    u2.r<-u2(z=z0,eta=eta,n=rn, rho.err=rho.err, v=rnorm(rn))
    cace0 <- mean(plogis(as.numeric(d1 * beta0 + z0 %*% kappa)+u2.r +(as.numeric(d1 * beta0 + z0 %*% kappa)+u2.r)^2/3)) - 
      mean( plogis(as.numeric(d2 * beta0 + z0 %*% kappa)+u2.r+ (as.numeric(d2 * beta0 + z0 %*% kappa)+u2.r)^2/3))
    asf0<- mean(plogis(as.numeric(d2 * beta0 + z0 %*% kappa)+u2.r +(as.numeric(d2 * beta0 + z0 %*% kappa)+u2.r)^2/3)) 
  }
  
  
  list(
    Z = Z, D = D, Y = Y,
    V0 = 1:s0, 
    eta=eta,
    kappa=kappa,
    asf0=asf0,
    cace0 = cace0,
    z0 = z0
  )
}

beta0 = 0.25
rho.z = 0
rho.err= 0.25
s0=5
d1 = -1
d2 = 2
J=7
z0<-c(rep(0,J-1), 0.1)
normal.Z=T #normal.Z=T or F
Nround <- 500 
mm=2 ##mm=1 or 2

set.seed(123)
for(n in c(1000, 500)){
  for(c.gam in c(0.8, 0.6, 0.4)){
    cov=0
    cov.oracle<- 0
    re.semi.mat = matrix(NA, ncol=5, nrow = Nround)  ##### result matrix
    ########
    maj.fail.temp=0
    for(round in 1:Nround){
      #generate data which satisfy the inside assumption
      dat0 <- gen.InSIDE.bidata(n=n, c.gam = c.gam, beta0=beta0,
                                sig.v.sq = 1, sig.u.sq = 1, rho.err = 0.25, J=J, 
                                s0=5, d1=d1, d2=d2, z0=z0, c.alp=0.8, mm=mm, normal.Z=normal.Z)
      Y<-as.numeric(dat0$Y)
      D <- dat0$D
      Z <- dat0$Z
      sim.re <- main.semi(Y, D, Z,  bootstrap = T, bs.Niter = 40,
                          d1 = d1, d2 = d2, z0 = z0, M=2, M.est=T)
      maj.fail.temp=maj.fail.temp + sim.re$Maj.fail
      oracle.re<-main.semi(Y, D, Z,  bootstrap = T, bs.Niter = 40, V=1:5,
                           d1 = d1, d2 = d2, z0 = z0, M=2, M.est=T)
      re.semi.mat[round,] <- as.numeric(c(dat0$cace0, sim.re$cace.hat, sim.re$cace.sd, oracle.re$cace.hat, oracle.re$cace.sd))
      cov<- cov+ (abs(re.semi.mat[round,1]-re.semi.mat[round,2])<=1.96* re.semi.mat[round,3])#oracle coverage
      cat('round=', round, re.semi.mat[round,],cov/round,'\n')
    }#round
    sum.tsth.cace <- sum.result(re.semi.mat[,1:3])
    sum.oracle.cace<- sum.result(re.semi.mat[,c(1,4,5)])
    result.semi <- c(Nround, n,  c.gam, 
                     sum.tsth.cace$mae, sum.tsth.cace$ci.cov.re, sum.tsth.cace$ci.len.re, maj.fail.temp, 
                     sum.oracle.cace$mae,sum.oracle.cace$ci.cov.re, sum.oracle.cace$ci.len.re)#inference results
    
  }##c.gam
  
}






