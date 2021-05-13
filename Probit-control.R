###control funciton method in probit outcome models with possibly invalid IVs###
## We deal with invalid IVs with majority rules
library(MASS)
library(mvtnorm)
library(foreach)   #on windows machine, you may need different packages for parallel computing,
library(doParallel)  #such as "library(parallel);detectCores()"
library(caret)
registerDoParallel(4)
source("SpotIV-main.R")

control.probit<- function(Y, D, Z, d1, d2 , z0, bootstrap=T, bs.Niter=40){
  pz<- ncol(Z)
  n <- length(Y)
  gam.hat<- lm(D ~ Z-1)$coef
  v.hat<-D-Z%*%gam.hat
  gram.Z <- t(Z)%*%Z/n
  sig.v.hat<- mean(v.hat^2)
  gam.cov<-sig.v.hat*solve(gram.Z)
  Gam.re<-glm(Y~cbind(Z,v.hat)-1, family=binomial(link='probit'))
  Gam.hat<-Gam.re$coef[1:pz]
  lam.hat<-Gam.re$coef[pz+1]
  Gam.cov<-n*vcov(Gam.re)[1:pz,1:pz]+lam.hat^2*gam.cov
  Cov.gGam<-rbind(cbind(gam.cov,lam.hat^2*gam.cov),
                  cbind(lam.hat^2*Gam.cov, Gam.cov))
  Select.re<-TSHT.VHat(n=n,ITT_Y=Gam.hat, ITT_D=gam.hat, Cov.gGam=Cov.gGam)
  SHat<-Select.re$SHat
  beta.hat<-median(Gam.hat[SHat]/gam.hat[SHat])
  pi.hat<- Gam.hat - gam.hat * beta.hat
  cace.hat = mean(pnorm(as.numeric(d1*beta.hat+z0%*%pi.hat) + v.hat*(lam.hat-beta.hat)))- 
    mean(pnorm(as.numeric(d2*beta.hat+z0%*%pi.hat) + v.hat*(lam.hat-beta.hat)))

  ####bootstrap####
  boot_b <- foreach(i=1:bs.Niter, .combine=c) %dopar% {
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    list(boot.fun(data=bootstrap_data, d1=d1,d2=d2, z0=z0, SHat=SHat))
  }
  #cace.sd<-sd(unlist(lapply(boot_b, function(x) x[1])), na.rm=T)
  cace.sd<-sqrt(mean((unlist(lapply(boot_b, function(x) x[1]))-cace.hat)^2))
  return(list(cace.hat=cace.hat, cace.sd= cace.sd))
}


boot.fun<-function(data, d1, d2,z0, SHat){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  pz<-ncol(Z)
  gam.bs<-lm(D~Z-1)$coef
  v.bs <- D- Z%*%gam.bs
  Gam.bs.re<-glm(Y~cbind(Z,v.bs)-1, family=binomial(link='probit'))
  Gam.bs<-Gam.bs.re$coef[1:pz]
  lam.bs<-Gam.bs.re$coef[pz+1]
  beta.bs<-median(Gam.bs[SHat]/gam.bs[SHat])
  pi.bs <- Gam.bs - gam.bs *beta.bs
  cace.bs = mean(pnorm(as.numeric(d1*beta.bs+z0%*%pi.bs) + v.bs*(lam.bs-beta.bs)))- 
    mean(pnorm(as.numeric(d2*beta.bs+z0%*%pi.bs) + v.hat*(lam.bs-beta.bs)))
  
  cace.bs
}
gen.probit.data <- function(n=1000, c.gam=0.8, beta0, sig.v.sq=1, sig.u.sq=1, rho.err=0.25, J=6, s0=5, 
                       d1=NULL, d2=NULL, z0=NULL,c.alp=0.4,normal.Z=T, rn=2000){
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
    pi0 <- c(rep(0, s0), c.alp, c.alp/2)
  }else{
    pi0 <- c(rep(0,s0), rep(c.alp, J-s0))
  }

  cov.noise<-matrix(c(sig.v.sq,sqrt(sig.v.sq*sig.u.sq)*rho.err, sqrt(sig.v.sq*sig.u.sq)*rho.err, sig.u.sq),ncol=2)
  noise.vec<-rmvnorm(n, rep(0,2), cov.noise)
  v.vec<-noise.vec[,1]
  D = Z %*% gam + v.vec
  
 u.vec<- noise.vec[,2]
 eta<-pi0/2
 kappa=pi0-eta
 Y<- (D * beta0 +Z%*%pi0 + u.vec>=0)
 
 u1.r<-rnorm(rn,0,sd=sqrt(sig.u.sq))
 cace0 <- mean((as.numeric(d1 * beta0 + z0 %*% pi0)+ u1.r )>=0) - 
                 mean((as.numeric(d2 * beta0 + z0 %*% pi0) + u1.r)>=0)
  list(
    Z = Z, D = D, Y = Y,
    V0 = 1:s0, 
    eta=eta,
    kappa=kappa,
    cace0 = cace0,
    z0 = z0
  )
}


####examples
beta0 = 0.25
rho.z = 0
rho.err= 0.25

d1 = -1
d2 = 2
J=7
z0<-c(rep(0,J-1), 0.1)
normal.Z=F
Nround <- 500 
set.seed(123)
for(n in c(1000,500)){
  for(c.gam in c(0.8, 0.6, 0.4)){
    re.probit.mat = matrix(NA, ncol=3, nrow = Nround)  ##### result matrix
    cov=0
    ########
    for(round in 1:Nround){
      dat0 <- gen.probit.data(n=n, c.gam = c.gam, beta0=beta0,
                                 sig.v.sq = 1, sig.u.sq = 1, rho.err = 0.25, J=J, 
                                 s0=5, d1=d1, d2=d2, z0=z0, c.alp=0.8, normal.Z=normal.Z)
      Y<-as.numeric(dat0$Y)
      D <- dat0$D
      Z <- dat0$Z
      probit.re <- control.probit(Y, D, Z,  bootstrap = T, bs.Niter = 40,
                          d1 = d1, d2 = d2, z0 = z0)
      re.probit.mat[round,] <- as.numeric(c(dat0$cace0, probit.re$cace.hat, probit.re$cace.sd))
      cov<- cov+ (abs(re.probit.mat[round,1]-re.probit.mat[round,2])<=1.96* re.probit.mat[round,3])
         cat('round=', round, re.probit.mat[round,], cov/round, '\n')
    }#round
    sum.cace <- sum.result(re.probit.mat[,1:3])
    result.semi <- c(Nround, n,  c.gam, 
                     sum.cace$mae, sum.cace$ci.cov.re, sum.cace$ci.len.re)#inference results

  }##c.gam
  
}






