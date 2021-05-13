###control funciton method which assumes all the IVs are valid###
library(MASS)
library(mvtnorm)
library(dr)
library(np)
library(foreach)
library(orthoDr)
library(doParallel)
library(caret)
registerDoParallel(4)


source("SpotIV-main.R")

control.semi<- function(Y, D, Z, bs.Niter=40, bootstrap = T, d1, d2 , z0){
  #M : number of directions
  pz<- ncol(Z)
  n <- length(Y)
  gam.hat<- lm(D ~ Z-1)$coef
  v.hat<-D-Z%*%gam.hat
  bw.z<- apply(cbind(D,v.hat), 2, function(x) 0.9*n^(-1/6)*min(sd(x), (quantile(x,0.75)-quantile(x,0.25))/1.34))
  #bw.z<-rep(silverman(2, n),2)
  asf.dw<- ASF.est(d1=d1,d2=d2, z0=z0, Y, D, v.hat=v.hat,bw.z=bw.z)
  cace.hat = asf.dw$cace.hat #cate
  bw.z=asf.dw$bw.z
  cat(bw.z,'\n')
  ####bootstrap####
  boot_b <- foreach(i=1:bs.Niter, .combine=c) %dopar% {
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    list(boot.fun(data=bootstrap_data, d1=d1,d2=d2, z0=z0, bw.z=bw.z))
  }
  #cace.sd<-sd(unlist(lapply(boot_b, function(x) x[1])), na.rm=T)
  cace.sd<-sqrt(mean((unlist(lapply(boot_b, function(x) x[1]))-cace.hat)^2))
  return(list(cace.hat=cace.hat, cace.sd= cace.sd))
}


boot.fun<-function(data, M, d1, d2,z0, bw.z=NULL){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  pz<-ncol(Z)
  v.bs <- lm(D~Z-1)$res
  asf.dw<-ASF.est(d1=d1,d2=d2, z0=z0, Y=Y, D=D, v.hat=v.bs, bw.z=bw.z)
  asf.dw$cace.hat
  
}

box.ker<-function(xx, X, h){
  apply(X, 1, function(xi) max(abs(xi-xx)/h)<1)/prod(h)
}

ASF.est <- function(d1, d2, z0, Y, D, v.hat, bw.z=NULL){
  n=length(D)
  v.hat <- as.matrix(v.hat,ncol=1)
  D<-as.matrix(D,ncol=1)
  index <- cbind(D,v.hat)
  if(is.null(bw.z)){
    bw.z0 <- matrix(0,ncol=2,nrow=10)
    bw.z0<- cbind(seq(0.4,1.6,length.out=10), seq(0.4,1.6,length.out=10))* (n*4/5)^(-1/6)
    mse.bs <- foreach(i=1:nrow(bw.z0), .combine=c) %dopar% {#5-fold CV
      g.hat(index=index, Y=Y, bw=bw.z0[i,])
    }
    bw.z<- bw.z0[which.min(mse.bs),] * (n*4/5)^(1/6)*n^(-1/5)
  }
  #cace
  index1.z <- cbind(rep(d1,n),v.hat) 
  index2.z <- cbind(rep(d2,n),v.hat) 
  asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] ])), na.rm=T)
  asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2]])), na.rm=T)

  cace.hat <- asf.dw1-asf.dw2
  list(cace.hat = cace.hat, bw.z=bw.z)
}
g.hat<-function(Y, index, bw, cv=T){
  samp.eval<-createFolds(1:length(Y), k=5,list=T)
  mse<- 0
  if(!cv & ncol(index)==2){
    apply(index,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw[1] & 
                                         abs(xx[2]-index[,2])<=bw[2]]))
  }else if(!cv & ncol(index)==3){
    apply(index,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw[1] & 
                                         abs(xx[2]-index[,2])<=bw[2] & abs(xx[3]-index[,3])<=bw[3] ]))
  }else{
    for(fold in 1:5){
      samp.te<-samp.eval[[fold]]
      Y.tr<-Y[-samp.te]
      if(ncol(index)==3){
        ghat.est.cur<-apply(index[samp.te,],1, function(xx) mean(Y.tr[abs(xx[1]-index[-samp.te,1])<=bw[1] & abs(xx[2]-index[-samp.te,2])<=bw[2] &
                                                                        abs(xx[3]-index[-samp.te,3])<=bw[3]], na.rm=T))
        mse=mse+mean((Y[samp.te]-ghat.est.cur)^2, na.rm=T)
      }else if(ncol(index)==2){
        ghat.est.cur<-apply(index[samp.te,],1, function(xx) mean(Y.tr[abs(xx[1]-index[-samp.te,1])<=bw[1] & 
                                                                        abs(xx[2]-index[-samp.te,2])<=bw[2]], na.rm=T))
        mse=mse+mean((Y[samp.te]-ghat.est.cur)^2, na.rm=T)
      }
    }
    mse
  }
  
}
# plugging in the scaled estimate of beta and pi, not the original
############

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
mm=1
set.seed(123)
for(n in c(1000,500)){
  for(c.gam in c(0.8, 0.6, 0.4)){
    re.semi.mat = matrix(NA, ncol=3, nrow = Nround)  ##### result matrix
    cov=0
    ########
    for(round in 1:Nround){
      dat0 <- gen.bidata(n=n, c.gam = c.gam, beta0=beta0,
                                 sig.v.sq = 1, sig.u.sq = 1, rho.err = 0.25, J=J, 
                                 s0=5, d1=d1, d2=d2, z0=z0, c.alp=0.8, mm=mm, normal.Z=normal.Z)
      Y<-as.numeric(dat0$Y)
      D <- dat0$D
      Z <- dat0$Z
      sim.re <- control.semi(Y, D, Z,  bootstrap = T, bs.Niter = 40,
                          d1 = d1, d2 = d2, z0 = z0)
      re.semi.mat[round,] <- as.numeric(c(dat0$cace0, sim.re$cace.hat, sim.re$cace.sd))
      cov<- cov+ (abs(re.semi.mat[round,1]-re.semi.mat[round,2])<=1.96* re.semi.mat[round,3])
         cat('round=', round, re.semi.mat[round,], cov/round, '\n')
    }#round
    sum.tsth.cace <- sum.result(re.semi.mat[,1:3])
    result.semi <- c(Nround, n,  c.gam, 
                     sum.tsth.cace$mae, sum.tsth.cace$ci.cov.re, sum.tsth.cace$ci.len.re)#inference results

  }##c.gam
  
}






