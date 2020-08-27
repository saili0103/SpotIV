library(MASS)
library(mvtnorm)
library(dr) ###for sliced-inverse regression
library(np)
library(foreach) ##for parallel computing of bootstrapped confidence intervals
library(doParallel) ##for parallel computing of bootstrapped confidence intervals
registerDoParallel(4) ##for parallel computing of bootstrapped confidence intervals
library(caret) ##for bandwith selection based on CV

####generating data for simulation
gen.bidata <- function(n=1000, c.gam=0.8, beta0, sig.v.sq=1, sig.u.sq=1, rho.err=0.25, J=6, s0=5, 
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
    pi0 <- c(rep(0, s0), - c.alp, -0.5, 0.5, c.alp)
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
  }else if (mm==2){ #logit-ii
    u.vec <- u2(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<- sapply(plogis((D * beta0 + Z%*%kappa+u.vec)+(D * beta0 + Z%*%kappa+u.vec)^2/3), function(x) rbinom(1,1,x))
    u2.r<-u2(z=z0,eta=eta,n=rn, rho.err=rho.err, v=rnorm(rn))
    cace0 <- mean(plogis(as.numeric(d1 * beta0 + z0 %*% kappa)+u2.r +(as.numeric(d1 * beta0 + z0 %*% kappa)+u2.r)^2/3)) - 
      mean( plogis(as.numeric(d2 * beta0 + z0 %*% kappa)+u2.r+ (as.numeric(d2 * beta0 + z0 %*% kappa)+u2.r)^2/3))
  } else if(mm == 3){#continuous g
    u.vec <- u1(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<-  D * beta0 +Z%*%kappa+ u.vec + (D * beta0 +Z%*%kappa+ u.vec)^2/3
    u1.r<-u1(z=z0, eta=eta, n=rn, rho.err=rho.err, v=rnorm(rn)) 
    cace0 <- mean(as.numeric(d1 * beta0 +z0%*%kappa)+ u1.r + (as.numeric(d1 * beta0 +z0%*%kappa)+ u1.r)^2/3)-
      mean(as.numeric(d2 * beta0 +z0%*%kappa)+ u1.r + (as.numeric(d2 * beta0 +z0%*%kappa)+ u1.r)^2/3)
  }else if(mm==4){
    u.vec <- u2(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<-  as.numeric(D * beta0 + Z %*% kappa)^3 * u.vec
    u2.r<-u2(z=z0,eta=eta,n=rn, rho.err=rho.err, v=rnorm(rn))
    cace0 <- mean(as.numeric(d1 * beta0 + z0 %*% kappa)^3 * u2.r) - 
      mean(as.numeric(d2 * beta0 + z0 %*% kappa)^3 *u2.r)
  }
  
  
  list(
    Z = Z, D = D, Y = Y,
    V0 = 1:s0, 
    eta=eta,
    kappa=kappa,
    cace0 = cace0,
    z0 = z0
  )
}


###main body of the function
###Y: outcome; D: exposure; Z: IVs; bs.Niter: number of bootstrap iteration;
###M: upper bound on the dimension; d1& d2& z0: output CATE(d1,d2|z=z0);
###tuning: tuning parameter for estimating Shat.
main.semi<- function(Y, D, Z, bs.Niter=50, M=3, bootstrap = T, M.est=T,
                     d1, d2 , z0, tuning=2.01){
  pz<- ncol(Z)
  n <- length(Y)
  gam.hat<- lm(D ~ Z-1)$coef
  v.hat<-D-Z%*%gam.hat
  SIR.re <-SIR.est(X.cov=cbind(Z,v.hat), Y, M= M, M.est=M.est)
  cat('M.hat', ncol(SIR.re$theta.hat),'\n')
  M <- ncol(SIR.re$theta.hat)
  Gam.hat<-as.matrix(SIR.re$theta.hat, ncol=M)
  Gam.hat<-as.matrix(SIR.re$theta.hat[1:pz,], ncol=M)
  cat(SIR.re$evalues[1:M],'\n')
  
  ##final estimator##
  gram.Z <- t(Z)%*%Z/n
  sig.v.hat<- mean(v.hat^2)
  SHat = which(abs(gam.hat) >= (sqrt(sig.v.hat*diag(solve(gram.Z))) * sqrt(tuning*log(pz)/n)))
  beta.hat<-sapply(1:M, function(m) median(Gam.hat[SHat,m]/gam.hat[SHat]))
  pi.hat<- Gam.hat - gam.hat %*% matrix(beta.hat,nrow=1,ncol=M)
  
  asf.dw<- ASF.est(d1=d1,d2=d2, z0=z0, beta.hat= beta.hat, pi.hat= pi.hat, 
                   Y, D, Z, v.hat=v.hat)
  cace.hat = asf.dw$cace.hat #cate
  bw.z=asf.dw$bw.z
  cat(bw.z,'\n')
  ####bootstrap####
  boot_b <- foreach(i=1:bs.Niter, .combine=c) %dopar% {
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    list(boot.fun(data=bootstrap_data, M=M, d1=d1,d2=d2, z0=z0, SHat=SHat, 
                  bw.z=bw.z))
  }
  cace.sd<-sd(unlist(lapply(boot_b, function(x) x[1])), na.rm=T)
  return(list(gam.hat=gam.hat, pi.hat=pi.hat, beta.hat=beta.hat,
              cace.hat=cace.hat, cace.sd= cace.sd))
}

### SIR function with M estimation
SIR.est<- function(X.cov,Y, nterms, M=3, M.est=T){
  p<- ncol(X.cov)
  n<-length(Y)
  SIR.re<-NULL
  nslice=ifelse(length(table(Y))==2,2,8)
  SIR.re<-dr(Y ~ X.cov -1, method='sir', numdir=M, nslices=nslice)
  if(M.est){
    evalues=SIR.re$evalues
    #cat(evalues[1:M],'\n')
    nobs.slice <- median(SIR.re$slice.info$slice.sizes)
    M <- which.max(
      sapply(1:M, function(m) sum(log(evalues[(m+1):p]+1)-evalues[(m+1):p])*n/2-log(n)*m*(2*p-m+1)/4/nobs.slice))
    SIR.re<-dr(Y ~ X.cov -1, method='sir', numdir=M)
  }
  list(theta.hat=as.matrix(SIR.re$evectors[,1:M]), evalues=SIR.re$evalues)
}
boot.fun<-function(data, M, d1, d2,z0, SHat, gam.bs=NULL, beta.bs=NULL, pi.bs=NULL, bw.z=NULL){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  pz<-ncol(Z)
  if(is.null(beta.bs)){
    gam.bs<-lm(D~Z-1)$coef
    v.bs <- D- Z%*%gam.bs
    SIR.bs.re <- SIR.est(cbind(Z,v.bs), Y, M=M, M.est=F)
    # Gam.bs<-as.matrix(SIR.bs.re$theta.hat, ncol=M)
    Gam.bs<-as.matrix(SIR.bs.re$theta.hat[1:pz,], ncol=M)
    beta.bs<-sapply(1:M, function(m) median(Gam.bs[SHat,m]/gam.bs[SHat]))
    pi.bs <- Gam.bs - gam.bs %*% matrix(beta.bs,nrow=1,ncol=M)
  }
  asf.dw<-ASF.est(d1=d1,d2=d2, z0=z0, beta.hat= beta.bs, pi.hat= pi.bs, 
                  Y=Y, D=D, Z=Z, v.hat=D-Z%*%gam.bs, bw.z=bw.z)
  asf.dw$cace.hat
  
}


####apply box kernel to estimate the function value at xx
box.ker<-function(xx, X, h){
  apply(X, 1, function(xi) max(abs(xi-xx)/h)<1)/prod(h)
}

#estimating the ASF(average structural function)
#beta.hat: coef of D; pi.hat: coef of Z; bw.z: bandwidth for kernel estimation
#v.hat: control variable
ASF.est <- function(d1, d2, z0, beta.hat, pi.hat,  Y, D, Z, v.hat, bw.z=NULL){
  beta.hat=as.vector(beta.hat)
  M=length(beta.hat)
  n=length(D)
  v.hat <- as.matrix(v.hat,ncol=1)
  D<-as.matrix(D,ncol=1)
  index <- cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat)
  
  if(is.null(bw.z)& M==1){
    bw.z0 <- matrix(0,ncol=2,nrow=10)
    bw.z0<- cbind(seq(0.4,2,length.out=10), seq(0.4,2,length.out=10))* (n*4/5)^(-1/6)
    mse.bs <- foreach(i=1:nrow(bw.z0), .combine=c) %dopar% {#5-fold CV
      g.hat(index=index, Y=Y, bw=bw.z0[i,])
    }
    bw.z<- bw.z0[which.min(mse.bs),] * (n*4/5)^(1/6)*n^(-1/5)
    
  }else if(is.null(bw.z)& M==2){
    bw.z0 <- matrix(0,ncol=3,nrow=10)
    bw.z0<- cbind(seq(0.4,2,length.out=10), seq(0.4,2,length.out=10), seq(0.4,2,length.out=10))* (n*4/5)^(-1/7)
    mse.bs <- foreach(i=1:ncol(bw.z0), .combine=c) %dopar% {
      g.hat(index=index, Y=Y, bw=bw.z0[i,])
    }
    bw.z<- bw.z0[which.min(mse.bs),]* (n*4/5)^(1/7)*n^(-1/6)
  }
  
  #cace
  index1.z <- cbind(apply(matrix(d1%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat) 
  index2.z <- cbind(apply(matrix(d2%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat) 
  q1<-quantile(index[,1], 0.975)
  q2<-quantile(index[,1], 0.025)
  
  if(sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0 & sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0){
    index1.z<- index1.z[index1.z[,1]<= q1 & index1.z[,1]>= q2,]
    index2.z<- index2.z[index2.z[,1]<= q1 & index2.z[,1]>= q2,]
  }
  if(M==1){
    asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] ])), na.rm=T)
    asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2]])), na.rm=T)
  }else if (M==2){
    asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
    asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
  }
  
  cace.hat <- asf.dw1-asf.dw2
  
  
  list(cace.hat = cace.hat, ace.hat=0, bw=bw.z, bw.z=bw.z)
}
############


#CV estimation of g function for bandwith selection
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

