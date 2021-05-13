library(Matrix)
gen.bidata <- function(n=1000, c.gam=0.8, beta0, sig.v.sq=1, sig.u.sq=1, rho.err=0.25, J=6, s0=5, 
                       d1=NULL, d2=NULL, z0=NULL, mm=1,c.alp=0.4,normal.Z=T, rn=2000){
  u1<-function(z,eta,n,rho.err, v){
    as.numeric( z%*% eta)+ v*rho.err + rnorm(n) 
  }
  u2<-function(z,eta,n,rho.err, v){
    exp(as.numeric(z%*% eta)+ rho.err*v) + runif(n,-1,1)
  }
  u3<-function(z,eta1,eta2,n,rho.err, v){
    exp(as.numeric(z%*% eta1)+ rho.err*v)+as.numeric(z%*%eta2+v) + runif(n,-1,1)
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
    pi0 <- c(rep(0, s0), c.alp, c.alp/2)
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
  } else if(mm == 3){#continuous iii
    u.vec <- u1(z=Z, eta=eta, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<-  D * beta0 +Z%*%kappa+ u.vec + (D * beta0 +Z%*%kappa+ u.vec)^2/3
    u1.r<-u1(z=z0, eta=eta, n=rn, rho.err=rho.err, v=rnorm(rn)) 
    cace0 <- mean(as.numeric(d1 * beta0 +z0%*%kappa)+ u1.r + (as.numeric(d1 * beta0 +z0%*%kappa)+ u1.r)^2/3)-
      mean(as.numeric(d2 * beta0 +z0%*%kappa)+ u1.r + (as.numeric(d2 * beta0 +z0%*%kappa)+ u1.r)^2/3)
    asf0=mean(as.numeric(d2 * beta0 +z0%*%kappa)+ u1.r + (as.numeric(d2 * beta0 +z0%*%kappa)+ u1.r)^2/3)
  }else if(mm==4){#continuous iv
    eta1 <- eta; eta2<- eta
    eta1[J]<-0; eta2[J-1] <- 0
    u.vec <- u3(z=Z, eta1=eta1, eta2=eta2, n=n, rho.err=rho.err, v=v.vec) 
    kappa <- pi0-eta
    Y<-  as.numeric(D * beta0 + Z %*% kappa) * u.vec
    u4.r<-u3(z=z0,eta1=eta1, eta2=eta2,n=rn, rho.err=rho.err, v=rnorm(rn))
    cace0 <- mean(as.numeric(d1 * beta0 + z0 %*% kappa) * u4.r) - 
      mean(as.numeric(d2 * beta0 + z0 %*% kappa) *u4.r)
    asf0= mean(as.numeric(d2 * beta0 + z0 %*% kappa) *u4.r)
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

#get mae, coverage, and ci length from iterations
sum.result <- function(re.mat, alpha=0.95){
  if(anyNA(re.mat[,2]) | anyNA(re.mat[,3])){
    err.loc <- which(is.na(re.mat[,3])| is.na(re.mat[,2]))
    re.mat<-re.mat[-err.loc,]
  }
  mae <- median(abs(re.mat[,2] -re.mat[,1]))
  ci.cov.re <- mean(re.mat[,1] < re.mat[,2] -qnorm(1/2-alpha/2) * re.mat[,3]& re.mat[,1] > re.mat[,2] +qnorm(1/2-alpha/2) * re.mat[,3])
  ci.len.re <- mean(re.mat[,3])
  return(list(mae = mae, ci.cov.re = ci.cov.re, ci.len.re = ci.len.re))
}

#main function
main.semi<- function(Y, D, Z, bs.Niter=32, M=2, bootstrap = T, M.est=T, V=NULL,
                     d1, d2 , z0){
  #M : number of directions
  Maj.fail=F # whether majority rule fails or not given by voting
  pz<- ncol(Z)
  n <- length(Y)
  #first-stage regression
  gam.hat<- lm(D ~ Z-1)$coef
  v.hat<-D-Z%*%gam.hat
  gram.Z <- t(Z)%*%Z/n
  sig.v.hat<- mean(v.hat^2)
  #voting and applying the majority rule
  if(is.null(V)){
    #get reduced-from
    SIR.re <-SIR.est(X.cov=cbind(Z,v.hat), Y, M= M, M.est=M.est)
    cat('M.hat=', ncol(SIR.re$theta.hat),'\n')
    M <- ncol(SIR.re$theta.hat)
    Gam.hat<-as.matrix(SIR.re$theta.hat[1:pz,], ncol=M) #estimate Gamma
    ##voting
    Select.re<-TSHT.VHat(n=n,ITT_Y=Gam.hat[,1], ITT_D=gam.hat, Cov.gGam=bdiag(sig.v.hat*solve(gram.Z),SIR.re$vGam[-(pz+1),-(pz+1)]))
    SHat<-Select.re$SHat
    if(length(Select.re$VHat)< length(SHat)/2){
      cat('Majority rule fails.','\n')
      Maj.fail=T
      #return(list(cace.hat=NA,cace.sd=NA, Maj.fail=Maj.fail))
    }
    beta.hat<-sapply(1:M, function(m) median(Gam.hat[SHat,m]/gam.hat[SHat]))
    pi.hat<- Gam.hat - gam.hat %*% matrix(beta.hat,nrow=1,ncol=M)
  }else{###oracle method
    SIR.re <-SIR.est(X.cov=cbind(Z%*%gam.hat,Z[,-V],v.hat), Y, M= M, M.est=M.est)
    M <- ncol(SIR.re$theta.hat)
    beta.hat<-SIR.re$theta.hat[1,]
    pi.hat<-matrix(0,nrow=pz, ncol=M)
    pi.hat[V,]<-0
    pi.hat[-V,]<-SIR.re$theta.hat[2:(pz-length(V)+1),]
  }
  ##estimate cace##
  bw.z<- apply(cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat), 2, function(x) 0.9*n^(-1/6)*min(sd(x), (quantile(x,0.75)-quantile(x,0.25))/1.34))
  #bw.z<-rep(silverman(M+1, n),M+1)
  asf.dw<-ASF.est(d1=d1,d2=d2, z0=z0, beta.hat= beta.hat, pi.hat= pi.hat, 
                  Y, D, Z, v.hat=v.hat, bw.z=bw.z)
  cace.hat = asf.dw$cace.hat #cace
  bw.z=asf.dw$bw.z
  ####bootstrap with paralell computing####
  boot_b <- foreach(i=1:bs.Niter, .combine=c) %dopar% {
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    
    list(boot.fun(data=bootstrap_data, M=M, d1=d1,d2=d2, z0=z0, SHat=SHat, 
                  bw.z=bw.z, V=V))
  }
 # cace.sd<-sd(unlist(lapply(boot_b, function(x) x[1])), na.rm=T)
  cace.sd<-sqrt(mean((unlist(lapply(boot_b, function(x) x[1]))-cace.hat)^2))
  return(list(gam.hat=gam.hat, pi.hat=pi.hat, beta.hat=beta.hat, Maj.fail=Maj.fail,
              cace.hat=cace.hat, cace.sd= cace.sd))
}

SIR.est<- function(X.cov,Y, M=2, M.est=T){
  p<- ncol(X.cov)
  n<-length(Y)
  if(M.est){
    nslice=ifelse(length(table(Y))==2,2,8)
    SIR.re<-dr(Y ~ X.cov -1, method='sir', numdir=M, nslices=nslice)
    evalues=SIR.re$evalues
    nobs.slice <- median(SIR.re$slice.info$slice.sizes)
    M <- which.max(
      sapply(1:M, function(m) sum(log(evalues[(m+1):p]+1)-evalues[(m+1):p])*n/2-log(n)*m*(2*p-m+1)/4/nobs.slice))
  }
  if(length(table(Y))==2){
    init.re<- glm(Y~X.cov-1,family=binomial(link='logit'))
  }else{
    init.re<-lm(Y~X.cov-1)
  }
  Gam.init<-init.re$coef
  vGam<-vcov(init.re)*n
  if(M==1){
    theta.hat <- orthoDr_reg(x=X.cov, y = Y, B.init=as.matrix(Gam.init,ncol=1), ndr=1)$B
  }else{
    theta.hat<-dr(Y ~ X.cov -1, method='sir', numdir=M)$evectors[,1:M]
    #theta.hat2 <- orthoDr_reg(x=X.cov, y = Y, B.init=theta.hat, ndr=2, maxitr=100)$B
  }

  list(theta.hat=theta.hat, vGam=vGam)
}

TSHT.VHat <- function(n, ITT_Y,ITT_D, Cov.gGam, tuning = 2.01, majority=T) {
  #ITT_Y=Gam.hat; ITT_D=gam.hat; Cov.gGam=rbind(cov.gGam1,cov.gGam2)
  Var.comp.est <- function(Cov.mat, gam.hat, j){
    diag(Cov.mat) + (gam.hat/gam.hat[j])^2 * Cov.mat[j,j] - 2*gam.hat/gam.hat[j] * Cov.mat[j,] 
  }
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)
  # Check Sigmas 
  stopifnot(!missing(Cov.gGam))
  
  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)
  
  # Constants
  pz = length(ITT_Y)
  
  # First Stage
  Var.gam.hat <- diag(Cov.gGam)[1:pz]
  SHat = (1:pz)[abs(ITT_D) >= (sqrt(Var.gam.hat) * sqrt(tuning*log(pz)/n))]
  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE
  
  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat
  for(j in SHat) {
    beta.j = ITT_Y[j] / ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    #compute three components in eq(33)
    Sig1.j <- Var.comp.est(Cov.gGam[1:pz,1:pz], ITT_D, j)
    Sig2.j <- Var.comp.est(Cov.gGam[(pz+1):(2*pz),(pz+1):(2*pz)], ITT_D, j)
    Sig3.j <-  Var.comp.est(Cov.gGam[1:pz,(pz+1):(2*pz)], ITT_D, j)
    sigmasq.j <- beta.j^2 *Sig1.j +  Sig2.j - 2* beta.j * Sig3.j
    PHat.bool.j <- abs(pi.j) <= sqrt(sigmasq.j) * tuning * sqrt(log(pz)/n)
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }
  
  # Voting
  diag(VHats.bool) <- rep(TRUE, nCand)
  VM = rowSums(VHats.bool)
  cat(VM,'\n')
  VHat = rownames(VHats.bool)[VM > (0.5 * length(SHat))] # Majority winners
  #VM.p = rownames(VHats.bool)[max(VM) == VM] #Plurality winners
  #VHat = as.numeric(union(VM.m,VM.p))
  #cat('m', VM.m, 'p',VM.p,'\n')
  # Error check
  # if(length(VHat) == 0){
  #   warning("VHat Warning: No valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.")
  #   warning("Defaulting to all IVs being valid")
  #   VHat = 1:pz
  # }
  
  return(list(VHat = VHat,SHat=SHat))
} 



boot.fun<-function(data, M, d1, d2,z0, SHat, gam.bs=NULL, beta.bs=NULL, pi.bs=NULL, bw.z=NULL, V=NULL){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  pz<-ncol(Z)
  if(is.null(beta.bs)){
    gam.bs<-lm(D~Z-1)$coef
    v.bs <- D- Z%*%gam.bs
    if(is.null(V)){
      SIR.bs.re <- SIR.est(cbind(Z,v.bs), Y, M=M, M.est=F)
      Gam.bs<-as.matrix(SIR.bs.re$theta.hat[1:pz,], ncol=M)
      beta.bs<-sapply(1:M, function(m) median(Gam.bs[SHat,m]/gam.bs[SHat]))
      pi.bs <- Gam.bs - gam.bs %*% matrix(beta.bs,nrow=1,ncol=M)
    }else{
      SIR.bs <-SIR.est(X.cov=cbind(Z%*%gam.bs,Z[,-V],v.bs), Y, M= M, M.est=F)
      M <- ncol(SIR.bs$theta.hat)
      beta.bs<-SIR.bs$theta.hat[1,]
      pi.bs<-matrix(0,nrow=pz, ncol=M)
      pi.bs[V,]<-0
      pi.bs[-V,]<-SIR.bs$theta.hat[2:(pz-length(V)+1),]
    }
    
  }
  asf.dw<-ASF.est(d1=d1,d2=d2, z0=z0, beta.hat= beta.bs, pi.hat= pi.bs, 
                  Y=Y, D=D, Z=Z, v.hat=D-Z%*%gam.bs, bw.z=bw.z)
  asf.dw$cace.hat
  
}

box.ker<-function(xx, X, h){
  apply(X, 1, function(xi) max(abs(xi-xx)/h)<1)/prod(h)
}

ASF.est <- function(d1, d2, z0, beta.hat, pi.hat,  Y, D, Z, v.hat, bw.z=NULL){
  beta.hat=as.vector(beta.hat)
  M=length(beta.hat)
  n=length(D)
  v.hat <- as.matrix(v.hat,ncol=1)
  D<-as.matrix(D,ncol=1)
  index <- cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat)
  
  if(is.null(bw.z)){
    bw.z0 <- matrix(0,ncol=M+1,nrow=10)
    bw.z0<- matrix(rep(seq(0.4,2,length.out=10),M+1), ncol=M+1) * (n*4/5)^(-1/(5+M))
    mse.bs <- foreach(i=1:nrow(bw.z0), .combine=c) %dopar% {#5-fold CV
      g.hat(index=index, Y=Y, bw=bw.z0[i,])
    }
    bw.z<- bw.z0[which.min(mse.bs),] * (n*4/5)^(1/(5+M))*n^(-1/5)
  }
  
  #cace
  index1.z <- cbind(apply(matrix(d1%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat) 
  index2.z <- cbind(apply(matrix(d2%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat) 
  #cat(index1.z[1,1],index2.z[1,1],'\n')
  q1<-quantile(index[,1], 0.975)
  q2<-quantile(index[,1], 0.025)
  
  if(sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0 & sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0){
    index1.z<- index1.z[index1.z[,1]<= q1 & index1.z[,1]>= q2,]
    index2.z<- index2.z[index2.z[,1]<= q1 & index2.z[,1]>= q2,]
  }
  asf.dw1<- NA
  asf.dw2<-NA
  bw.z<-bw.z/1.5
  while(is.na(asf.dw1) | is.na(asf.dw2)){
    bw.z<-bw.z*1.5
    if(M==1){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] ])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2]])), na.rm=T)
    }else if(M==2){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
    }else if(M==3){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3] & abs(xx[4]-index[,4])<=bw.z[4]])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3] & abs(xx[4]-index[,4])<=bw.z[4]])), na.rm=T)
    }
  }

  
  cace.hat <- asf.dw1-asf.dw2
  
  
  list(cace.hat = cace.hat, ace.hat=0, bw=bw.z, bw.z=bw.z)
}

# plugging in the scaled estimate of beta and pi, not the original
############

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

