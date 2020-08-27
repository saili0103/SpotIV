source("~/SpotIV-main.r")
##estimation and construct confidence intervals for CATE(d1,d2|z0)

set.seed(123)
mm=2 #1 & 2 are binary outcomes; 3& 4 are continuous outcomes
J=7 # number of IVs
d1=0
d2=2
z0=c(rep(0,J-2), -0.2,0.2)
dat0 <- gen.bidata(n=1000, c.gam = 1.2, beta0=0.25,
                   sig.v.sq = 1, sig.u.sq = 1, rho.err = 0.25, J=J, 
                   s0=5, d1=d1, d2=d2, z0=z0, c.alp=0.8, mm=mm, normal.Z=F)
Y<-as.numeric(dat0$Y)
D <- dat0$D
Z <- dat0$Z
sim.re <- main.semi(Y, D, Z,  bootstrap = T, bs.Niter = 50,
                    d1 = d1, d2 = d2, z0 = z0, M=3, M.est=T)

c(dat0$cace0, sim.re$cace.hat, sim.re$cace.sd) ###true CATE(d1,d2|z0), estimated CATE(d1,d2|z0), standard error of the estimate

ci<-c(sim.re$cace.hat-1.96*sim.re$cace.sd, sim.re$cace.hat+1.96*sim.re$cace.sd) ##95% two-sided CI
ci
