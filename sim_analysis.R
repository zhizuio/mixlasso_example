#================================================================================================================
# This script is to run simulated data by mix-lasso or tree lasso model
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 17-Jun-2022
#================================================================================================================

rm(list = ls())
library(mixlasso)
source("sim_fun.R")

# use idx to set different random seeds
idx <- 1
set.seed(idx)
n <- 300
m <- 24*5
p <- 1000
t0 <- 10

## simulate learning dataset
# simulation scenario 1
sim <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(.5, .5,rep(.5,t0/2-2),rep(.5,t0/2)),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)
# simulation scenario 2
sim <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(-.5, -.5,rep(.5,t0/2-2),rep(.5,t0/2)),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)
# simulation scenario 3
sim <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(rep(.4,t0/5),rep(.6,t0/5),rep(.8,t0/5),rep(1.0,t0/5),rep(1.2,t0/5)),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)
# simulation scenario 4
sim <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(-.7,-.5,-.3,.2,.4,.6,.8,1,1.2,1.4),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)

y <- sim$Y
x <- sim$X
z <- sim$Z
t.idx <- sim$t.idx

## simulate validation dataset
# simulation scenario 1
sim_test <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(.5, .5,rep(.5,t0/2-2),rep(.5,t0/2)),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)
# simulation scenario 2
sim_test <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(-.5, -.5,rep(.5,t0/2-2),rep(.5,t0/2)),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)
# simulation scenario 3
sim_test <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(rep(.4,t0/5),rep(.6,t0/5),rep(.8,t0/5),rep(1.0,t0/5),rep(1.2,t0/5)),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)
# simulation scenario 4
sim_test <- sim1(p=p, n=n, m=m, rho=0.4, t0=t0, B.mix=matrix(rep(c(-.7,-.5,-.3,.2,.4,.6,.8,1,1.2,1.4),length(p)),ncol=length(p)), stratum.predictor=T,t.prob=rep(1/t0,t0),sigma.u=1,tg.cov=NULL)

y_test <- sim_test$Y
x_test <- sim_test$X
z_test <- sim_test$Z
t.idx_test <- sim_test$t.idx

y.mis <- matrix(0, nrow=dim(y)[1], ncol=dim(y)[2])
y.mis_test <- matrix(0, nrow=dim(y_test)[1], ncol=dim(y_test)[2])

##=======================
## add 10% missing values for y
mis.p <- .05 * 2
for(i in 1:t0){
  y.mis[t.idx==i,] <- matrix(c(sample(1:0, floor(sum(t.idx==i)*m/2), replace=T, prob=c(mis.p,1-mis.p)),rep(0,sum(t.idx==i)*m-floor(sum(t.idx==i)*m/2))), ncol=m, nrow=sum(t.idx==i), byrow=T)
  y.mis_test[t.idx_test==i,] <- matrix(c(sample(1:0, floor(sum(t.idx_test==i)*m/2), replace=T, prob=c(mis.p,1-mis.p)),rep(0,sum(t.idx_test==i)*m-floor(sum(t.idx_test==i)*m/2))), ncol=m, nrow=sum(t.idx_test==i), byrow=T)
}
  
y[y.mis==1] <- NA
y_test[y.mis_test==1] <- NA

# construct tree structure of response variables
tree.parm <- tree.parms(y=y, h=0.7)
tree.parm$pre0 <- pre.grad(tree.parm$Tree, tree.parm$Tw)

set.seed(1837)
foldid <- NULL
for(i in 1:max(t.idx))
  foldid <- c(foldid, sample(rep(seq(3),length=sum(t.idx==i))))

# set searching intervals for penalty parameters. If using log scale, mixlasso function needs set the corresponding "bound.scale" as "exp" otherwise "NA".
bounds <- t(data.frame(lambda=log(c(.5,20)), gamma=c(.01, 1)))
colnames(bounds) <- c("lower", "upper")

# run mix-lasso model with argument "method='IPF-tree-lasso-re'" 
beta.hat <- mixlasso(x=x, y=y, z=z, x_test=x_test, y_test=y_test, z_test=z_test, p=p, foldid=foldid, num.nonpen=0, method="IPF-tree-lasso-re", 
                          alpha=.9, tree.parm=tree.parm, y.mis=y.mis, y.mis_test=y.mis_test, 
                          search.path=FALSE, EI.eps=0.01, fminlower=-10,intercept=T, bounds=bounds, bound.scale=c("exp",NA),
                          threshold=1e-4, tol=1e-5, mu=0.2, NoVar=50, N=20, min.iter=10, maxiter=1000,seed=1234,parallel=T, verbose=TRUE, 
                          t.idx=t.idx, t.idx_test=t.idx_test, t.glasso=T, predict.re=T)

save(beta.hat, file="sim_mix_lasso1.rda")


## run tree-lasso after removing missing observations
foldid <- foldid[!is.na(rowSums(y))]
x <- x[!is.na(rowSums(y)),]
z <- z[!is.na(rowSums(y)),]
y <- y[!is.na(rowSums(y)),]
x_test <- x_test[!is.na(rowSums(y_test)),]
z_test <- z_test[!is.na(rowSums(y_test)),]
y_test <- y_test[!is.na(rowSums(y_test)),]

fit <- mixlasso(x=x, y=y, x_test=x_test, y_test=y_test, p=p, foldid=foldid, num.nonpen=0, method="tree-lasso", 
                     tree.parm=tree.parm, lambda=seq(10,100,length=10), intercept=T,
                     threshold=1e-4, tol=1e-5, mu=0.2, maxiter=5000,seed=1234, parallel=T)

save(fit, file="sim_tree_lasso1.rda")

