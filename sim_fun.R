#================================================================================================================
# This script is to simulate heterogeneous samples, and one or multiple heterogeneous data sources (covariates) with multivariate response variables.
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 17-Jun-2022
#================================================================================================================

library(Matrix)
library(mvnfast)

sim <- function(p=500,n=100,m=24,rho=.4,B.elem=0.2,B.mix=NULL,t0=0,stratum.predictor=FALSE,t.prob=c(.1,.2,.3,.4),sigma.u=4,shift=0,tg.cov=NULL){
    b<-10
    # generate covariance matrix for X
    Ap1<-matrix(rep(rho,(p[1]/b)^2),nrow=p[1]/b)
    diag(Ap1)<-rep(1,p[1]/b)
    Xsigma1<-Ap1
    for(i in 2:b){
      Xsigma1<-bdiag(Xsigma1,Ap1)
    }
    
    Xsigma<-Xsigma1
    if(is.null(tg.cov)){
      X<-	mvnfast::rmvn(n,mu=rep(0,p),sigma=Xsigma)
    }
    
    # generate uncorrelated error term
    esd<-diag(m)
    e<-mvnfast::rmvn(n,mu=rep(0,m),sigma=esd)
    
    t.idx<-1
    Beta.mix<-NULL
    B.elem <- B.mix[1,]
    repeat{
      ## generate beta1 matrix
      Beta1<-matrix(0,nrow=m,ncol=p[1])
      Beta1[,shift*(t.idx-1)+1]<-B.elem[1]
      for(i in 1:2){
        Beta1[((i-1)*m/2+1):(i*m/2),shift*(t.idx-1)+(1+(i-1)*2+1):(1+i*2)]<-B.elem[1]
      }
      for(i in 1:4){
        Beta1[((i-1)*m/4+1):(i*m/4),shift*(t.idx-1)+(1+2*2+(i-1)*4+1):(1+2*2+i*4)]<-B.elem[1]
      }
      for(i in 1:8){
        Beta1[((i-1)*m/8+1):(i*m/8),shift*(t.idx-1)+(1+2*2+4*4+(i-1)*8+1):(1+2*2+4*4+i*8)]<-B.elem[1]
      }
      
      Beta<-t(Beta1)
      
      if(stratum.predictor & t.idx<=t0){
        Beta.mix <- rbind(Beta.mix, Beta)
        
        t.idx <- t.idx + 1
        if(t.idx > t0){
          Beta <- Beta.mix
          break
        }
        B.elem <- B.mix[t.idx,]
      }else{
        break
      }
    }
  
  if(t0 == 0){
    b.re <- 0
    Z <- 0
    Y <- X%*%Beta+e
  }else{
    # add random effects
    Z <- t( rmultinom(n, 1 ,t.prob) )
    t.idx <- sort(apply(Z, 1, function(xx) which(xx==1)))
    b.re <- NULL
    
    if(!is.null(tg.cov)){
      Xsigma <- matrix(rho, ncol=n, nrow=n)
      for(i in 1:t0){
        tmp <- matrix(tg.cov, nrow=sum(t.idx==i), ncol=sum(t.idx==i))
        diag(tmp) <- 8
        Xsigma[sum(t.idx<i)+1:sum(t.idx==i),sum(t.idx<i)+1:sum(t.idx==i)] <- tmp
      } 
      
      X<-	t(mvnfast::rmvn(p,mu=rep(0,n),sigma=Xsigma))
    }
    
    if(stratum.predictor){
      Y <- matrix(ncol=m,nrow=n)
      for(i in 1:t0){
        UEsigma <- matrix(sigma.u^2, nrow=sum(t.idx==i), ncol=sum(t.idx==i))
        diag(UEsigma) <- sigma.u^2+1
        for(j in 1:m)
          Y[sum(t.idx<i)+1:sum(t.idx==i),j] <- mvnfast::rmvn(1,mu=X[sum(t.idx<i)+1:sum(t.idx==i),]%*%Beta[(i-1)*sum(p)+1:sum(p),j],sigma=UEsigma)
      }
      b.re <- NULL
    }else{
      for(i in 1:t0){
        UEsigma <- matrix(sigma.u^2, nrow=sum(t.idx==i), ncol=sum(t.idx==i))
        diag(UEsigma) <- diag(UEsigma) * 1.1
        b.re <-	rbind(b.re, t( mvnfast::rmvn(m,mu=rep(0,sum(t.idx==i)),sigma=UEsigma) ))
      }
      Y <- X%*%Beta + b.re
    }
    
  }
  
  return(list(Y=Y, X=X, Beta=Beta, e=e, p=p, Z=Z, b.re=b.re, t.idx=t.idx))
}
