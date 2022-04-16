#================================================================================================================
# This script is to run the preprocessed CTRP data by mix-lasso or tree lasso model
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 16-Apr-2022
#================================================================================================================

rm(list = ls())
library(mixlasso)

# load preproccessed CTRP data
load("cellminer_ctrp_auc_preselect.RData")
load("cellminer_ctrp_tissue_preselect.RData")
load("cellminer_ctrp_rna_preselect.RData")
load("cellminer_ctrp_cnv_preselect.RData")
load("cellminer_ctrp_mut_preselect.RData")

# create index for cancer types
t.idx <- numeric(dim(tissue.dum)[1])
for(i in 1:dim(tissue.dum)[2]) t.idx[tissue.dum[,i]==1] <- i

## create indicator matrices for missingness
auc.ctrp.mis <- matrix(0, nrow=dim(auc.ctrp)[1], ncol=dim(auc.ctrp)[2])
auc.ctrp.mis[is.na(auc.ctrp)] <- 1
rna.ctrp.mis <- matrix(0, nrow=dim(rna.ctrp)[1], ncol=dim(rna.ctrp)[2])
rna.ctrp.mis[is.na(rna.ctrp)] <- 1
cnv.ctrp.mis <- matrix(0, nrow=dim(cnv.ctrp)[1], ncol=dim(cnv.ctrp)[2])
cnv.ctrp.mis[is.na(cnv.ctrp)] <- 1
mut.ctrp.mis <- matrix(0, nrow=dim(mut.ctrp)[1], ncol=dim(mut.ctrp)[2])
mut.ctrp.mis[is.na(mut.ctrp)] <- 1
tissue.dum.mis <- matrix(0, nrow=dim(tissue.dum)[1], ncol=dim(tissue.dum)[2])
tissue.dum.mis[is.na(tissue.dum)] <- 1

## re-order all cell lines
auc.ctrp <- auc.ctrp[,order(t.idx)]
rna.ctrp <- rna.ctrp[,order(t.idx)]
cnv.ctrp <- cnv.ctrp[,order(t.idx)]
mut.ctrp <- mut.ctrp[,order(t.idx)]
tissue.dum <- tissue.dum[order(t.idx),]
auc.ctrp.mis <- auc.ctrp.mis[,order(t.idx)]
cnv.ctrp.mis <- cnv.ctrp.mis[,order(t.idx)]
mut.ctrp.mis <- mut.ctrp.mis[,order(t.idx)]
tissue.dum.mis <- tissue.dum.mis[order(t.idx),]
t.idx <- sort(t.idx)

## create foldid for cross-validation
set.seed(1)
foldid <- NULL
k <- 3+1
for(i in 1:ncol(tissue.dum))
  foldid <- c(foldid, sample(rep(seq(k),length=sum(t.idx==i))))

## use the kth fold data for validation and the other 4 folds for CV to choose lambda
t.idx_test <- t.idx[foldid==k]
auc.ctrp_test <- auc.ctrp[,foldid==k]
rna.ctrp_test <- rna.ctrp[,foldid==k]
cnv.ctrp_test <- cnv.ctrp[,foldid==k]
mut.ctrp_test <- mut.ctrp[,foldid==k]
tissue.dum_test <- tissue.dum[foldid==k,]
auc.ctrp.mis_test <- auc.ctrp.mis[,foldid==k]
rna.ctrp.mis_test <- rna.ctrp.mis[,foldid==k]
cnv.ctrp.mis_test <- cnv.ctrp.mis[,foldid==k]
mut.ctrp.mis_test <- mut.ctrp.mis[,foldid==k]
tissue.dum.mis_test <- tissue.dum.mis[foldid==k,]

t.idx <- t.idx[foldid!=k]
auc.ctrp <- auc.ctrp[,foldid!=k]
rna.ctrp <- rna.ctrp[,foldid!=k]
cnv.ctrp <- cnv.ctrp[,foldid!=k]
mut.ctrp <- mut.ctrp[,foldid!=k]
tissue.dum <- tissue.dum[foldid!=k,]
auc.ctrp.mis <- auc.ctrp.mis[,foldid!=k]
rna.ctrp.mis <- rna.ctrp.mis[,foldid!=k]
cnv.ctrp.mis <- cnv.ctrp.mis[,foldid!=k]
mut.ctrp.mis <- mut.ctrp.mis[,foldid!=k]
tissue.dum.mis <- tissue.dum.mis[foldid!=k,]
foldid <- foldid[foldid!=k]

# construct tree structre for drugs
tree.parm0 <- tree.parms(y=t(auc.ctrp))
drug.idx <- rownames(auc.ctrp) %in% tree.parm0$y.colnames

auc.tree <- auc.ctrp[drug.idx,]
auc.tree <- auc.tree[,!is.na(colSums(auc.tree))]
tree.parm0 <- tree.parms(y=t(auc.tree))
tree.parm0$pre0 <- pre.grad(tree.parm0$Tree, tree.parm0$Tw)

auc.ctrp <- auc.ctrp[drug.idx,]
auc.ctrp_test <- auc.ctrp_test[drug.idx,]
auc.ctrp.mis <- auc.ctrp.mis[drug.idx,]
auc.ctrp.mis_test <- auc.ctrp.mis_test[drug.idx,]

# prepare multiple data sources, i.e. GEX/CNV/MUT, or GEX/MUT
all.ctrp <- t(rbind(rna.ctrp , cnv.ctrp , mut.ctrp))
all.ctrp_test <- t(rbind(rna.ctrp_test , cnv.ctrp_test , mut.ctrp_test))
all.ctrp.mis <- t(rbind(rna.ctrp.mis , cnv.ctrp.mis , mut.ctrp.mis))
all.ctrp.mis_test <- t(rbind(rna.ctrp.mis_test , cnv.ctrp.mis_test , mut.ctrp.mis_test))
auc.ctrp <- t(auc.ctrp)
auc.ctrp.mis <- t(auc.ctrp.mis)
auc.ctrp_test <- t(auc.ctrp_test)
auc.ctrp.mis_test <- t(auc.ctrp.mis_test)
p <- c( nrow(rna.ctrp), nrow(cnv.ctrp), nrow(mut.ctrp) )

# set searching interval for penalty parameters
bounds <- t(data.frame(lambda=log(c(1,10)), gamma=c(50, 200), ipf1=c(0.5,1.1), ipf2=log(c(0.01,0.5))))
colnames(bounds) <- c("lower", "upper")

fit <- mixPenaltyReg(x=all.ctrp, y=auc.ctrp, z=tissue.dum, x_test=all.ctrp_test, y_test=auc.ctrp_test, z_test=tissue.dum_test, p=p, foldid=foldid, num.nonpen=0, 
                          method="IPF-tree-lasso-re", alpha=.9, tree.parm=tree.parm0,
                          y.mis=auc.ctrp.mis, y.mis_test=auc.ctrp.mis_test, x.mis=all.ctrp.mis, x.mis_test=all.ctrp.mis_test, 
                          strata.surv=NULL, search.path=FALSE, EI.eps=0.01, fminlower=-10,intercept=T, bounds=bounds, bound.scale=c("exp",NA,NA,"exp"),
                          threshold=1e-4, tol=1e-5, mu=0.2, NoVar=50, N=7, min.iter=2, maxiter=500,seed=1234,parallel=T, verbose=TRUE, t.idx=t.idx, t.idx_test=t.idx_test, t.glasso=T, predict.re=T)
save(fit, file="ctrp_mix_lasso1.rda")

# prepare multiple data sources, i.e. GEX/CNV/MUT, or GEX/MUT for tree-lasso by removing missing observations
t.idx <- t.idx[!is.na(rowSums(auc.ctrp))]
all.ctrp <- all.ctrp[!is.na(rowSums(auc.ctrp)),]
tissue.dum <- tissue.dum[!is.na(rowSums(auc.ctrp)),]
foldid <- foldid[!is.na(rowSums(auc.ctrp))]
auc.ctrp <- auc.ctrp[!is.na(rowSums(auc.ctrp)),]
all.ctrp <- cbind(tissue.dum, all.ctrp)#[,1:1000]
all.ctrp_test <- all.ctrp_test[!is.na(rowSums(auc.ctrp_test)),]
tissue.dum_test <- tissue.dum_test[!is.na(rowSums(auc.ctrp_test)),]
t.idx_test <- t.idx_test[!is.na(rowSums(auc.ctrp_test))]
auc.ctrp_test <- auc.ctrp_test[!is.na(rowSums(auc.ctrp_test)),]
all.ctrp_test <- cbind(tissue.dum_test, all.ctrp_test)#[,1:1000]

##select colorectal(n=17,t.idx=4)/glioma(n=16,t.idx=7)/melanoma(n=23,t.idx=12)/ovary(n=17,t.idx=14)
i <- 4
foldid <- foldid[t.idx==i]
auc.ctrp <- auc.ctrp[t.idx==i,]
auc.ctrp_test <- auc.ctrp_test[t.idx_test==i,]
all.ctrp <- all.ctrp[t.idx==i,-c(1:max(t.idx))]
all.ctrp_test <- all.ctrp_test[t.idx_test==i,-c(1:max(t.idx))]
t.idx <- 0

beta.hat <- mixlasso(x=all.ctrp, y=auc.ctrp, x_test=all.ctrp_test, y_test=auc.ctrp_test, p=p, foldid=foldid, num.nonpen=max(t.idx), method="tree-lasso", 
                     tree.parm=tree.parm0, lambda=seq(10,100,length=10), intercept=T, threshold=1e-4, tol=1e-6, mu=0.2, maxiter=5000,seed=1234, parallel=T)
save(fit, file="ctrp_tree_lasso1.rda")
