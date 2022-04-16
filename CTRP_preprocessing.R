#================================================================================================================
# This script is to preproccess CTRP data
#
# author: Zhi Zhao (zhi.zhao@medisin.uio.no)
# date: 16-Apr-2022
#================================================================================================================

rm(list = ls())
library(PharmacoGx)
library(data.table)

## load CTRP pharmacological data downloaded from https://ocg.cancer.gov/programs/ctd2/data-portal
auc.ctrp0 <- read.table("v20.data.curves_post_qc.txt", header=T)
auc.ctrp00 <- auc.ctrp0[,c(1,17,16)]
compound <- read.table("v20.meta.per_compound.txt", fill=T)
compound <- compound[,1:2]
colnames(compound) <- compound[1,]
compound <- compound[-1,]
auc.ctrp1 <- merge(auc.ctrp00, compound, by="master_cpd_id")

experiment <- read.table("v20.meta.per_experiment.txt", header=T)
experiment <- experiment[,c(1,9)]
auc.ctrp2 <- merge(auc.ctrp1, experiment, by="experiment_id")

cell_line <- read.table("v20.meta.per_cell_line.txt", fill=T)
cell_line <- cell_line[,c(1,2)]
colnames(cell_line) <- cell_line[1,]
cell_line <- cell_line[-1,]
auc.ctrp3 <- merge(auc.ctrp2, cell_line, by="master_ccl_id")
auc.ctrp4 <- auc.ctrp3[,c(6,5,4)]
auc.ctrp5 <- auc.ctrp4[!duplicated(auc.ctrp4[,1:2]),]
y <- reshape(auc.ctrp5, v.names="area_under_curve", timevar="ccl_name", idvar="cpd_name", direction="wide")
y1 <- y
y1 <- y1[!y1[,1] %like% ":",]
rownames(y1) <- y1[,1]; y1 <- y1[,-1]
colnames(y1) <- substr(colnames(y1), 18, nchar(colnames(y1)))

auc.ctrp <- y1

### genomic data are from CCLE with the same measured cell lines in CTRP
availablePSets()
CCLE <- downloadPSet("CCLE_2015")

rna.ctrp <- assay(summarizeMolecularProfiles(CCLE, mDataType="rna"))
cnv.ctrp <- assay(summarizeMolecularProfiles(CCLE, mDataType="cnv"))
mutation.ctrp <- assay(summarizeMolecularProfiles(CCLE, mDataType="mutation", summary.stat="and"))
mut.ctrp <- matrix(as.numeric(mutation.ctrp), nrow=nrow(mutation.ctrp), ncol=ncol(mutation.ctrp))
rownames(mut.ctrp) <- rownames(mutation.ctrp)
colnames(mut.ctrp) <- colnames(mutation.ctrp)

## load CCLE tumor categorization
tissue.ccle <- read.csv("CCLE_tumor_types.csv", header=T, sep=";")
tissue.ccle$Cell.line.primary.name <- gsub("-","", tissue.ccle$X...Cell.line.primary.name)
tissue.ccle$Cell.line.primary.name <- gsub(" ","", tissue.ccle$Cell.line.primary.name)

tissue.ctrp <- data.frame(Name=colnames(auc.ctrp),cell.type=rep(NA, ncol(auc.ctrp)))
tissue.ctrp$Name <- gsub("-",".", tissue.ctrp$Name)
tissue.ctrp$Name <- gsub(" ",".", tissue.ctrp$Name)
tissue.ctrp <- tissue.ctrp[tissue.ctrp$Name %in% tissue.ccle$Cell.line.primary.name, ]
cell.type <- rep(NA, nrow(tissue.ctrp))
cell.type2 <- rep(NA, nrow(tissue.ctrp))
for(i in 1:length(cell.type)){
  cell.type[i] <- as.vector(tissue.ccle$CCLE.tumor.type[tissue.ccle$Cell.line.primary.name==tissue.ctrp$Name[i]])
  cell.type2[i] <- as.vector(tissue.ccle$Hist.Subtype1[tissue.ccle$Cell.line.primary.name==tissue.ctrp$Name[i]])
}
cell.type[cell.type=="lung_NSC"] <- cell.type2[cell.type=="lung_NSC"]

tissue.ctrp <- data.frame(name=tissue.ctrp$Name, cell.type=cell.type)
tissue.ctrp <- tissue.ctrp[tissue.ctrp$cell.type!="",]
# choose each tissue at least 20 cell lines?
tissue.ctrp <- tissue.ctrp[ tissue.ctrp$cell.type %in% names(table(tissue.ctrp$cell.type))[table(tissue.ctrp$cell.type)>15], ]
tissue.ctrp$cell.type <- as.vector(tissue.ctrp$cell.type)

### remove completely missing cell lines in CTRP
rna.ctrp <- rna.ctrp[,colSums(rna.ctrp,na.rm=T)!=0]
cnv.ctrp <- cnv.ctrp[,colSums(cnv.ctrp,na.rm=T)!=0]
mut.ctrp <- mut.ctrp[,colSums(mut.ctrp,na.rm=T)!=0]

## select cell lines
colnames(rna.ctrp) <- gsub("-","", colnames(rna.ctrp)); colnames(rna.ctrp) <- gsub(" ","", colnames(rna.ctrp))
colnames(cnv.ctrp) <- gsub("-","", colnames(cnv.ctrp)); colnames(cnv.ctrp) <- gsub(" ","", colnames(cnv.ctrp))
colnames(mut.ctrp) <- gsub("-","", colnames(mut.ctrp)); colnames(mut.ctrp) <- gsub(" ","", colnames(mut.ctrp))
name.cell <- intersect( intersect( intersect( intersect( colnames(auc.ctrp), colnames(rna.ctrp)), colnames(cnv.ctrp)), colnames(mut.ctrp)), tissue.ctrp$name)
auc.ctrp <- auc.ctrp[,colnames(auc.ctrp) %in% name.cell]
rna.ctrp <- rna.ctrp[,colnames(rna.ctrp) %in% name.cell]
cnv.ctrp <- cnv.ctrp[,colnames(cnv.ctrp) %in% name.cell]
mut.ctrp <- mut.ctrp[,colnames(mut.ctrp) %in% name.cell]
tissue.ctrp <- tissue.ctrp[tissue.ctrp$name %in% name.cell,]

## sort cell lines
tissue.ctrp <- tissue.ctrp[order(tissue.ctrp$cell.type),]
auc.ctrp <- auc.ctrp[,order(match(colnames(auc.ctrp), tissue.ctrp$name))]
rna.ctrp <- rna.ctrp[,order(match(colnames(rna.ctrp), tissue.ctrp$name))]
cnv.ctrp <- cnv.ctrp[,order(match(colnames(cnv.ctrp), tissue.ctrp$name))]
mut.ctrp <- mut.ctrp[,order(match(colnames(mut.ctrp), tissue.ctrp$name))]

### remove completely missing genomic features in CTRP
rna.ctrp <- rna.ctrp[!is.na(rowSums(rna.ctrp)),]
cnv.ctrp <- cnv.ctrp[!is.na(rowSums(cnv.ctrp)),]
mut.ctrp <- mut.ctrp[!is.na(rowSums(mut.ctrp)),]

## create tissue design matrix (dummy variables)
t.idx <- factor(tissue.ctrp$cell.type)
tissue.dum <- model.matrix(~ t.idx + 0)
colnames(tissue.dum) <- levels(t.idx)
levels(t.idx) <- 1:length(levels(t.idx))
t.idx <- as.numeric(t.idx)

# pre-select GEX and CNV by univariate analysis
pvalues1 <- rep(NA, nrow(rna.ctrp))
for(j in 1:nrow(rna.ctrp)){
  pvalues1[j] <- min(apply(auc.ctrp, 1, function(xx){cor.test(xx,rna.ctrp[j,])$p.value})*ncol(auc.ctrp), na.rm=T)
}
rna.ctrp1 <- rna.ctrp[pvalues1*nrow(rna.ctrp) < 0.1,]
pvalues2 <- rep(NA, nrow(cnv.ctrp))
for(j in 1:nrow(rna.ctrp)){
  pvalues2[j] <- min(apply(cnv.ctrp, 1, function(xx){cor.test(xx,cnv.ctrp[j,])$p.value})*ncol(cnv.ctrp), na.rm=T)
}
cnv.ctrp1 <- rna.ctrp[pvalues2*nrow(cnv.ctrp) < 0.1,]

# pre-select GEX's with 50% variations across cell lines
rna.ctrp.var <- apply(rna.ctrp, 1, var)
var.sort <- sort(rna.ctrp.var, decreasing=TRUE)
sum.var <- cumsum(var.sort)
half.var <- var.sort[which(sum.var>(sum.var[length(rna.ctrp.var)]*0.5))[1]]
rna.ctrp <- rna.ctrp[rna.ctrp.var>=half.var,]

## pre-select CNVs with explaining top 50% variations across cell lines
cnv.ctrp.var <- apply(cnv.ctrp, 1, var, na.rm=T)
var.sort <- sort(cnv.ctrp.var, decreasing=TRUE)
sum.var <- cumsum(var.sort)
half.var <- var.sort[which(sum.var>(sum.var[length(cnv.ctrp.var)]*0.5))[1]]
cnv.ctrp <- cnv.ctrp[cnv.ctrp.var>=half.var,]

## pre-select mutation features with at least one mutated/non-mutated cell lines of each tissue type
mut.census <- read.csv("Census_all2021.csv", header=T, sep=",")
mut.census <- as.character(mut.census[,1])

## select deleterious mutations based on CCLE
mut.deleterious0 <- read.csv("CCLE_mutations.csv", header=T)
mut.deleterious0 <- mut.deleterious0[,c("Hugo_Symbol","isDeleterious")]
mut.deleterious <- unique(mut.deleterious0[mut.deleterious0$isDeleterious=="True",1])

## select pathogenic mutations based on COSMIC
load("mut.pathogenic.rda")

mut.var <- intersect(mut.census, mut.pathogenic)
mut.ctrp <- mut.ctrp[rownames(mut.ctrp) %in% mut.var,]
mut.ctrp <- mut.ctrp[apply(mut.ctrp, 1, sd)>0,]

## map gene/prob ID to gene name
## https://chapmandu2.github.io/post/2016/11/25/pharmacogenetics-using-pharmacogx/
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
gene_info <- ensembldb::genes(EnsDb.Hsapiens.v75, filter=GeneNameFilter(c(''), condition="!="), return.type='data.frame') %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::mutate(probe_id=paste0(gene_id, '_at'))

map.idx <- gene_info$gene_id %in% rownames(rna.ctrp)
map.genes <- gene_info[map.idx,2:3]
map.genes <- map.genes[sort.list(map.genes$probe_id),]
rownames(rna.ctrp)[1:nrow(map.genes)] <- map.genes$gene_name
## "ENSG00000275234" "ENSG00000275342" "ENSG00000276644" "ENSG00000277443" "ENSG00000277586" "ENSG00000278291" "ENSG00000278828"
rownames(rna.ctrp)[-(1:nrow(map.genes))] <- c("AC010503.4", "PRAG1", "DACH1", "MARCKS", "NEFL", "AL161772.1", "H3C10")

save(auc.ctrp, file="cellminer_ctrp_auc_preselect.RData")
save(tissue.dum, file="cellminer_ctrp_tissue_preselect.RData")
save(rna.ctrp, file="cellminer_ctrp_rna_preselect.RData")
save(cnv.ctrp, file="cellminer_ctrp_cnv_preselect.RData")
save(mut.ctrp, file="cellminer_ctrp_mut_preselect.RData")


