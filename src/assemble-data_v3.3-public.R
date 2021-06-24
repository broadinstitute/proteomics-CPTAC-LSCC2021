options(stringsAsFactors = F)
rm(list=ls())
#library(pacman)
library(cmapR)
library(dplyr)
library(glue)

setwd('c:/Users/karsten/Dropbox/Devel/CPTAC-LSCC2021/')

## import gct files
data.dir <- 'g:/Shared drives/CPTAC3.0_LSCC_Data/lscc-v3.3-public-data-freeze/'
label <- 'v3.3-public-tumor-over-nat'

ord.column <- 'NMF.cluster'


gct.str <- c(glue("{data.dir}lscc-v3.3-public-acetylome-ratio-norm-NArm.gct"),
             glue("{data.dir}lscc-v3.3-public-gene-level-cnv-gistic2-log-ratio-all_data_by_genes.gct"),
             glue("{data.dir}lscc-v3.3-public-phosphoproteome-ratio-norm-NArm.gct"),
             glue("{data.dir}lscc-v3.3-public-proteome-ratio-norm-NArm.gct"),
             glue("{data.dir}lscc-v3.3-public-rnaseq-uq-fpkm-log2-NArm.gct"),
             glue("{data.dir}lscc-v3.3-public-ubiquitylome-batch-corrected-ratio-norm-NArm.gct")
             )

gct <- lapply(gct.str, parse_gctx)
names(gct) <- c('5_acK', '1_CNA', '4_pSTY', '3_Protein', '2_RNAseq', '6_ubK')


###########################################
## calculate tumor-normal ratios
gct.expr <- vector('list', length(gct))
names(gct.expr) <- names(gct)
gct.rdesc <- gct.cdesc <- gct.rid <- gct.cid <- gct.expr
  
for(i in names(gct)){
  cdesc <- gct[[i]]@cdesc
  rdesc <- gct[[i]]@rdesc
  rid <- gct[[i]]@rid %>% gsub(' ', '', .)
  cid <- gct[[i]]@cid
  mat <- gct[[i]]@mat
  
  if('NAT' %in% cdesc$Type){
    nat.idx <- which(cdesc$Type == 'NAT')
    tumor.idx <- which(cdesc$Type == 'Tumor')
    
    ## separate tumor and nat
    mat.n <- mat[, nat.idx]
    mat.t <- mat[, tumor.idx]
    
    ## - exclude tumors without nat
    ## - bring into same order
    n.cid <- colnames(mat.n)
    t.cid <- colnames(mat.t)
    
    common <- intersect( t.cid, sub('\\.N$','', n.cid) )
    t.cid <- common
    n.cid <- paste0(common, '.N') 
      
    mat.t <- mat.t[, t.cid]
    mat.n <- mat.n[, n.cid]
    
    ## tumor-nat
    mat <- mat.t - mat.n
    cid <- t.cid
    cdesc <- cdesc[cid, ]
  }
  
  gct.expr[[i]] <- mat
  gct.rdesc[[i]] <- rdesc
  gct.cdesc[[i]] <- cdesc
  gct.rid[[i]] <- rid
  gct.cid[[i]] <- cid
}

###############################################################
##  - samples to include - all samples with T/N pairs
##  - use the proteome T/N table as reference
samp.common <- gct.cid[[grep('_Protein', names(gct.cid))]] 

## number of samples (T/N-pairs) per data table
gct.N.samp <- sapply(gct.cid, function(x) sum(x %in% samp.common))

###################
## rdesc 
gct.gene <- lapply(gct.rdesc, function(x) try(x[ ,grep('geneSymbol', colnames(x), value=T, ignore.case = T)[1]], silent=T))
if(sum(sapply(gct.gene, length) < 2) > 0){
  idx <- which( sapply(gct.gene, length) < 2 )
  for(i in idx)
    gct.gene[[i]] <- gct.rid[[i]]
}

## data type
gct.data.type <- lapply(names(gct.gene), function(x) rep(x, nrow(gct.expr[[x]])))
names(gct.data.type) <- names(gct.gene)

## ids
gct.id <- lapply(gct.rdesc, function(x) gsub(' ', '', x[, 'id']) )

if(sum(sapply(gct.id, length) < 2) > 0){
  idx <- which( sapply(gct.id, length) < 2 )
  for(i in idx)
    gct.id[[i]] <- gct.rid[[i]]
}

## combine
gct.rdesc <- lapply(names(gct.gene), function(x) data.frame(ID=gct.id[[x]], geneSymbol=gct.gene[[x]], DataType=gct.data.type[[x]]))

rdesc <- Reduce(f = rbind, gct.rdesc)
row.anno <- rdesc

###########################################
##               sample annotation
## use cdesc from proteome
cdesc <- gct[['3_Protein']]@cdesc


## convert data types
cdesc$Smoking.score.wxs <- as.numeric(cdesc$Smoking.score.wxs)
cdesc$ESTIMATE.ImmuneScore.rna <- as.numeric(cdesc$ESTIMATE.ImmuneScore.rna)
cdesc$TSNet.Purity.rna <- as.numeric(cdesc$TSNet.Purity.rna)
cdesc$Total.Mutation.Count.per.Mb.wxs <- as.numeric(cdesc$Total.Mutation.Count.per.Mb.wxs)
cdesc$ESTIMATE.StromalScore.rna <- as.numeric(cdesc$ESTIMATE.StromalScore.rna)

#################################
## single data frame
column.anno <- cdesc[samp.common, ]

for(i in names(gct.expr)){
  
  if(i == names(gct.expr)[1]) {
    
    tab.expr.all <- gct.expr[[i]][, samp.common]
    
  } else {
    
    tab_to_add <- gct.expr[[i]]
    
    if( sum( !samp.common %in% colnames(tab_to_add)) > 0) {
      
      samp_missing <- setdiff(samp.common, colnames(tab_to_add))
      dummy <- matrix(NA, nrow=nrow(tab_to_add), ncol=length(samp_missing), dimnames = list(rownames(tab_to_add), samp_missing))

      tab_to_add <- cbind(tab_to_add, dummy)
      tab_to_add <- tab_to_add[, samp.common]
      
    } else{
      tab_to_add <- tab_to_add[, samp.common]  
    }
    
    tab.expr.all <- rbind(tab.expr.all, tab_to_add)
  }
}

#tab.expr.all <- Reduce(f=rbind, gct.expr)
rownames(tab.expr.all) <- make.unique(rownames(tab.expr.all)) 
rownames(row.anno) <- rownames(tab.expr.all)

## reorder
ord.idx <- order(column.anno[, ord.column])
column.anno <- column.anno[ord.idx,]
tab.expr.all <- tab.expr.all[, ord.idx]

## reorder rows
row.idx <- with(row.anno, order(geneSymbol, DataType))
tab.expr.all <- tab.expr.all[row.idx, ]
row.anno <- row.anno[row.idx, ]

## export
save(column.anno, row.anno, tab.expr.all, gct.N.samp, file = glue('data/data-lscc-{label}.RData'))

