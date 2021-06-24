#################################################################
## Filename: global.r
## Created: March 13, 2020
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC LUSC discovery cohort
## This file imports the underlying data and contains global functions
## and variables used by 'ui.R' and 'server.R'
#################################################################
library(pacman)
p_load(BiocManager)
p_load(scales)
p_load(gtable)
p_load(ComplexHeatmap)
p_load(RColorBrewer)
p_load(circlize)
p_load(RColorBrewer)
p_load(gplots)
p_load(WriteXLS)
p_load(grid)
p_load(bcrypt)
p_load(glue)

## import the data
load('data/data-lscc-v3.2-tumor-over-nat-ubK.RData')

## global parameters
GENE.COLUMN <<- 'geneSymbol' 
DATATYPE.COLUMN <<- 'DataType'

GENEMAX <<- 20
TITLESTRING <<- '<font size="5" face="times"><i><b>"A proteogenomic portrait of lung squamous cell carcinoma"</b></i> (<a href="https://twitter.com/shankha_ss" target="_blank_">Satpathy</a> <i>et al.</i> 2021. accepted in principle at Cell)</font><br>'
TITLESTRING.HEATMAP <<- 'Suppl. data: "A proteogenomic portrait of lung squamous cell carcinoma", Satpathy et al. Cell. 2021 - accepted in principle'

WINDOWTITLE <<- 'CPTAC-LSCC-2021'
GAPSIZEROW <<- 20
FILENAMESTRING <<- 'CPTAC-LSCC2021'

cellwidth <<- 8
cellheight <<- 10

## if TRUE, CNA data will not be z-scored
## (for thresholded CNA values)
exclude_cna_from_zscore <- FALSE

################################################
## list of gens shown at start
GENESSTART <<- c('SLK', 'PIK3CA', 'PIK3CD', 'SOX2', 'TP63', 'VIM')

##########################################
## annotation tracks shown in heatmap 
anno.all <- rev(c('Multi.omic.subtype'='NMF.cluster', 
                  'RNA.subtype.TCGA'='Subtype.TCGA.rna',
                  'TSNet.Purity'='TSNet.Purity.rna',
                  'ESTIMATE.ImmuneScore'='ESTIMATE.ImmuneScore.rna',
                  'ESTIMATE.StromalScore'='ESTIMATE.StromalScore.rna',
                  'Immune.cluster'='Immune.Cluster.rna',
                  'CIMP.status'='CIMP.status.meth',
                  'CDKN2A.mutation'='CDKN2A.mutation.status',
                  "PTEN.mutation"="PTEN.mutation.status" ,
                  "KEAP1.mutation"="KEAP1.mutation.status", 
                  "KMT2D.mutation"="KMT2D.mutation.status",
                  'TMB.WXS'='Total.Mutation.Count.per.Mb.wxs'
))

##############################
## color mappings for 'anno.all'
column.anno.col <<- list(
  Multi.omic.subtype=c('EMT enriched'='yellow', 
                       'Basal inclusive'='red', 
                       'Classical'='#FCB462', 
                       'Inflamed secretory'='#BB81B8',
                       'Proliferative primitive'='#82B2D4'),
  RNA.subtype.TCGA=c( 'Basal'='red', 
                      'Classical'='#FCB462', 
                      'Secretory'='#BB81B8',
                      'Primitive'='#82B2D4',
                      'Unknown'='grey'),
  Immune.cluster=c('Hot Tumor'='orange', 
                   'Cold Tumor'='skyblue', 
                   'NAT-enriched'='darkgreen'),
  CIMP.status=c('CIMP-high'='grey20',
                'CIMP-intermediate'='grey50',
                'CIMP-low'='grey80'
                ),

  CDKN2A.mutation=c('0'='white', '1'='black'),
  PTEN.mutation=c('0'='white', '1'='black'),
  KEAP1.mutation=c('0'='white', '1'='black'),
  KMT2D.mutation=c('0'='white', '1'='black'),
  
  
  Smoking.Score.WGS=circlize::colorRamp2( c(0, 0.5, 1), c('grey90', 'grey50','black')),
  TSNet.Purity=circlize::colorRamp2( c(0, 0.5, 1), c('grey90', 'grey50','black')),
  ESTIMATE.ImmuneScore=circlize::colorRamp2( c(676, 4700, 6800), 
                                             c(rgb(255,245,240, maxColorValue = 255),  
                                               rgb(251,106,74, maxColorValue = 255),  
                                               rgb(165,21,22, maxColorValue = 255))),
  ESTIMATE.StromalScore=circlize::colorRamp2( c(-1388, 1854, 4157), 
                                             c(rgb(255,245,240, maxColorValue = 255),  
                                               rgb(251,106,74, maxColorValue = 255),  
                                               rgb(165,21,22, maxColorValue = 255))),
  
  TMB.WXS=circlize::colorRamp2( c(0.4, 9, 24), c('grey90', 'grey50','black'))
)


#############################
## columns used for sorting
columns.to.sort <- rev(anno.all)


##################################################################
## 21060613 bcrypt
authenticateUser <- function(passphrase){
  if(nchar(as.character(passphrase)) > 0){
    return(checkpw(as.character(passphrase), "$2a$12$yJXCOp04vsk4SoUFpKHLauJ31qpOiiziFDbUVwyjU/uNJSPOgfY.G"))
  } else {
    return(FALSE)
  }
}


###################################################################
## heatmap using ComplexHeatmap package
MyComplexHeatmap <- function(m, rdesc, cdesc, cdesc.color, max.val, column2sort, main=''){
  
  if(is.null(max.val)){
    m.max <- ceiling(max(abs(m), na.rm=T))
  } else {
    m.max <- max.val
    m[m > m.max] <- m.max
    m[m < -m.max] <- -m.max
  }
  cdesc$TMB.WXS <- as.numeric(cdesc$TMB.WXS)
  
  ## #####################################
  ## colorscale for heatmap
  col.hm <- colorRamp2(seq(-m.max, m.max, length.out=11), rev(brewer.pal (11, "RdBu")))
  
  #########################################
  ## annotations
  cdesc.ha <- HeatmapAnnotation(df=cdesc, 
                                col=cdesc.color,
                                show_legend = T, show_annotation_name = T, 
                                annotation_name_side = 'left',
                               # na_col = "white",
                                annotation_legend_param=list(
                                  direction='horizontal'
                                  )
  )
  ####################################
  ## heatmap
  ## determine height
  n.entries <- nrow(m)
  n.genes <- length(unique(rdesc$geneSymbol))

  cat('hm height:', dynamicHeightHM(n.entries, n.genes), '\n')

  cat(column2sort,'\n')
  
  hm <- Heatmap(m, col=col.hm,
                
                cluster_columns = F,
                cluster_rows = F,
                
                top_annotation = cdesc.ha,
                
                row_split = rdesc$geneSymbol, 
                row_title_rot=0,
                #column_split=column2sort,
                
                name='relative abundance',
                show_row_names = T,
                show_column_names = F,
                #use_raster = FALSE,
                
                height = unit(0.5 ,'cm') * n.entries,
                
                column_title = main
                )
  ## plot
  draw(hm, annotation_legend_side='bottom')
}


##################################################################
## function to extract gene names from a string of
## character
extractGenes <- function(genes.char){

    gene.max=GENEMAX

    if(is.null(genes.char))
      return(NULL)

    if( length(genes.char) == 0 ){
        return(NULL)
    }
    ## extract genes
    genes.vec= unlist(strsplit(genes.char, ','))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ' '))
    if(length(genes.vec)==1)
        genes.vec=unlist(strsplit(genes.char, ';'))

    ## unique gene names
    genes.vec <- unique(genes.vec)

    ## limit to 'gene.max' genes
    if(length(genes.vec) > gene.max){
        warning(paste('more than', gene.max,'gene ids submitted! Showing results for the first 20 genes in the list.\n'))
        genes.vec <- genes.vec[1:gene.max]
    }
    return(genes.vec)
}

##################################################################
## function to dynamically determine the height (in px) of the heatmap
## depending on the number of genes
dynamicHeightHM <- function(n.entries, n.genes){
  
 height <- 0.3937*n.entries + 0.7423 ## height in inch
  #height <- 5*n.entries + 18.85
  height <- height + 12*0.3937 + 0.7423            ## add annotation tracks
  height <- height * 48             ## inch  to pixel
  
  
  return(height)
}


#######################################################
## find all entries with associated gene name in
## the dataset. returns vector of indices.
findGenesInDataset <- function(gene,                                 ## vector of gene symbols 
                               show.sites=c('most variable', 'all')  
                               ){
  
  show.sites <- match.arg(show.sites)
  
  ## remove spaces
  gene <- gsub(' ', '', gene )
  gene <- unique(gene)
  
  ## remove empty strings
  gene.nchar=which(nchar(gene) == 0)
  if(length(gene.nchar) > 0)
    gene <- gene[-gene.nchar]
  
  if(length(gene) == 0) return()
  
  ## check whether the genes are present in the dataset
  gene.idx <- grep( paste(paste('(^|,)', gene, '($|,)', sep=''), collapse='|'), gsub(' ', '', row.anno[, GENE.COLUMN]) )
  if( length(gene.idx) == 0 ){
    stop('None of the gene ids you have entered could be found in the dataset!\n')
  }

  ## use row names
  gene.idx <- rownames(tab.expr.all)[gene.idx]
  
  ## extract data and remove empty rows
  data.tmp <- tab.expr.all[gene.idx, ]
  row.anno.tmp <- row.anno[gene.idx, ]
  
  ###################################
  ## most variable site
  if(show.sites=='most variable'){
    
    ## extract MS phospho
    vm.idx <- grep('_pSTY', row.anno.tmp[, DATATYPE.COLUMN])
    if( length(vm.idx) > 0 ){
      
      vm.sd <- apply(data.tmp[ vm.idx, ], 1, sd, na.rm=T)
      
      rm.idx <- tapply(vm.sd, row.anno.tmp[vm.idx, GENE.COLUMN],  function(x) names(x)[which.max(x)])
      rm.idx <- setdiff( names(vm.sd), unlist(rm.idx) )
      
      gene.idx <- setdiff(gene.idx, rm.idx)
      
      data.tmp <- data.tmp[gene.idx, ]
      row.anno.tmp <- row.anno.tmp[gene.idx, ]
    }
    
    ## extract MS acetyl
    vm.idx <- grep('_acK', row.anno.tmp[, DATATYPE.COLUMN])
    if( length(vm.idx) > 0 ){
      
      vm.sd <- apply(data.tmp[ vm.idx, ], 1, sd, na.rm=T)
      
      rm.idx <- tapply(vm.sd, row.anno.tmp[vm.idx, GENE.COLUMN],  function(x) names(x)[which.max(x)])
      rm.idx <- setdiff( names(vm.sd), unlist(rm.idx) )
      
      gene.idx <- setdiff(gene.idx, rm.idx)
      
      data.tmp <- data.tmp[gene.idx, ]
      row.anno.tmp <- row.anno.tmp[gene.idx, ]
    }
    
    ## extract MS ubiquityl
    vm.idx <- grep('_ubK', row.anno.tmp[, DATATYPE.COLUMN])
    if( length(vm.idx) > 0 ){
      
      vm.sd <- apply(data.tmp[ vm.idx, ], 1, sd, na.rm=T)
      
      rm.idx <- tapply(vm.sd, row.anno.tmp[vm.idx, GENE.COLUMN],  function(x) names(x)[which.max(x)])
      rm.idx <- setdiff( names(vm.sd), unlist(rm.idx) )
      
      gene.idx <- setdiff(gene.idx, rm.idx)
      
      data.tmp <- data.tmp[gene.idx, ]
      row.anno.tmp <- row.anno.tmp[gene.idx, ]
    }
    
    
  }

  return(gene.idx)
}


#################################################################
## draw the actual heatmap
##
#################################################################
makeHM <- function(gene, filename=NA, expr=tab.expr.all, 
                   column.anno=column.anno, row.anno=row.anno, zscore="none", 
                   anno.class='PAM50', sort.dir, 
                   show.sites='all', min.val=-3, max.val=3, 
                   main='', ...){

    n.bins=12

    ## #############################
    
    ## reorder
    cat(sort.dir, '\n')
    ord.idx <- order(column.anno[, anno.class], decreasing = ifelse(sort.dir == 'descending', T, F))
    #ord.idx <- with(column.anno, order( eval(parse( text=anno.class )), decreasing = ifelse(sort.dir == 'descending', T, F) ))
    column.anno <- column.anno[ord.idx,]
    expr <- expr[, ord.idx]
    
    ################################
    ## find genes
    gene.idx <- findGenesInDataset(gene, show.sites) 
    
    #################################
    ## extract sample ids
    sampleIDs <- colnames(expr)

    #################################
    ## extract genes of interest
    expr.select <- expr[gene.idx, ]
    row.anno.select <- row.anno[gene.idx, ]
    
    #####################################################
    ## labels for the rows in the heatmap
    featureIDs.anno.select <- paste(row.anno.select[ , 'geneSymbol'],
                                gsub('6_ubK', 'ubK' ,    
                                    gsub( '5_acK', 'acK',
                                      gsub( '4_pSTY', 'pSTY', 
                                          gsub('3_Protein', 'Protein', 
                                               gsub('2_RNAseq', 'RNA-Seq', 
                                                    gsub('1_CNA', 'CNA', row.anno.select[ , 'DataType'])
                                                    ) 
                                               )
                                          )
                                       )
                                    )
                                  )

    #####################################################
    ## add phosphosite annotation
    #####################################################

    ## MS
    ms.psty.idx <- grep('pSTY', featureIDs.anno.select)
    if(length(ms.psty.idx)>0){
        featureIDs.anno.select[ ms.psty.idx ] <- paste( sub('pSTY', '', 
                                                            featureIDs.anno.select[ ms.psty.idx ]), 
                                                        paste('p', sub('.*_([S|T|Y][0-9]*)[s|t|y].*', '\\1', row.anno.select[ ms.psty.idx, 'ID']), sep=''), sep='' )
    }
    ms.ack.idx <- grep('acK', featureIDs.anno.select)
    if(length(ms.ack.idx)>0){
      featureIDs.anno.select[ ms.ack.idx ] <- paste( sub('acK', '', 
                                                         featureIDs.anno.select[ ms.ack.idx ]), 
                                                     paste('ac', sub('.*_([K][0-9]*)[k].*', '\\1', row.anno.select[ ms.ack.idx, 'ID']), sep=''), sep='' )
    }
    ms.ubk.idx <- grep('ubK', featureIDs.anno.select)
    if(length(ms.ubk.idx)>0){
      featureIDs.anno.select[ ms.ubk.idx ] <- paste( sub('ubK', '', 
                                                         featureIDs.anno.select[ ms.ubk.idx ]), 
                                                     paste('ub', sub('.*_([K][0-9]*)[k].*', '\\1', row.anno.select[ ms.ubk.idx, 'ID']), sep=''), sep='' )
    }
    
    
    #################################
    ## apply zscore
    rownames(expr.select) <- featureIDs.anno.select

    if(zscore == "row"){
      
      if(exclude_cna_from_zscore){
        ## exclude CNA data from Z-scoring
        expr.select.zscore.tmp <- lapply( rownames(expr.select), function(xx){x=expr.select[xx,];
            if( length(grep( 'CNA', xx)) == 0)return((x-mean(x, na.rm=T))/sd(x, na.rm=T));
            if( length(grep( 'CNA', xx)) > 0)return(x);
        })
       
      } else {
        expr.select.zscore.tmp <- apply(expr.select, 1, function(x)(x-mean(x, na.rm=T))/sd(x, na.rm=T))
      }
      
      expr.select.zscore <- matrix(unlist(expr.select.zscore.tmp), ncol=ncol(expr.select), byrow=T, dimnames=list(rownames(expr.select), colnames(expr.select)))
      
    } else {
        expr.select.zscore <- expr.select
    }

    ## cap at -3/3
    expr.select.zscore[which(expr.select.zscore < min.val)] <- min.val
    expr.select.zscore[which(expr.select.zscore > max.val)] <- max.val

    ##############################
    ## column annotation
    column.anno.fig <- column.anno[, anno.all]
    colnames(column.anno.fig) <- names(anno.all)

    
    ################################
    ## gaps
    ################################
    #if(sort.dir == 'descending')
    #  gaps.column=cumsum(c( rev(table(column.anno[, anno.class])) ))
    #if(sort.dir == 'ascending')
    #  gaps.column=cumsum(c(  table(column.anno[, anno.class]) ))
    #gaps.row=cumsum(table(sub(' .*', '', featureIDs.anno.select)))


    ################################
    ## colors misc
    #color.breaks = seq(min.val, max.val, length.out=n.bins)
    #color.hm =  colorRampPalette( c('blue', 'grey', 'red'))(length(color.breaks))
    #color.border = 'white'

    #legend_breaks=seq(-3, 3, 1)
    #legend_labels=c('-3               ', '-2', '-1', ' 0', '+1', '+2' ,'+3')


    ###############################
    ## heatmap
    column2sort <- column.anno[, anno.class]
    ## disable column split if sorting vector is numeric
    if(mode(column.anno.col[[names(anno.class)]]) == 'function')
      column2sort <- NULL
    
    save(expr.select.zscore, row.anno.select, column.anno.fig, column.anno.col,  featureIDs.anno.select, ms.ubk.idx,
          max.val, column2sort, anno.class, file='debug.RData')
        
    MyComplexHeatmap(data.matrix(expr.select.zscore), 
                     row.anno.select, column.anno.fig, column.anno.col, 
                     max.val=max.val, column2sort=column2sort,
                     main=main
                     )
    

    #########################################################################################
    ##
    ## - return part of table that is shown in the heatmap and that can be downloaded
    ## - change the CNA values (-3, -1, 0, 1, 3) back to the original values (-1, -.3, 0, .3, 1)
    ##
    #########################################################################################
    # cna.idx <- grep('CNA', rownames(expr.select))
    # if(length(cna.idx) > 0){
    #     expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == -3 ] <- -1
    #     expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == -1 ] <- -.3
    #     expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == 1 ] <- .3
    #     expr.select[ cna.idx, ][expr.select[ cna.idx, ]  == 3 ] <- 1
    # }
    # 
    ## add row annotation
    mat.row <- as.data.frame( cbind( row.anno[gene.idx, ], expr.select, deparse.level=0 ), stringsAsFactors=F )
    rownames(mat.row) <- make.names(featureIDs.anno.select, unique = T)

    ## column annotation
    mat.col <- as.data.frame( cbind(matrix('', nrow=ncol(column.anno.fig), ncol=ncol(row.anno)), t(column.anno.fig), deparse.level=0), stringsAsFactors=F)
    colnames(mat.col) <- colnames(mat.row)

    ## put everything together
    mat <- rbind( mat.col, mat.row)

    return(mat)
}
