setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/1.calculate the result")
metadata<-read.csv("metadata.final2.csv",row.names = 1)
table<-read.csv("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/table/table.nolichenandhost.csv",row.names = 1)
fungi.ind<-table
fungi.ind<- ifelse(fungi.ind>0,1,0)
fungi.ind<-as.data.frame(t(fungi.ind))
fungi.ind$sample_names<-row.names(fungi.ind)
fungi2<-merge(fungi.ind,metadata,by="sample_names",all=T)
row.names(fungi2)<-fungi2$sample_names
fungi2<-fungi2[,-1]
fungi2<-fungi2[!is.na(fungi2$e57bceb79df4a766f8fb9d049798c662ba89dc9e),]
###
LV<-fungi2[fungi2$site %in% "LV",]
LV <- LV[, c(colSums(LV[, -(ncol(LV)-1):-(ncol(LV))]) > 2, TRUE, TRUE)]
LV <-LV [,-ncol(LV)]
LW<-fungi2[fungi2$site %in% "LW",]
LW <- LW[, c(colSums(LW[, -(ncol(LW)-1):-(ncol(LW))]) > 2, TRUE, TRUE)]
LW <-LW [,-ncol(LW)]
LZ<-fungi2[fungi2$site %in% "LZ",]
LZ <- LZ[, c(colSums(LZ[, -(ncol(LZ)-1):-(ncol(LZ))]) > 2, TRUE, TRUE)]
LZ <-LZ [,-ncol(LZ)]
###
LZ.ind<-list(LZ)
LV.ind<-list(LV)
LW.ind<-list(LW)
calc_phi <- function(x, perm = 10000){
  setDF(x)
  res <- list()
  res$PHI <- try( multipatt(
    x = x[, -ncol(x)],
    cluster = x[,ncol(x)],
    func = "r.g",
    control = how(nperm = perm),
    duleg = TRUE) )
  ## Calculate complementary fidelity values for the tree species (UbinB metric)
  res$III <- try( multipatt(
    x = x[, -ncol(x)],
    cluster = x[,ncol(x)],
    func = "IndVal.g",
    control = how(nperm = perm),
    duleg = TRUE) )
  
  return(res)
}

library(plyr)
library(doFuture)
library(data.table)
library(indicspecies)
registerDoFuture()
plan(multisession, workers = 4)      # for RStudio
options(future.globals.maxSize = 5e9)  # 5GB; default = 500 * 1024 ^ 2 = 500 MiB

LZ.PIRES <- llply(.data = LZ.ind, .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
LV.PIRES <- llply(.data = LV.ind, .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
LW.PIRES <- llply(.data = LW.ind, .fun = function(x){
  calc_phi(x)
}, .parallel = TRUE)
##indicators
##using phi to analysis
load("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi/phi.separate.sites.RData")
###extract phi
library(metagMisc)
extract<-function(fungi2.PIRES,fungi2){
  
  ## Extract indvals
  IIIs <- llply(.data = fungi2.PIRES, .fun = function(x){
    
    iii <- metagMisc::dfRowName(x = x$III$sign, name = "OTU")
    setDT(iii)
    
    treesp <- colnames(iii)
    treesp <- treesp[ ! treesp %in% c("OTU", "index", "stat", "p.value") ]
    
    iii[, BestHost := treesp[ index ] ]
    setcolorder(x = iii, neworder = c("OTU", "BestHost", "index", "stat", "p.value"))
    
    return(iii)
  })
  ## Extract summary in tabular form (based on `summary.multipatt`)
  get_stats_by_cluster <- function(x,
                                   alpha = 0.05, minstat = NULL, At = NULL, Bt = NULL, indvalcomp = FALSE, ...) {
    
    ncomb = ncol(x$str)
    ncolsign = ncol(x$sign)
    nsps = nrow(x$sign)
    
    cat("\n Multilevel pattern analysis")
    cat("\n ---------------------------\n")
    cat("\n Association function:", x$func)
    cat("\n Significance level (alpha):", alpha)
    if (x$func == "IndVal" || x$func == "IndVal.g") {
      if (!is.null(At)) { cat("\n Minimum positive predictive value (At):", At) }
      if (!is.null(Bt)) { cat("\n Minimum sensitivity (Bt):", Bt) }
    }
    cat("\n\n Total number of species:", nsps)
    
    sel = !is.na(x$sign$p.value) & x$sign$p.value <= alpha
    if (!is.null(minstat)) { sel = sel & (x$sign$stat >= minstat) }
    if (!is.null(Bt) && !is.null(x$B)) {
      for (i in 1:nrow(x$sign)){ sel[i] = sel[i] && (x$B[i, x$sign$index[i]] >= Bt) }
    }
    if (!is.null(At) && !is.null(x$A)) {
      for (i in 1:nrow(x$sign)) { sel[i] = sel[i] && (x$A[i, x$sign$index[i]] >= At) }
    }
    
    a = x$sign[sel, ]
    
    cat("\n Selected number of species:", nrow(a), "\n")
    
    cols = (ncolsign - 1):ncolsign
    
    if (indvalcomp && !is.null(x$B) && !is.null(x$A)) {
      As = numeric(nrow(x$sign))
      Bs = numeric(nrow(x$sign))
      for (i in 1:nrow(x$sign)) {
        As[i] = x$A[i, x$sign$index[i]]
        Bs[i] = x$B[i, x$sign$index[i]]
      }
      y = cbind(x$sign, As, Bs)
      cols = c(ncol(y) - 1, ncol(y), cols)
      names(y) = c(names(x$sign), "A", "B")
    }
    else y = x$sign
    
    for (k in 1:(ncolsign - 4)) {
      cat(" Number of species associated to", k, if (k == 1) 
        "group:"
        else "groups:", sum(rowSums(a[, 1:(ncolsign - 3)]) == k), "\n")
    }
    
    
    res <- list()
    
    for (i in 1:ncomb) {
      sel = x$sign$index == i & !is.na(x$sign$p.value) & x$sign$p.value <= alpha
      if (!is.null(minstat)) { sel = sel & (x$sign$stat >= minstat) }
      if (!is.null(Bt) && !is.null(x$B)) {
        for (j in 1:nrow(x$sign)) { sel[j] = sel[j] && (x$B[j, x$sign$index[j]] >= Bt) }
      }
      if (!is.null(At) && !is.null(x$A)) {
        for (j in 1:nrow(x$sign)) { sel[j] = sel[j] && (x$A[j, x$sign$index[j]] >= At) }
      }
      m = y[sel, ]
      if (nrow(m) > 0) {
        m <- m[ order(m$stat, decreasing = TRUE), cols ]
        m <- metagMisc::dfRowName(x = m, name = "OTU")
        setDT(m)
        res[[ colnames(x$comb)[i] ]] <- m
      }
    }
    return(res)
  }
  IIS <- llply(.data = fungi2.PIRES,
               .fun = function(x){ 
                 get_stats_by_cluster(x$III, alpha = 1, indvalcomp = TRUE)
               })
  
  library(ggplot2)
  x <- fungi2
  occ <- data.table(Sp = colnames(x)[-ncol(x)], Occ = colSums(x[, -ncol(x)]))
  r <- IIIs[[1]]
  r <- merge(r, occ, by.x = "OTU", by.y = "Sp")
  r[, invp := 1/p.value]
  r[, sign := p.value < 0.05 ]
  ##non-separated
  r$BestHost<-gsub("^s.","",r$BestHost)
  p<-ggplot(data = r, aes(x = Occ, y = stat)) + 
    geom_point(aes(size = invp, color = sign,alpha=0.5),shape=1) + 
    facet_wrap(~ BestHost)+
    labs(x="Occ",y="phi value")+
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  ###
  p2<-ggplot(data = r, aes(x = p.value, y = stat)) + #Vasula[Occ < 8 & p.value < 0.5 ]
    geom_vline(xintercept=0.05, color="darkgrey", linetype = "longdash") + 
    geom_point(aes(size = Occ, color = Occ), alpha = 0.5) + 
    scale_size(range = c(0, 6)) +
    scale_color_distiller(direction =0.5 ) + 
    scale_shape_manual(values = c(10, 18)) +
    facet_wrap(~ BestHost, nrow = 3) + 
    scale_y_continuous(breaks=c(0.1, 0.3, 0.5,0.7,0.9)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
    labs(x = "p-value", y = "Phi", color = "Number of\noccurrences", size = "Number of\noccurrences") +
    labs(x="p-value",y="phi value")+
    theme_light() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  return(list(p,p2,r))  
  
}
library(plyr)
library(data.table)
LW.result<-extract(LW.PIRES,LW)
LZ.result<-extract(LZ.PIRES,LZ)
LV.result<-extract(LV.PIRES,LV)
LV.result[[1]]
LV.result[[2]]
##growth forms
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/third/")
tax<-read.csv("tax.final.csv",row.names = 1)
tax$Phylum[is.na(tax$Phylum)]<-"p__unclassified"
tax$Class[is.na(tax$Class)]<-"c__unclassified"
tax$Order[is.na(tax$Order)]<-"o__unclassified"
tax$Family[is.na(tax$Family)]<-"f__unclassified"
tax$Genus[is.na(tax$Genus)]<-"g__unclassified"
tax$Species[is.na(tax$Species)]<-"s__unclassified"
tax2<-tax
tax2$Genus<-gsub(".*_","",tax2$Genus)
##
r2<-merge(LZ.result[[3]],tax2,by.x="OTU",by.y="qseqid",all=T)
fungi2.r2<-r2[!is.na(r2$BestHost),]
setwd("C:/Users/meirong/Desktop/PhD project/second preference/final.table/organised code/1.phi/")
LZ.r2<-fungi2.r2
LZ.r<-fungi2.r2[fungi2.r2$p.value<0.05,]
write.csv(LZ.r,"LZ.r2.tax.csv")

r2<-merge(LW.result[[3]],tax2,by.x="OTU",by.y="qseqid",all=T)
fungi2.r2<-r2[!is.na(r2$BestHost),]

LW.r2<-fungi2.r2
LW.r<-fungi2.r2[fungi2.r2$p.value<0.05,]
write.csv(LW.r,"LW.r2.tax.csv")

r2<-merge(LV.result[[3]],tax2,by.x="OTU",by.y="qseqid",all=T)
fungi2.r2<-r2[!is.na(r2$BestHost),]

LV.r2<-fungi2.r2
LV.r<-fungi2.r2[fungi2.r2$p.value<0.05,]
write.csv(LV.r,"LV.r2.tax.csv")

