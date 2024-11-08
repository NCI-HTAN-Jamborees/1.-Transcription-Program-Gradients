# Usage: enrich.gsets(fg, gsets, bg, nc=1L, overlap.cutoff=0, padj.cutoff=1.1)
# * fg: a vector of genes
# * gsets: an R list of gene sets with names, each element of the list is a vector of genes within the gene set
# * bg: a vector of "background" genes (the "universe")
# * nc: number of cores for parallelization
# * overlap.cutoff: only look at gene sets with number of overlapping genes > this number, default 0
# * padj.cutoff: return enrichment results with BH-adjusted P < this number, default 1.1 means return all
# * the only thing you need to carefully define when using this is the bg -- fg is trivial, gsets can be just default ones downloaded from databases
# code by Kuoyuan Cheng


library(parallel)
library(data.table)

enrich.gsets <- function(fg, gsets, bg, nc=1L, overlap.cutoff=0, padj.cutoff=1.1) {
  # the over-representation type of enrichment test (with Fisher's exact test here.)
  # fg: query genes; gsets: gene sets as a list object; bg: background genes; overlap.cutoff: only select those gene sets with >this value overlap with fg genes; padj.cutoff: fdr threshold.
  # nc: number of cores

  fg1 <- intersect(fg, bg)
  tmp <- sapply(gsets, function(x) sum(x %in% fg1))
  gsets <- gsets[tmp>overlap.cutoff]

  if (length(fg)==0) {
    warning("The number of query genes is zero, NULL returned.\n")
    return(NULL)
  }

  enrich.gset0 <- function(fg, gset, bg) {
    res <- enrich.test(qset=fg, refset=gset, uset=bg, alternative="greater")
    data.table(overlap.size=res$table[1,1], gene.set.size=res$table[1,1]+res$table[2,1], odds.ratio=res$estimate, pval=res$p.value, overlap.genes=list(unique(intersect(intersect(fg, gset),bg))))
  }
  res <- mclapply(gsets, enrich.gset0, fg=fg, bg=bg, mc.cores=nc)
  res <- rbindlist(res, idcol="gene.set")
  if (ncol(res)==0) return(NULL)
  res[, padj:=p.adjust(pval, method="BH")]
  res <- res[order(padj,pval)][padj<padj.cutoff]
  setcolorder(res, c("gene.set","odds.ratio","pval","padj","gene.set.size","overlap.size","overlap.genes"))
  res
}


enrich.test <- function(qset=NULL, refset=NULL, uset=NULL, confus.mat=NULL, ...) {

  # Fisher's exact test of enrichment of a reference set 
  # (or actual positive set, refset) in a query set 
  # (or predicted positive set, qset), with the background being the 
  # universal set (uset), or from a given confusion matrix as confus.mat.
  # ... passed to fisher.test()

  if (is.null(confus.mat)) {
    conf <- make.confus.mat(qset, refset, uset, margins=FALSE)
  } else conf <- confus.mat[1:2, 1:2]

  # fisher's exact test
  res <- fisher.test(conf, ...)
  res$table <- conf
  return(res)
}

make.confus.mat <- function(qset, refset, uset, margins=TRUE) {

  # make a confusion matrix of TP/FP/TN/FN given a reference set (or actual positive set, refset) and a query set (or predicted positive set, qset), with the background being the universal set (uset).
  # set margins=TRUE to add margins

  # if qset is empty or NA, return NA
  if (length(qset)==0) {
    warning("In make.confus.mat: qset has zero length. NA returned.\n")
    return(NA)
  } else if (length(qset)==1 && is.na(qset)) {
    warning("In make.confus.mat: qset is NA. NA returned.\n")
    return(NA)
  }
  # make sure uset has unique items
  uset <- unique(uset)
  # logical vector denoting whether each item is in qset or not
  qsetl <- factor(uset %in% qset, levels=c(TRUE, FALSE), labels=c("Positive","Negative"))
  # logical vector denoting whether each item is in refset or not
  refsetl <- factor(uset %in% refset, levels=c(TRUE, FALSE), labels=c("Positive","Negative"))
  # a confusion matrix as a table
  res <- table(`Query/Prediction`=qsetl, `Reference/Actual`=refsetl)
  # add margins
  if (margins) res <- addmargins(res)
  # return
  return(res)
}
