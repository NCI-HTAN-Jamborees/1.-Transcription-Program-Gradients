source("./enrich.funcs.R")

share_NMFs <- readRDS("./shared_NMFs.rds")

# gene <- c("VSIR", "TIGIT", "ICOS", "EOMES", "HAVCR2", "PDCD1", "BTLA", "CD244", "LAG3", "CD160", "CTLA4", "CD96")

over_rep_basal <- enrichGO(gene = share_NMFs$shared_basal,
                                   universe = names(rownames(nmf_info$feature.loadings)),
                                   OrgDb = organism, 
                                   keyType = 'SYMBOL',
                                   readable = T,
                                   ont = "BP",
                                   pvalueCutoff = 0.05, 
                                   qvalueCutoff = 0.10)

over_rep_luminal <- enrichGO(gene = share_NMFs$shared_luminal,
                                   universe = names(rownames(nmf_info$feature.loadings)),
                                   OrgDb = organism, 
                                   keyType = 'SYMBOL',
                                   readable = T,
                                   ont = "BP",
                                   pvalueCutoff = 0.05, 
                                   qvalueCutoff = 0.10)

over_rep_basal@result
