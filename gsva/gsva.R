library(GSVA)
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(msigdbr)
library(limma)
library(edgeR)
library(pheatmap)


exp <- read.table('exp.txt')


gmt <- read.gmt('./msigdb.v7.5.1.entrez.gmt')


gmt_BP <- gmt[grep('^GOBP', gmt$term),]
gmt_BP <- split(gmt_BP$gene, gmt_BP$term, drop = TRUE)


#### msigdbr
# GO_df_all <- msigdbr(species = "Homo sapiens",
#                     category = "C5")  
# GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
# GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
# go_list <- split(GO_df$gene_symbol, GO_df$gs_name) ##Based on gs_name, classify gene_symbol



## Clean the expression matrix
exp <- rownames_to_column(exp, 'ensembl')
gene_info <- bitr(exp$ensembl, fromType = 'ENSEMBL', toType = 'ENTREZID',
                  OrgDb = 'org.Hs.eg.db')
exp_new <- data.frame(gene_info, exp[match(gene_info$ENSEMBL, exp$ensembl),])
exp_new <- exp_new[,-c(1,3)]
exp_new <- exp_new[!duplicated(exp_new$ENTREZID),]
rownames(exp_new) <- exp_new$ENTREZID
exp_new <- exp_new[,-1]


# gsva
gsva_gobp <- gsva(expr = as.matrix(exp_new), gmt_BP,
                  kcdf = "Poisson",
                  verbose = T)


# limma
group_info <- factor(c(rep('R',10), rep('S', 9)))
design_limma <- model.matrix(~group_info)

fit <- lmFit(gsva_gobp,design_limma)
fit <- eBayes(fit)
limma_results <- topTable(fit, n = Inf)

# Heatmap 
keep <- rownames(limma_results[limma_results$P.Val < 0.05 & abs(limma_results$logFC)>0.4, ])
heat_dat <- gsva_gobp[keep,]

pheatmap(heat_dat)

































