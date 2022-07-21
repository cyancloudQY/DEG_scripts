library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(enrichplot)

### Obtain the gene name and log2FC
gsea_list <- read.table('./DESeq_results.txt')

gsea_file <- rownames_to_column(gsea_list, 'gene_name')
gsea_file <- dplyr::select(gsea_file, gene_name, log2FoldChange)
gsea_file <- gsea_file[order(gsea_file$log2FoldChange, decreasing = TRUE),]

bitr <- bitr(gsea_file$gene_name,
             fromType = "ENSEMBL",
             toType = "ENTREZID",
             OrgDb = 'org.Hs.eg.db')

gsea_file <- left_join(gsea_file, bitr, by = c('gene_name' = 'ENSEMBL'))

gsea_file <- gsea_file[complete.cases(gsea_file), ]


### Define gene_list ###
geneList <- gsea_file$log2FoldChange 
names(geneList) <- gsea_file$ENTREZID

# GO
gsea_go <- gseGO(geneList = geneList,
                 OrgDb = org.Hs.eg.db,
                 ont = 'ALL',
                 verbose = FALSE,
                 seed = FALSE,
                 by = 'fgsea')


View(gsea_go@result)

gseaplot2(gsea_go, 1, color = "firebrick", rel_heights=c(1, .2, .6))
features <- c('GO:0090049', 'GO:0051276', 'GO:0030334', 'GO:0000793')
gseaplot2(gsea_go, features, color = "firebrick", rel_heights=c(1, .2, .6))


# KEGG
gsea_kegg <- gseKEGG(geneList = geneList,
                 organism = 'hsa',
                 verbose = FALSE,
                 seed = FALSE,
                 by = 'fgsea')


View(gsea_kegg@result)
gseaplot2(gsea_kegg, 1, color = "firebrick", rel_heights=c(1, .2, .6))
gseaplot2(gsea_kegg, 2, color = "firebrick", rel_heights=c(1, .2, .6))



##################### Self-define ######################
gmt <- read.gmt('./msigdb.v7.5.1.entrez.gmt')

## Select desired pathways
gmt_mir <- gmt[grep('^MIR',gmt$term),]
mir <- GSEA(geneList,TERM2GENE = gmt_mir)
View(mir@result)
gseaplot2(mir, 1, color = "firebrick", rel_heights=c(1, .2, .6))








