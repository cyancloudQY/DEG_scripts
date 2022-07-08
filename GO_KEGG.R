library(optparse)
library(tidyverse)
library(clusterProfiler)

option_list <- list(
    make_option(c("-i","--input")),
	make_option(c("-w","--work"))
)

opt <- parse_args(OptionParser(option_list=option_list))

input_file <- opt$input
work_dir <- opt$work


gene_list <- read.table(input_file, header = FALSE)
setwd(work_dir)

gene_list <- gene_list$V1



######## GO ###########

ego_all <- enrichGO(gene = gene_list, 
                    OrgDb = 'org.Hs.eg.db', 
                    ont = "ALL", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    keyType = 'ENSEMBL',
                    readable = TRUE)


pdf('./go_barplot.pdf')
barplot(ego_all,showCategory=10)
dev.off()


pdf('./go_dotplot.pdf')
dotplot(ego_all,showCategory=10)
dev.off()


results_GO <- ego_all@result
write.table(results_GO, 'results_GO.txt',quote = FALSE)




####### KEGG ##########
entrez_list <- bitr(gene_list, fromType = 'ENSEMBL',toType = 'ENTREZID', OrgDb="org.Hs.eg.db")




kegg_all <- enrichKEGG(gene = entrez_list$ENTREZID,
                       organism = "hsa",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)
					   
pdf('./kegg_barplot.pdf')
barplot(kegg_all,showCategory=10)
dev.off()


pdf('./kegg_dotplot.pdf')
dotplot(kegg_all,showCategory=10)
dev.off()					 

results_KEGG <- kegg_all@result
write.table(results_KEGG, 'results_KEGG.txt',quote = FALSE)









