## DESeq + Limma + EdgeR ##
library(VennDiagram)
library(DESeq2)
library(tidyverse)
library(stringr)
library(pheatmap)
library(edgeR)
library(limma)
library(optparse)
library(FactoMineR)

# getopt
option_list <- list(
    make_option(c("-i","--input")),
	make_option(c("-g","--group")),
	make_option(c("-w","--work")),
	make_option(c("-f","--log2fc"), default = 2),
	make_option(c("-p","--pvalue"), default = 0.05)	
)
opt <- parse_args(OptionParser(option_list=option_list))

exp_matrix <- opt$input
work_dir <- opt$work
group_file <- opt$group
fold_change <- as.numeric(opt$log2fc)
p_value <- as.numeric(opt$pvalue)



### Read in the clean count matrix ### 
# 1. Columns Should be arranged based on the order of the biological duplication.
# 2. Should be tab seperated.
# 3. Should be raw counts without any normalization.
# 4. Row names should be the gene names.

raw_counts <- read.table(exp_matrix, header = TRUE, row.names = 1)

### Read in the group info ###
# 1. Table with the information of the number of samples in each group.
# 2. Table should be created based on the order of the colnames of count matrix.
# R 10
# S 9
g_file <- read.table(group_file)
group_info <- factor(c(rep(g_file$V1[1], g_file$V2[1]),rep(g_file$V1[2], g_file$V2[2])))

print(group_info)

setwd(work_dir)
### DESeq Analysis ###
DESeq_object <- DESeqDataSetFromMatrix(raw_counts, colData = DataFrame(group_info), design= ~group_info)
DESeq_object <- DESeq(DESeq_object)
DESeq_results <- as.data.frame(results(DESeq_object))

DESeq_results$type <- ifelse((DESeq_results$log2FoldChange > fold_change) & (DESeq_results$pvalue < p_value), 'UP',
							 ifelse((DESeq_results$log2FoldChange < -fold_change) & (DESeq_results$pvalue < p_value), 'DOWN','NOT'))

DESeq_results <- DESeq_results[order(DESeq_results[,5]),]

DESeq_up <- filter(DESeq_results, DESeq_results$log2FoldChange > fold_change & DESeq_results$pvalue < p_value)
DESeq_down <- filter(DESeq_results, DESeq_results$log2FoldChange < -fold_change & DESeq_results$pvalue < p_value)
DESeq_diff <- filter(DESeq_results, (DESeq_results$log2FoldChange < -fold_change | DESeq_results$log2FoldChange > fold_change )
  					& DESeq_results$pvalue < p_value)


DESeq_up <- rownames(DESeq_up)
DESeq_down <- rownames(DESeq_down)


write.table(DESeq_diff, './DEGSeq_diff.txt', quote = FALSE, sep = '\t')				
write.table(DESeq_results, './DEGSeq_results.txt', quote = FALSE, sep = '\t')			
					
					
###	EdgeR Analysis ###
### https://zhuanlan.zhihu.com/p/250974725 
### Option to Integrate more complex model (design with multiple factors)

edgeR_object <- DGEList(counts = raw_counts, group = group_info)			
keep <- filterByExpr(edgeR_object, group = group_info)	
edgeR_object <- edgeR_object[keep,,keep.lib.sizes=FALSE]
edgeR_object <- calcNormFactors(edgeR_object, method = 'TMM')
edgeR_object <- estimateDisp(edgeR_object)
edgeR_results <- exactTest(edgeR_object)
edgeR_results <- edgeR_results$table

edgeR_results$type <- ifelse((edgeR_results$logFC > fold_change) & (edgeR_results$PValue < p_value), 'UP',
							 ifelse((edgeR_results$logFC < -fold_change) & (edgeR_results$PValue < p_value), 'DOWN','NOT')) 
edgeR_results <- edgeR_results[order(edgeR_results[,3]),] 
 
 
edgeR_up <- filter(edgeR_results, edgeR_results$logFC > fold_change & edgeR_results$PValue < p_value)
edgeR_down <- filter(edgeR_results, edgeR_results$logFC < -fold_change & edgeR_results$PValue < p_value)
edgeR_diff <- filter(edgeR_results, (edgeR_results$logFC < -fold_change | edgeR_results$logFC > fold_change )
  					& edgeR_results$PValue < p_value) 
 
 
edgeR_up <- rownames(edgeR_up)
edgeR_down <- rownames(edgeR_down) 


write.table(edgeR_diff, './edgeR_diff.txt', quote = FALSE, sep = '\t')
write.table(edgeR_results, './edgeR_results.txt', quote = FALSE, sep = '\t')
				
					
					
### limma analysis ###
### voom transformation is for samples with much difference between library sizes

design_limma <- model.matrix(~group_info)

limma_object <- DGEList(counts = raw_counts)
limma_object <- calcNormFactors(limma_object)
limma_object <- voom(limma_object,design_limma)

fit <- lmFit(limma_object,design_limma)
fit <- eBayes(fit)
limma_results <- topTable(fit, n = Inf)

limma_results$type <- ifelse((limma_results$logFC > fold_change) & (limma_results$P.Value < p_value),'UP',
							ifelse((limma_results$logFC < -fold_change) & (limma_results$P.Value < p_value),'DOWN', 'NOT'))

limma_results <- limma_results[order(limma_results[,4]),]



limma_up <- filter(limma_results, limma_results$logFC > fold_change & limma_results$P.Value < p_value)
limma_down <- filter(limma_results, limma_results$logFC < -fold_change & limma_results$P.Value < p_value)
limma_diff <- filter(limma_results, (limma_results$logFC < -fold_change | limma_results$logFC > fold_change)
					& limma_results$P.Value < p_value)


limma_up <- rownames(limma_up)
limma_down <- rownames(limma_down)


write.table(limma_results, './limma_results.txt', quote = FALSE, sep = '\t')
write.table(limma_diff, './limma_diff.txt', quote = FALSE, sep = '\t')



### CPM normalization 
cpm_dat <- log2(cpm(raw_counts)+1)
write.table(cpm_dat, './CPM_exp.txt', sep = '\t', quote = FALSE)


### Intersection Analysis
venn.plot <- venn.diagram(
  x = list(
    DESeq = DESeq_up,
    edgeR = edgeR_up,
    limma = limma_up
  ),
  filename = "Venn_up.png", imagetype = 'png',
  col = "transparent",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
)

venn.plot <- venn.diagram(
  x = list(
    DESeq = DESeq_down,
    edgeR = edgeR_down,
    limma = limma_down
  ),
  filename = "Venn_down.png", imagetype = 'png',
  col = "transparent",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
)


intersection_diff_genes <- Reduce(intersect,list(rownames(DESeq_diff), rownames(limma_diff), rownames(edgeR_diff)))
write.table(intersection_diff_genes, 'intersection_gene_list.txt', sep = '\t', quote = FALSE, row.names = FALSE)



















