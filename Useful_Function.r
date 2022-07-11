library(pheatmap)
library(ggplot2)


#### Volcano Plot ####
# 1. voldata need to have columns contains information about the p-value and log2FoldChange value.
# 2. Default input voldata could be directly from the DESeq. If not, need to change the column names.
# 3. p,fc are the threshold of p-value and fold change value.
vol_plot <- function(voldata,p,fc){
  voldata$color <- ifelse(voldata$pvalue <p & abs(voldata$log2FoldChange)> fc,
                          ifelse(voldata$log2FoldChange > fc,'red','blue'),'gray') 
  
  color <- c(red = "red",gray = "gray",blue = "blue")
  volplot <- ggplot(voldata, aes(log2FoldChange, -log10(pvalue), col = color)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = color) +
    labs(x="log2 (fold change)",y="-log10 (p-value)") +
    geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
    geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
    theme(legend.position = "none",
          panel.grid=element_blank(),
          axis.title = element_text(size = 16),
         axis.text = element_text(size = 14))
    return (volplot)
    }

#### Heatmap ####
# 1. exp should be the expression matrix that is already normalized (cpm).
# 2. exp should contain the column names and row names. 

heat_map <- function(exp){
  pheatmap(exp, legend = FALSE, cluster_cols = FALSE, scale = 'row',clustering_method = "average", 
           show_rownames = FALSE, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  
}

#### PCA ####
# 1. exp should be normalized (cpm) expression matrix.
# 2. group_info should be factor list. 
# eg: group_info <- factor(c(rep("R",10), rep('S',9)))
library(FactoMineR)
> PCA_exp <- function(exp, group_info){
+     pca_data <- t(exp)
+     gene.pca <- PCA(pca_data, ncp = 2, scale.unit = TRUE, graph = FALSE)
+     pca_sample <- data.frame(gene.pca$ind$coord)
+     group_df <- as.data.frame(group_info)
+     pca_sample <- cbind(pca_sample, group_df)
+     write.table(pca_sample, './pca_info.txt', quote = FALSE)
+     
+     PCA_g <- plot(gene.pca)
+     return(PCA_g)
+ }




