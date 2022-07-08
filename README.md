# DEG_scripts
Scripts that can be used in differentially expressed genes analysis.

# LiDeEd
**Runs limma, DESeq, edgeR all together**
## Input
1. **-i** input gene expression matrix

a. Should be tab delimited file, with sample names as column names and gene names as row names.

b. Should be raw counts without any normalization.

c. Columns should be arranged in a manner that the same biological duplication is in the same side. (C1 C2 C3 T1 T2 T3)

2. **-g** group information

a. Table with information of the number of samples in each group.

b. Should be ordered using the same manner as the column names of the expression matrix. 

c. An example.


| Group | Number of Samples |
| :---: | :-----------: |
| R | 10|
| S | 9 |

3. **-w** output directory

4. **-f** threshold of log2fold change

5. **-p** threshold of p-value.

## Output
1. CPM_exp : Expression matrix after CPM normalization. 
2. DESeq_diff, DESeq_results, edgeR_diff, edgeR_results, limma_diff, limma_results : results file contains the full output of three DEG software, diff file contains the genes that only pass the given threshold.
3. intersection_gene_list: Genes that are identified by all three software.
4. Venn_up, Venn_down: Venn diagram of the three softwares.  


