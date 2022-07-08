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

print(gene_list)






