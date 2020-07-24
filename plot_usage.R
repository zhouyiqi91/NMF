library(Seurat)
library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds file")
argv <- add_argument(argv,"--usage_file", help="usage_file")
argv <- add_argument(argv,"--sample", help="sample name")
argv <- add_argument(argv,"--k", help="k selection")
argv <- add_argument(argv,"--reduction", help="tsne or umap")
argv <- parse_args(argv)


rds = readRDS(argv$rds)
usage.norm.file = argv$usage_file
sample.name = argv$sample
k = argv$k
reduction = argv$reduction

usage.norm = read_tsv(usage.norm.file)
usage.norm = tibble::column_to_rownames(usage.norm,"index")
usage.norm = as.data.frame(usage.norm)
meta = rds@meta.data
new.meta = transform(merge(usage.norm,meta,by=0,all=T),row.names=Row.names, Row.names=NULL)
rds@meta.data = new.meta
cols = paste("Usage",c(1:k),sep="_")
pdf.name = paste(sample.name,"usage_plot.pdf",sep="_")
pdf(pdf.name,height=14,width=10)
FeaturePlot(rds,features.plot = cols,no.legend = F,nCol=2,reduction.use = reduction)
dev.off()