library(Seurat)
library(tidyverse)
argv = commandArgs()
print(argv[2])
rds = readRDS(argv[2])
usage.norm.file = argv[3]
sample.name = argv[4]
k = argv[5]

usage.norm = read_tsv(usage.norm.file)
usage.norm = tibble::column_to_rownames(usage.norm,"index")
usage.norm = as.data.frame(usage.norm)
meta = rds@meta.data
new.meta = transform(merge(usage.norm,meta,by=0,all=T),row.names=Row.names, Row.names=NULL)
rds@meta.data = new.meta
cols = paste("Usage",c(1:k),sep="_")
pdf.name = paste(sample.name,"usage_plot.pdf",sep="_")
pdf(pdf.name,height=14,width=10)
FeaturePlot(rds,features.plot = cols,no.legend = F)
dev.off()