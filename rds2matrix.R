library(Seurat)
argv = commandArgs()
rds = readRDS(argv[2])
raw.data = rds@raw.data

dir.create("matrix")
setwd("matrix")
setwd("matrix")
writeMM(mtx, "matrix.mtx")
genes = as.data.frame(rownames(mtx))
genes[,2] = "fake_ID"
genes = genes[,c(2,1)]
cell_barcodes = as.data.frame(colnames(mtx))
write.table(genes,"genes.tsv",row.names = F,col.names=F,sep="\t")
write.table(cell_barcodes,"barcodes.tsv",row.names = F,col.names=F,sep="\t")