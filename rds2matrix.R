library(Seurat)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="rds file")
argv <- parse_args(argv)
print ("loading rds")
rds = readRDS(argv$rds)
mtx = rds@raw.data

dir.create("matrix")
setwd("matrix")
writeMM(mtx, "matrix.mtx")
genes = as.data.frame(rownames(mtx))
genes[,2] = "fake_ID"
genes = genes[,c(2,1)]
cell_barcodes = as.data.frame(colnames(mtx))
write.table(genes,"genes.tsv",row.names = F,col.names=F,sep="\t")
write.table(cell_barcodes,"barcodes.tsv",row.names = F,col.names=F,sep="\t")