library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
source('logisticRegression.R')


colon_canc=read.csv("/home/jovyan/Dediff/sc_crc/GSE81861_CRC_tumor_all_cells_COUNT.csv")
rownames(colon_canc)=colon_canc$X
colon_canc=colon_canc[,2:ncol(colon_canc)]
colon_canc_mtx=round(as.matrix(colon_canc), digits = 0)
colon_matrix_sparse=Matrix(colon_canc_mtx, sparse = T)
rownames(colon_matrix_sparse)=make.unique(sapply(strsplit(rownames(colon_canc), "_"), '[',2))
colon_canc_annot=sapply(strsplit(colnames(colon_matrix_sparse), "__"), '[',2)

comm_genes=intersect(rownames(colon_matrix_sparse), rownames(srat_comb_nocc))

colon_train_canc=trainModel(srat_comb_nocc@assays$RNA@counts[comm_genes,], srat_comb_nocc@meta.data$Cluster, workers=NULL)
colon_pred_canc=predictSimilarity(colon_train_canc, colon_matrix_sparse[comm_genes,], colon_canc_annot)
similarityHeatmap(colon_pred_canc)

saveRDS(col_canc_pred, "/home/jovyan/Dediff/col_canc_pred.rds")

x=sapply(unique(colon_seurat@meta.data$exp_annot), function(x){
  colnames(colon_seurat)[which(colon_seurat@meta.data$exp_annot==x)][1:1000]
})

similarityHeatmap(confusion_pred, row_order=unique(colon_seurat@meta.data$cell_group_age), column_order=unique(colon_seurat@meta.data$cell_group_age))
unlist(x, use.names = F)
x=as.vector(x)
x<-x[!is.na(x)]
x=unique(x)

colon_tr_expanded=trainModel(colon_seurat@assays$RNA@counts[comm_genes,], colon_seurat@meta.data$exp_annot, workers=NULL)
col_exp_canc_pred=predictSimilarity(colon_tr_expanded, colon_matrix_sparse[comm_genes,], colon_canc_annot)

similarityHeatmap(col_exp_canc_pred)


rownames(colon_seurat)[1:100]
obs$gene_ids


similarityHeatmap(confusion_pred,row_order=ann_ord, column_order=ann_ord)
