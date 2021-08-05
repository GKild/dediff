####other fetal lung reference#######
library(Seurat)
library(Matrix)
fetal_lung_mtx=read.table('fetal_lung_elife/Combined_raw.txt', sep='\t', header=T)
fetal_lung_mtx$X=make.unique(fetal_lung_mtx$X)
rownames(fetal_lung_mtx)=fetal_lung_mtx$X
fetal_lung_mtx=fetal_lung_mtx[,2:ncol(fetal_lung_mtx)]

fetal_lung_mtx=as.matrix(fetal_lung_mtx)
fetal_lung_mtx=Matrix(fetal_lung_mtx, sparse = T)

fetal_lung_annot=read.csv('fetal_lung_elife/Cell_Type_annotation.csv')

fetal_lung_srat=CreateSeuratObject(fetal_lung_mtx)
fetal_lung_srat@meta.data$orig_clust=fetal_lung_annot$Cluster
fetal_lung_srat@meta.data$annot=fetal_lung_annot$Type

fetal_lung_srat_processed=process_seurat(fetal_lung_srat)