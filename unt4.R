library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
source('logisticRegression.R')

#Regev data
epi_mtx=readMM('regev_adult_gut/gene_sorted-Epi.matrix.mtx')
epi_barcodes=read.table('regev_adult_gut/Epi.barcodes2.tsv')
epi_genes=read.table('regev_adult_gut/Epi.genes.tsv')
meta=read.table('regev_adult_gut/all.meta2.txt', sep = "\t", header = T)
rownames(meta)=meta$NAME
rownames(epi_mtx)=epi_genes$V1
colnames(epi_mtx)=epi_barcodes$V1
epi_meta=meta[colnames(epi_mtx), ]
healthy_epi_meta=epi_meta[epi_meta$Health=="Healthy",]
healthy_epi_mtx=epi_mtx[,rownames(healthy_epi_meta)]
epi_seurat=CreateSeuratObject(healthy_epi_mtx, meta.data = healthy_epi_meta)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
process_seurat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

epi_seurat=process_seurat(epi_seurat)
#regev data - healthy patients, remove cycling cells
epi_seurat_nocc=subset(epi_seurat, cells=rownames(epi_seurat@meta.data)[epi_seurat@meta.data$Phase=="G1"])

epi_seurat_nocc=subset(epi_seurat_nocc, cells=rownames(epi_seurat@meta.data)[epi_seurat@meta.data$Cluster!="Cycling TA"])



#Rasa's data - fetal. Only epithelial cells. Cycling cells removed
li_seurat=readRDS('large_intestine_rasa.rds')

li_seurat_adult=subset(li_seurat, Age%in%c("A34", "A32", "A38", "A26", "A39"))
li_seurat_fetal=subset(li_seurat, Age%in%c("F78", "F66", "F67", "F72", "F73"))
lif_nocc_seurat=subset(li_seurat_fetal, Phase=="G1")

lif_nocc_seurat=subset(li_seurat_fetal, exp_annot%in%c("fetal_Goblet cell","fetal_Enterocyte",
                                                       "fetal_BEST4 enterocyte","fetal_crypt",
                                                       "fetal_Enteroendocrine","fetal_Tuft cell", 
                                                       "fetal_TA"))




