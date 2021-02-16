library(Matrix)
library(Seurat)
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

epi_seurat_nocc=subset(epi_seurat, cells=rownames(epi_seurat@meta.data)[epi_seurat@meta.data$Phase=="G1"])

source('logisticRegression.R')

li_seurat=readRDS('large_intestine_rasa.rds')

li_seurat_adult=subset(li_seurat, Age%in%c("A34", "A32", "A38", "A26", "A39"))
li_seurat_fetal=subset(li_seurat, Age%in%c("F78", "F66", "F67", "F72", "F73"))
lif_nocc_seurat=subset(li_seurat_fetal, Phase=="G1")

li_seurat_adult_nocc=subset(li_seurat_adult, Phase=="G1")

comm_genes=intersect(rownames(li_seurat), rownames(epi_seurat))

epi_cc_m=trainModel(epi_seurat@assays$RNA@counts[comm_genes,], epi_seurat@meta.data$Cluster, workers=NULL)

epi_cc_ps=predictSimilarity(epi_cc_m, li_seurat_adult@assays$RNA@counts[comm_genes,], li_seurat_adult@meta.data$exp_annot)

similarityHeatmap(epi_cc_ps)

epi_nocc_m=trainModel(epi_seurat_nocc@assays$RNA@counts[comm_genes,], epi_seurat_nocc@meta.data$Cluster, workers=NULL)
epi_nocc_ps=predictSimilarity(epi_nocc_m, li_seurat_adult_nocc@assays$RNA@counts[comm_genes,], li_seurat_adult_nocc@meta.data$exp_annot)
similarityHeatmap(epi_nocc_ps)

obs_shared=obs[which(obs$index%in%rownames(epi_seurat_nocc)),]

comm_ensg=intersect(obs_shared$gene_ids, rownames(bulk_mat))

epi1=epi_seurat_nocc@assays$RNA@counts[obs_shared$index,]
rownames(epi1)=obs_shared$gene_ids                      

crc_bulk_mat=bulk_mat[comm_ensg, mDat_col$UniqueSampleID]
rownames(len)=make.unique(as.character(len$geneName))
len=len[which(rownames(len)%in%comm_ensg),]
gene_lengths=len$geneLengths

tr_for_bulk=trainModel(epi1[comm_ensg,], epi_seurat_nocc$Cluster)

predict_bulk=predictSimilarity(tr_for_bulk, crc_bulk_mat, lengths = gene_lengths)

similarityHeatmap(predict_bulk)


