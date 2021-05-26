library(Matrix)
library(Seurat)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
ref_mtx=readMM('cellsig_gut/sc_combined/combined.mtx')
ref_cols=read.table('cellsig_gut/sc_combined/combined_columnNames.tsv', header = F, sep = "\t")
ref_rows=read.table('cellsig_gut/sc_combined/combined_rowNames.tsv', header = F, sep = "\t")

colnames(ref_mtx)=ref_cols$V1
rownames(ref_mtx)=ref_rows$V1

annot=sapply(strsplit(ref_cols$V1, ":"), "[", 1)

ref_srat=CreateSeuratObject(ref_mtx)

Idents(ref_srat)=annot

process_10x =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  #srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat, features=rownames(srat))
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

ref_srat=process_10x(ref_srat)

Idents(ref_srat)=annot

avg_exp=AverageExpression(ref_srat, slot="scale.data")

cordat=cor(avg_exp$RNA)

Heatmap(cordat)
dim(avg_exp$RNA)


dim(ref_srat)
