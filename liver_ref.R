library(Matrix)
library(Seurat)
library(ComplexHeatmap)
library(ggplot2)
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

#Fetal liver - from Muzz's lab
l_f_mtx=readMM('fetal_liver/matrix.mtx')
l_f_genes=read.csv('fetal_liver/obs.csv')
l_f_meta=read.csv('fetal_liver/var.csv')
rownames(l_f_meta)=l_f_meta$index

rownames(l_f_mtx)=l_f_genes$index
colnames(l_f_mtx)=l_f_meta$index

fetal_liver_srat=CreateSeuratObject(l_f_mtx, meta.data = l_f_meta)

fetal_liver_srat_processed=process_seurat(fetal_liver_srat)

DimPlot(fetal_liver_srat_processed, group.by = "cell.labels")

adult_liver=read.csv('GSE115469_adult_liver/FilteredCounts.csv', row.names = 1)
adult_liver=as.matrix(adult_liver)
adult_liver=Matrix(adult_liver, sparse=T)
adult_liver_clusts=read.table('GSE115469_adult_liver/Cell_clusterID_cycle.txt', sep = "\t", header = T)
rownames(adult_liver_clusts)=adult_liver_clusts$CellName

adult_liver_srat=CreateSeuratObject(adult_liver, meta.data = adult_liver_clusts)

adult_liver_srat_processed=process_seurat(adult_liver_srat)

DimPlot(adult_liver_srat_processed, group.by = "Cluster.")

clusts = 1:20
annot= c("Hep_1","ab_Tcells","Hep_2","Inflammatory_macs", "Hep_3",
         "Hep_4", "Plasma_cells", "NK-like_cells", "gd_Tcells_1", "Non-inflammatory_macs", 
         "Periportal_LSECs", "Central_venous_LSECs", "Portal_endothelium", "Hep_5", "Hep_6",
         "Mature_Bcells", "Cholangiocytes", "gd_Tcells_2", "Erythroid_cells", "Hepatic_stellate_cells")

adult_liver_srat@meta.data$Cluster.[adult_liver_srat@meta.data$Cluster. %in% clusts] <- annot[match(adult_liver_srat@meta.data$Cluster., clusts, nomatch = 0)]

adult_liver_srat_processed$annot=adult_liver_srat@meta.data$Cluster.

VlnPlot(adult_liver_srat_processed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

adult_liver_srat_processed[["percent.mt"]] = PercentageFeatureSet(adult_liver_srat_processed, pattern = "^MT-")

fetal_liver_srat_processed@meta.data$annot=fetal_liver_srat_processed@meta.data$cell.labels
fetal_liver_srat_processed@meta.data$annot=paste0("Fetal_",fetal_liver_srat_processed@meta.data$annot)

adult_liver_srat_processed@meta.data$annot=paste0("Adult_",adult_liver_srat_processed@meta.data$annot)

merged_liver_srat=merge(fetal_liver_srat_processed, adult_liver_srat_processed)

merged_liver_srat_processed=process_seurat(merged_liver_srat)

saveRDS(merged_liver_srat_processed, 'merged_liver/merged_liver_processed.rds')
merged_liver_srat_processed=readRDS('merged_liver/merged_liver_processed.rds')

merged_liver_srat_processed@meta.data$annot[grep("Adult_Hep_", merged_liver_srat_processed@meta.data$annot)]="Adult_Hepatoctes"

merged_liver_srat_processed_nocc=subset(merged_liver_srat_processed, Phase=="G1")
DimPlot(merged_liver_srat_processed, group.by = "Phase", label=T) +theme(legend.position ="none")
#################################
#Matts cor
dat = merged_liver_srat_processed@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),merged_liver_srat_processed@meta.data$annot),
                           function(e) rowMeans(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)
#################################
#Gerda's cor
Idents(merged_liver_srat_processed)=merged_liver_srat_processed@meta.data$annot
avg_exp=AverageExpression(merged_liver_srat_processed, slot="scale.data")

cordat=cor(avg_exp$RNA)
Heatmap(cordat)
#################################

ensembl_ids=read.csv('onecs/obs.csv')


inter=intersect(rownames(merged_liver_srat_processed),ensembl_ids$index)

ensembl_inter=ensembl_ids[inter,]
#function to make a Cell Signal analysis-compatible files from a seurat object that doesn't have ENSEMBL IDs already. 
#annotation needs to be under annot column in the metadata
write_cellsig_files=function(srat, ensembl_id_file, out_location){
  rownames(ensembl_id_file)=ensembl_id_file$index
  inter=intersect(rownames(srat),ensembl_id_file$index)
  ensembl_inter=ensembl_id_file[inter,]
  mat=srat@assays$RNA@counts[inter,]
  rownames(mat)=ensembl_inter$gene_ids
  colnames(mat)=paste0(srat@meta.data$annot,":",colnames(mat))
  writeMM(mat, paste0(out_location, "liver.mtx"))
  write.table(rownames(mat), paste0(out_location, "liver_rowNames.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)
  write.table(colnames(mat), paste0(out_location, "liver_columnNames.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)
  
}

write_cellsig_files(merged_liver_srat_processed, ensembl_ids, '/home/jovyan/Dediff/cellsig_liver/liver_ref_cycling/')
write_cellsig_files(merged_liver_srat_processed_nocc, ensembl_ids, '/home/jovyan/Dediff/cellsig_liver/liver_ref_noncycling/')



