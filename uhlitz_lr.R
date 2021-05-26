library(Matrix)
library(Seurat)
source('logisticRegression.R')
library(loomR)
workers=NULL
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

ref_comb=readRDS('/lustre/scratch117/casm/team274/gk14/crc/rasa_and_smiley_combined_ref.rds')
uhlitz_meta=read.delim('/home/jovyan/Dediff/crc_uhlitz/GSE166555_meta_data.tsv', sep = "\t", header = T)
lfile = connect(filename = "/home/jovyan/Dediff/crc_uhlitz/sccrc_rawcounts.loom", mode = "r+")

gene_names=lfile[["row_attrs/Gene"]][]
cell_names=lfile[["col_attrs/CellID"]][]
uhlitz_mtx=Matrix(t(full.matrix), sparse = T)

rownames(uhlitz_mtx)=gene_names
colnames(uhlitz_mtx)=cell_names

rownames(uhlitz_meta)=uhlitz_meta$cell


uhlitz_srat=CreateSeuratObject(uhlitz_mtx, meta.data = uhlitz_meta)

saveRDS(uhlitz_srat, '/lustre/scratch117/casm/team274/gk14/Dediff/uhlitz_srat.rds')

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
ref_comb=process_seurat(ref_comb)

ref_comb_merge_ct=ref_comb
ref_comb@meta.data$merged_idents=ref_comb@meta.data$celltype
ref_comb@meta.data$merged_idents[which(ref_comb@meta.data$celltype%in%c("Stem", "TA 1", "TA 2"))]="Stem/TA"
ref_comb@meta.data$merged_idents[which(ref_comb@meta.data$celltype%in%c("Immature Enterocytes 1", "Immature Enterocytes 2"))]="Immature Enterocytes"

DimPlot(ref_comb, cells.highlight = rownames(ref_comb@meta.data)[which(ref_comb@meta.data$celltype%in%c("Stem", "TA 1", "TA 2"))])
uhlitz_srat=process_seurat(uhlitz_srat)

DimPlot(uhlitz_srat,cells.highlight = colnames(uhlitz_srat)[which(uhlitz_srat@meta.data$cell_type_epi_custom%in%c("TC1","TC2", "TC3","TC4"))])
DimPlot(uhlitz_srat, group.by = "sct_cell_type", label=T)
uhlitz_epi=subset(uhlitz_srat, main_cell_type=="Epithelial")

uhlitz_epi_tumor=subset(uhlitz_epi, sample_origin=="Tumor")



comm_genes=intersect(rownames(uhlitz_epi), rownames(ref_comb))

ref_train=trainModel(ref_comb@assays$RNA@counts[comm_genes,], ref_comb@meta.data$celltype, workers=NULL)
ref_train_merged=trainModel(ref_comb@assays$RNA@counts[comm_genes,], ref_comb@meta.data$merged_idents, workers=NULL)

col_ord=c("Stem","fetal_crypt","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
          "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
          "Immature Goblet","Secretory TA","fetal_Goblet cell",
          "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
          "fetal_Enterocyte","NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
          "Macrophages","CD4+ T-cells","Tregs",
          "Follicular","Inflammatory Monocytes","ILCs","Plasma")
col_ord_merged=c("Stem/TA","fetal_crypt","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","M cells",
                 "Enterocytes","Best4+ Enterocytes","Immature Enterocytes",
                 "Immature Goblet","Secretory TA","fetal_Goblet cell",
                 "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
                 "fetal_Enterocyte","NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
                 "Macrophages","CD4+ T-cells","Tregs",
                 "Follicular","Inflammatory Monocytes","ILCs","Plasma")

ps_uhlitz_tum=predictSimilarity(ref_train, uhlitz_epi_tumor@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor@meta.data$cell_type_epi_custom)
ps_uhlitz_tum_merged=predictSimilarity(ref_train_merged, uhlitz_epi_tumor@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor@meta.data$cell_type_epi_custom)
ps_uhlitz=predictSimilarity(ref_train, uhlitz_srat@assays$RNA@counts[comm_genes,], uhlitz_srat@meta.data$cell_type_epi)

similarityHeatmap(ps_uhlitz, column_order=col_ord)
similarityHeatmap(ps_uhlitz_tum, use_raster=F, column_order=col_ord, row_order=row_ord)
similarityHeatmap(ps_uhlitz_tum_merged, use_raster=F, column_order=col_ord_merged)


ps_uhlitz_tum=predictSimilarity(ref_train, uhlitz_epi_tumor@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor@meta.data$cell_type_epi_custom)

similarityHeatmap(ps_uhlitz_tum, use_raster=F)
ref_markers=getMarkers(ref_train)

ref_markers[ref_markers$class=="fetal_crypt",]


ref_stem_only=subset(ref_comb, merged_idents%in%c("fetal_crypt", "Stem/TA"))
ref_stem_train=trainModel(ref_stem_only@assays$RNA@counts[comm_genes,], ref_stem_only@meta.data$merged_idents, workers=NULL)

ref_stem_ps=predictSimilarity(ref_stem_train, uhlitz_epi_tumor_mt_filt@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor_mt_filt@meta.data$cell_type_epi_custom)

similarityHeatmap(ref_stem_ps, row_order=row_ord)

unique(uhlitz_epi_tumor@meta.data$cell_type_epi_custom)

obs=read.csv("onecs/obs.csv")

intersect(obs$index, rownames(uhlitz_epi_tumor))



obs_mut=obs[rownames(uhlitz_epi_tumor),]

identical(rownames(uhlitz_epi_tumor), obs_mut$index)


tum_mtx=uhlitz_epi_tumor@assays$RNA@counts
rownames(tum_mtx)=obs_mut$gene_ids

uhlitz_epi_tumor@meta.data$cell_type_epi_custom=str_replace(uhlitz_epi_tumor@meta.data$cell_type_epi_custom, "/", "_")
for (x in unique(uhlitz_epi_tumor@meta.data$cell_type_epi_custom)) {
  rs=rowSums(tum_mtx[,rownames(uhlitz_epi_tumor@meta.data)[which(uhlitz_epi_tumor@meta.data$cell_type_epi_custom==x)]])
  df=data.frame(geneName=rownames(tum_mtx))
  df$geneLengths=1
  df$bla=rs
  colnames(df)=c("geneName", "geneLengths", x)
  write.table(df, paste0('/home/jovyan/Dediff/cellsig_gut/pseudobulk_in/', x, ".tsv"), sep = "\t", quote=F, row.names = F)
}


pseudobulk_paths=paste0('/home/jovyan/Dediff/cellsig_gut/pseudobulk_in/', unique(uhlitz_epi_tumor@meta.data$cell_type_epi_custom), ".tsv")
write.table(pseudobulk_paths, '/home/jovyan/Dediff/cellsig_gut/pseudobulk_paths.txt', sep = "\t", quote=F, row.names = F, col.names = F)                       


cellsig='cellsig_gut/pseudobulk_out/OutRun_fitExposures.tsv'
fit=normaliseExposures(cellsig)
plotExposures(fit)

uhlitz_epi_tumor=process_seurat(uhlitz_epi_tumor)
uhlitz_epi_tumor_nocc=subset(uhlitz_epi_tumor, Phase=="G1")

ps_epi_tumor_nocc=predictSimilarity(ref_train, uhlitz_epi_tumor_nocc@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor_nocc@meta.data$cell_type_epi_custom )
similarityHeatmap(ps_epi_tumor_nocc)
uhltiz_rerun=predictSimilarity(ref_train, uhlitz_epi_tumor@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor@meta.data$RNA_snn_res.1)

similarityHeatmap(uhltiz_rerun, column_order=col_ord)


dat.anchors = FindTransferAnchors(reference = ref_comb, query = uhlitz_epi_tumor, 
                                  dims = 1:30)
predictions = TransferData(anchorset = dat.anchors, refdata = ref_comb@meta.data$celltype, 
                           dims = 1:30)
uhlitz_epi_tumor = AddMetaData(uhlitz_epi_tumor, metadata = predictions)

g1=DimPlot(uhlitz_epi_tumor, group.by = "predicted.id", label = T)
g2=DimPlot(uhlitz_epi_tumor, group.by = "cell_type_epi_custom", label = T)

FeaturePlot(uhlitz_epi_tumor, "prediction.score.Stem")
g1=DimPlot(uhlitz_epi_tumor, group.by = "predicted.id", label = T)
g2=DimPlot(uhlitz_epi_tumor, group.by = "cell_type_epi_custom", label = T)