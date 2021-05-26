uhlitz_srat=readRDS('/lustre/scratch117/casm/team274/gk14/Dediff/uhlitz_srat.rds')
source('logisticRegression.R')

uhlitz_epi=subset(uhlitz_srat, main_cell_type=="Epithelial")

uhlitz_epi_tumor=subset(uhlitz_epi, sample_origin=="Tumor")

process_seurat2 =function(srat){
  srat=subset(srat, subset = pMT < 0.5)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 10)
  srat = FindNeighbors(srat, dims=1:10)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:10, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

uhlitz_epi_new=process_seurat2(uhlitz_epi)
uhlitz_epi_new@meta.data$new_ident=paste0(uhlitz_epi_new@meta.data$sample_origin, "_", uhlitz_epi_new@meta.data$cell_type_epi_simple)

comm_genes=intersect(rownames(rasa_updated_srat_processed), rownames(uhlitz_epi_new))
rasa_ref_train=trainModel(rasa_updated_srat_processed@assays$RNA@counts[comm_genes,], rasa_updated_srat_processed@meta.data$annotation_2, workers=NULL)


rasa_ref_train_uhlitz=predictSimilarity(rasa_ref_train, uhlitz_epi_new@assays$RNA@counts[comm_genes,], uhlitz_epi_new@meta.data$new_ident)

similarityHeatmap(rasa_ref_train_uhlitz)

low_ct_filt=subset(rasa_updated_srat_processed_g1,annotation_2%in%c("Fetal_Microfold","Fetal_Paneth",
                                                                    "Fetal_CLDN10+ cells","Fetal_Tuft"), invert=T)

rasa_ref_g1_train=trainModel(low_ct_filt@assays$RNA@counts[comm_genes,],
                          low_ct_filt@meta.data$annotation_2, workers=NULL)
rasa_ref_g1_train_uhlitz=predictSimilarity(rasa_ref_g1_train, uhlitz_epi_new@assays$RNA@counts[comm_genes,], uhlitz_epi_new@meta.data$new_ident)
similarityHeatmap(rasa_ref_g1_train_uhlitz)

gse146771=read.delim('GSE146771/GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt', sep = " ")
gse146771_meta=read.delim('GSE146771/GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt')

gse146771=as.matrix(gse146771)
gse146771=Matrix(gse146771,sparse = T)


gse_comm_genes=intersect(rownames(rasa_updated_srat_processed), rownames(gse146771))
rasa_ref_train_gse=trainModel(rasa_updated_srat_processed@assays$RNA@counts[gse_comm_genes,], rasa_updated_srat_processed@meta.data$annotation_2, workers=NULL)


rasa_ps_gse=predictSimilarity(rasa_ref_train_gse, gse146771[gse_comm_genes,], gse146771_meta$Global_Cluster)

similarityHeatmap(rasa_ps_gse)

