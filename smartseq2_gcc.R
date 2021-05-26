#read the log tpm data and metadata for HCC samples
library(Matrix)
library(Seurat)
library(ComplexHeatmap)
source('logisticRegression.R')
workers=NULL
hcc_sc=read.table('smart_seq2_hcc/HCC_log_tpm_expression_matrix.txt', sep = "\t", header = T, row.names = 1)
hcc_meta=read.table('smart_seq2_hcc/HCC_cell_metadata.txt', sep = "\t", header = T)


hcc_sc_tpm=apply(hcc_sc, 2, function(x){x^2})

hcc_sc_tpm_mtx=Matrix(hcc_sc_tpm, sparse = T)

identical(colnames(hcc_sc_tpm_mtx), hcc_meta$name)

table(hcc_meta$cell_type)


table(merged_liver_srat_processed@meta.data$annot)

comm_genes=intersect(rownames(hcc_sc_tpm_mtx), rownames(merged_liver_srat_processed))

liver_model_cycling=trainModel(merged_liver_srat_processed@assays$RNA@counts[comm_genes,],
                               merged_liver_srat_processed@meta.data$annot, workers=NULL)

ps_hcc=predictSimilarity(liver_model_cycling, hcc_sc_tpm_mtx[comm_genes,], hcc_meta$cell_type)

similarityHeatmap(ps_hcc)

just_malignant_cells_meta=hcc_meta[grep('_Tumor', hcc_meta$cell_type),]
just_malignant_cells_meta$type_celltype=paste0(just_malignant_cells_meta$cell_type, "_", just_malignant_cells_meta$HCC_type)

just_malignant_mtx=hcc_sc_tpm_mtx[,just_malignant_cells_meta$name]

ps_hcc_tum=predictSimilarity(liver_model_cycling, just_malignant_mtx[comm_genes,], just_malignant_cells_meta$type_celltype)

similarityHeatmap(ps_hcc_tum)



liver_model_nocc=trainModel(merged_liver_srat_processed_nocc@assays$RNA@counts[comm_genes,],
                            merged_liver_srat_processed_nocc@meta.data$annot, workers=NULL)

ps_hcc_tum_nocc=predictSimilarity(liver_model_nocc, just_malignant_mtx[comm_genes,], just_malignant_cells_meta$type_celltype)
similarityHeatmap(ps_hcc_tum_nocc)

DimPlot(merged_liver_srat_processed,
        cells.highlight = colnames(merged_liver_srat_processed)[grep('Adult_Cholangiocytes', merged_liver_srat_processed@meta.data$annot)])

ps_hcc_nocc=predictSimilarity()


