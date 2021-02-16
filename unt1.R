library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
source('logisticRegression.R')

colon_seurat=readRDS("/lustre/scratch117/casm/team274/gk14/human_colon.rds")

colon_seurat@meta.data$age_cat="fetal"
colon_seurat@meta.data$age_cat[which(colon_seurat@meta.data$Age%in%c("10","12","4","6","9"))]="pediatric"
colon_seurat@meta.data$age_cat[which(colon_seurat@meta.data$Age%in%c("A26", "A32", "A34", "A38", "A39"))]="adult"
colon_seurat@meta.data$age_annot=paste0(colon_seurat$age_cat, "_", colon_seurat@meta.data$annotation_rasa)

rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
obs=read.csv("onecs/obs.csv")
mDat = colData(rc2)

mDat$UniqueSampleID = rownames(mDat)
mDat$age = mDat$cgc_case_age_at_diagnosis*12
mDat$biopsy = mDat$gdc_cases.samples.sample_type
mDat$tissue = as.character(mDat$xml_tumor_tissue_site)

mDat_col=mDat[which(mDat$tissue%in%c("Rectum", "Colon")),]
mDat_col=mDat_col[,colnames(mDat_col)[which(!colnames(mDat_col)%in%names(which(apply(mDat_col, 2, function(x){sum(is.na(x))})==715)))]]

bulk_mat=assays(rc2)$counts
rownames(bulk_mat)=sapply(strsplit(rownames(bulk_mat), "\\."), "[", 1)
comm_ensg=intersect(rownames(bulk_mat), obs$gene_ids)

crc_bulk_mat=bulk_mat[comm_ensg, mDat_col$UniqueSampleID]
rownames(len)=make.unique(as.character(len$geneName))
len=len[which(rownames(len)%in%comm_ensg),]
gene_lengths=len$geneLengths



kylie_annot=read.csv("Full_ob_230720.csv")

colon_seurat@meta.data$cell_group=kylie_annot$Cell_group[match(colnames(colon_seurat), kylie_annot$index)]

keep_cells=colnames(colon_seurat)[which(!as.character(colon_seurat@meta.data$cell_group)%in%c("doublet", "empty"))]

colon_seurat=subset(colon_seurat, cells=keep_cells)

colon_seurat@meta.data$cell_group=as.character(colon_seurat@meta.data$cell_group)
colon_seurat@meta.data$cell_group_age=paste0(as.character(colon_seurat$age_cat), "_", as.character(colon_seurat@meta.data$cell_group))


colon_seurat@meta.data$exp_annot=as.character(colon_seurat@meta.data$cell_group_age)
colon_seurat@meta.data$exp_annot[which(colon_seurat@meta.data$cell_group=="Epithelial")]=colon_seurat@meta.data$age_annot[which(colon_seurat@meta.data$cell_group=="Epithelial")]
DimPlot(colon_seurat, group.by = "exp_annot")

sort(table(colon_seurat@meta.data$exp_annot))
colon_seurat@meta.data$exp_annot[which(colon_seurat@meta.data$exp_annot=="adult_BEST1 enterocyte")]="adult_Enterocyte"

sc_mtx=colon_seurat@assays$RNA@counts

rownames(sc_mtx)=obs$gene_ids

train_for_bulk=trainModel(sc_mtx[comm_ensg,], colon_seurat@meta.data$exp_annot, workers=NULL, maxCells = 20000)
predict_bulk=predictSimilarity(train_for_bulk, crc_bulk_mat, lengths = gene_lengths)


ann_ord=c("fetal_Enterocyte","fetal_Goblet cell","fetal_crypt",
          "fetal_Mesothelium","fetal_BEST4 enterocyte","fetal_TA",
          "fetal_Enteroendocrine","fetal_Tuft cell", "fetal_LI Epithelial progenitor",
          "fetal_Epithelial progenitor","fetal_Germ cells", "pediatric_Enterocyte",
          "pediatric_crypt/TA","pediatric_Tuft cell","pediatric_Goblet cell",
          "pediatric_CHGA Enteroendocrine", "pediatric_NTS Enteroendocrine","pediatric_BEST4 enterocyte",
          "pediatric_Paneth","adult_Enterocyte","adult_crypt",
          "adult_BEST4 enterocyte", "adult_TA", "adult_Goblet cell",
          "adult_Tuft cell","adult_Enteroendocrine","adult_Paneth","adult_Mesothelium", 
          "fetal_Endothelial","fetal_Mesenchymal","fetal_Erythocyte","fetal_Neuronal", "fetal_Myeloid",
          "fetal_T/NK","fetal_B/Plasma","pediatric_B/Plasma","pediatric_T/NK","pediatric_Mesenchymal",
          "pediatric_Myeloid" ,"pediatric_Endothelial", "pediatric_Neuronal","adult_Endothelial","adult_Mesenchymal",
          "adult_B/Plasma","adult_T/NK","adult_Myeloid","adult_Neuronal")


similarityHeatmap(predict_bulk, row_split=mDat_col$gdc_cases.samples.sample_type, column_order=ann_ord)


just_large_cells=rownames(colon_seurat@meta.data)[which(colon_seurat@meta.data$Organ%in%c("TCL", "ACL", "SCL", "FLI", "REC", "DCL", "APD"))]
just_small_cells=rownames(colon_seurat@meta.data)[which(!colon_seurat@meta.data$Organ%in%c("TCL", "ACL", "SCL", "FLI", "REC", "DCL", "APD", "FMLN", "MLN"))]

li_seurat=subset(colon_seurat, cells=just_large_cells)
si_seurat=subset(colon_seurat, cells=just_small_cells)

saveRDS(li_seurat, "large_intestine_rasa.rds")

li_seurat@meta.data$exp_annot[which(li_seurat@meta.data$exp_annot=="adult_BEST1 enterocyte")]="adult_Enterocyte"

si_seurat@meta.data$exp_annot[which(si_seurat@meta.data$exp_annot=="adult_BEST1 enterocyte")]="adult_Enterocyte"


li_train=trainModel(li_seurat@assays$RNA@counts, li_seurat@meta.data$exp_annot, workers=NULL)
li_confusion_pred=predictSimilarity(li_train, li_seurat@assays$RNA@counts, li_seurat@meta.data$exp_annot)

similarityHeatmap(li_confusion_pred, row_order=unique(li_seurat@meta.data$exp_annot),column_order=unique(li_seurat@meta.data$exp_annot))


li_train_canc=trainModel(li_seurat@assays$RNA@counts[comm_genes,], li_seurat@meta.data$exp_annot, workers=NULL)
li_pred_canc=predictSimilarity(li_train_canc, colon_matrix_sparse[comm_genes,], colon_canc_annot)
similarityHeatmap(li_pred_canc)


li_mtx=li_seurat@assays$RNA@counts

rownames(li_mtx)=obs$gene_ids

train_for_bulk=trainModel(li_mtx[comm_ensg,], li_seurat@meta.data$exp_annot, workers=NULL)
predict_bulk=predictSimilarity(train_for_bulk, crc_bulk_mat, lengths = gene_lengths)

similarityHeatmap(predict_bulk,row_split=mDat_col$gdc_cases.samples.sample_type)

li_seurat_nocc=subset(li_seurat, cells=colnames(li_seurat)[which(li_seurat@meta.data$Phase=="G1")])
li_mtx_nocc=li_seurat_nocc@assays$RNA@counts

rownames(li_mtx_nocc)=obs$gene_ids

train_for_bulk_nocc=trainModel(li_mtx_nocc[comm_ensg,], li_seurat_nocc@meta.data$exp_annot, workers=NULL)
predict_bulk_nocc=predictSimilarity(train_for_bulk_nocc, crc_bulk_mat, lengths = gene_lengths)

similarityHeatmap(predict_bulk_nocc,row_split=mDat_col$gdc_cases.samples.sample_type)


si_mtx=si_seurat@assays$RNA@counts

rownames(si_mtx)=obs$gene_ids

si_train_for_bulk=trainModel(si_mtx[comm_ensg,], si_seurat@meta.data$exp_annot, workers=NULL)
si_predict_bulk=predictSimilarity(si_train_for_bulk, crc_bulk_mat, lengths = gene_lengths)

similarityHeatmap(si_predict_bulk,row_split=mDat_col$gdc_cases.samples.sample_type)


