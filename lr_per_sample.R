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
len=read.table('cellSignalAnalysis/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
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

just_large_cells=rownames(colon_seurat@meta.data)[which(colon_seurat@meta.data$Organ%in%c("TCL", "ACL", "SCL", "FLI", "REC", "DCL", "APD"))]

li_seurat=subset(colon_seurat, cells=just_large_cells)

li_train=trainModel(li_seurat@assays$RNA@counts, li_seurat@meta.data$exp_annot, workers=NULL, minCells = 300)
li_confusion_pred=predictSimilarity(li_train, li_seurat@assays$RNA@counts, li_seurat@meta.data$exp_annot)


li_seurat@meta.data$Donor=as.character(li_seurat@meta.data$Donor)
li_seurat=subset(li_seurat, cells=rownames(li_seurat@meta.data)[which(li_seurat@meta.data$Phase=="G1")])
remove_small_clusts=function(x){
  keep_idents=names(which(table(x@meta.data$exp_annot)>=5))
  x=subset(x, cells=rownames(x@meta.data)[which(x@meta.data$exp_annot%in%keep_idents)])
  return(x)
}

for (i in unique(as.character(li_seurat@meta.data$Donor))) {
  srat=subset(li_seurat, cells= rownames(li_seurat@meta.data)[which(li_seurat@meta.data$Donor==i)])
  srat=remove_small_clusts(srat)
  li_mtx=srat@assays$RNA@counts
  rownames(li_mtx)=obs$gene_ids
  tr=trainModel(li_mtx[comm_ensg,], srat@meta.data$exp_annot, workers=NULL)
  ps=predictSimilarity(tr,crc_bulk_mat, lengths = gene_lengths)
  pdf(paste0("logreg_nocc_",i, ".pdf"), width = 10, height = 10)
  print(similarityHeatmap(ps, row_split=mDat_col$gdc_cases.samples.sample_type))
  dev.off()
  
}


for (i in unique(as.character(li_seurat@meta.data$Donor))) {
  srat=subset(li_seurat, cells= rownames(li_seurat@meta.data)[which(li_seurat@meta.data$Donor==i)])
  srat=remove_small_clusts(srat)
  li_mtx=srat@assays$RNA@counts
  tr=trainModel(li_mtx, srat@meta.data$exp_annot, workers=NULL)
  pdf(paste0(i,"ag_all_samples.pdf"))
  for (x in unique(as.character(li_seurat@meta.data$Donor))) {
    srat=subset(li_seurat, cells= rownames(li_seurat@meta.data)[which(li_seurat@meta.data$Donor==x)])
    ps=predictSimilarity(tr,srat@assays$RNA@counts, srat@meta.data$exp_annot)
    inter=intersect(rownames(ps), colnames(ps))
    row_ord=c(inter, rownames(ps)[which(!rownames(ps)%in%inter)])
    col_ord=c(inter, colnames(ps)[which(!colnames(ps)%in%inter)])
    print(similarityHeatmap(ps, row_title=x, row_order=row_ord, column_order=col_ord))}
  dev.off()
}

