library(Matrix)
library(Seurat)
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



fetal_lung_mtx=readMM('/lustre/scratch117/casm/team274/gk14/fetal_lung/matrix.mtx')
fetal_lung_genes=read.csv('/lustre/scratch117/casm/team274/gk14/fetal_lung/obs.csv')
fetal_lung_meta=read.csv('/lustre/scratch117/casm/team274/gk14/fetal_lung/var.csv')

rownames(fetal_lung_meta)=fetal_lung_meta$X


adult_lung_mtx=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/krasnov_lung/krasnow_hlca_10x_UMIs.csv')
adult_lung_meta=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/krasnov_lung/krasnow_hlca_10x_metadata.csv')
rownames(adult_lung_meta)=adult_lung_meta$X

rownames(adult_lung_mtx)=adult_lung_mtx$X
adult_lung_mtx=adult_lung_mtx[,2:ncol(adult_lung_mtx)]

adult_lung_mtx=as.matrix(adult_lung_mtx)
adult_lung_mtx=Matrix(adult_lung_mtx,sparse = T)


rownames(fetal_lung_mtx)=fetal_lung_genes$gene_short_name
colnames(fetal_lung_mtx)=fetal_lung_meta$X

lung_comm_genes=intersect(rownames(fetal_lung_mtx), rownames(adult_lung_mtx))


fetal_lung_mtx=fetal_lung_mtx[lung_comm_genes,]
adult_lung_mtx=adult_lung_mtx[lung_comm_genes,]


fetal_lung_srat=CreateSeuratObject(fetal_lung_mtx, meta.data = fetal_lung_meta)
fetal_lung_srat = NormalizeData(fetal_lung_srat)
fetal_lung_srat =FindVariableFeatures(fetal_lung_srat, selection.method = "vst", nfeatures = 2000)
fetal_lung_srat = CellCycleScoring(fetal_lung_srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
fetal_lung_srat = ScaleData(fetal_lung_srat)
fetal_lung_srat = RunPCA(fetal_lung_srat, npcs = 50)
fetal_lung_srat = FindNeighbors(fetal_lung_srat, dims=1:50)
fetal_lung_srat = FindClusters(fetal_lung_srat, resolution = 1)
fetal_lung_srat = RunUMAP(fetal_lung_srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)

adult_lung_srat=CreateSeuratObject(adult_lung_mtx, meta.data=adult_lung_meta)
processed_adult_lung=process_seurat(adult_lung_srat)

DimPlot(processed_adult_lung, group.by = "annot")


DimPlot(fetal_lung_srat, group.by = "annot")

fetal_lung_srat@meta.data$annot=fetal_lung_srat@meta.data$Main_cluster_name
fetal_lung_srat@meta.data$annot=paste0("Fetal_", fetal_lung_srat@meta.data$annot)

processed_adult_lung@meta.data$annot=processed_adult_lung@meta.data$free_annotation
processed_adult_lung@meta.data$annot=paste0("Adult_", processed_adult_lung@meta.data$annot)


merged_lung=merge(processed_adult_lung, fetal_lung_srat)

#Matts cor
dat = merged_lung@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),merged_lung@meta.data$annot),
                           function(e) rowMeans(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)
#################################
#Gerda's cor

scale_dat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  return(srat)
}

merged_lung=scale_dat(merged_lung)
Idents(merged_lung)=merged_lung@meta.data$annot
avg_exp=AverageExpression(merged_lung, slot="scale.data")

cordat=cor(avg_exp$RNA)
Heatmap(cordat)
#################################


rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
obs=read.csv("onecs/obs.csv")
mDat = colData(rc2)

mDat$tissue = as.character(mDat$xml_tumor_tissue_site)

table(mDat$tissue)

mDat_lung=mDat[which(mDat$tissue=="Lung"),]

paths=paste0("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/TCGA/frag/TCGA_", rownames(mDat_lung), ".tsv")

write.table(paths, "cellsig_lung/lung_tcga_paths.txt", sep = "\t", row.names = F, col.names = F, quote = F)

sort(table(mDat_lung$xml_histological_type))

table(mDat_lung$gdc_cases.samples.sample_type)

elo_lung=readMM('/lustre/scratch117/casm/team274/gk14/Dediff/elo_lung/matrix.mtx')
elo_lung_genes=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/elo_lung/obs.csv')
elo_lung_meta=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/elo_lung/var.csv')


rownames(elo_lung)=elo_lung_genes$X
colnames(elo_lung)=elo_lung_meta$X

rownames(elo_lung_meta)=elo_lung_meta$X

elo_lung_srat=CreateSeuratObject(elo_lung, meta.data = elo_lung_meta)

elo_lung_srat=subset(elo_lung_srat, Material%in%c('cells'))


elo_lung_srat=process_seurat(elo_lung_srat)

DimPlot(elo_lung_srat, group.by = "Celltypes_master_high")

epi_cts=elo_lung_srat@meta.data[which(elo_lung_srat@meta.data$Celltypes_master_higher=="Epithelia"),]

elo_lung_srat@meta.data$annot=elo_lung_srat@meta.data$Celltypes_master_high
elo_lung_srat@meta.data$annot[which(elo_lung_srat@meta.data$Celltypes_master_higher=="Epithelia")]=elo_lung_srat@meta.data$Celltypes_int[which(elo_lung_srat@meta.data$Celltypes_master_higher=="Epithelia")]
DimPlot(elo_lung_srat, group.by = "annot")


lung_comm_genes=intersect(rownames(fetal_lung_mtx), rownames(elo_lung))


fetal_lung_mtx=fetal_lung_mtx[lung_comm_genes,]
elo_lung_mtx=elo_lung[lung_comm_genes,]


fetal_lung_srat=CreateSeuratObject(fetal_lung_mtx, meta.data = fetal_lung_meta)
elo_lung_srat=CreateSeuratObject(elo_lung_mtx, meta.data = elo_lung_meta)

fetal_lung_srat@meta.data$annot=fetal_lung_srat@meta.data$Main_cluster_name
fetal_lung_srat@meta.data$annot=paste0("Fetal_", fetal_lung_srat@meta.data$annot)

fetal_lung_srat@meta.data$annot

elo_lung_srat=subset(elo_lung_srat, Material%in%c('cells'))
elo_lung_srat@meta.data$annot=elo_lung_srat@meta.data$Celltypes_master_high
elo_lung_srat@meta.data$annot[which(elo_lung_srat@meta.data$Celltypes_master_higher=="Epithelia")]=elo_lung_srat@meta.data$Celltypes_int[which(elo_lung_srat@meta.data$Celltypes_master_higher=="Epithelia")]
elo_lung_srat@meta.data$annot=paste0("Adult_", elo_lung_srat@meta.data$annot)


merged_lung=merge(fetal_lung_srat, elo_lung_srat)


#Gerda's cor

scale_dat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  return(srat)
}

merged_lung=scale_dat(merged_lung)
Idents(merged_lung)=merged_lung@meta.data$annot
avg_exp=AverageExpression(merged_lung, slot="scale.data")

cordat=cor(avg_exp$RNA)
Heatmap(cordat)
#################################

#Matts cor
dat = merged_lung@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),merged_lung@meta.data$annot),
                           function(e) rowMeans(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)


rownames(fetal_lung_genes)=fetal_lung_genes$gene_short_name

fetal_lung_genes_overlap=fetal_lung_genes[lung_comm_genes,]
ens_ids=sapply(strsplit(fetal_lung_genes_overlap$gene_id, "\\."), "[[", 1)


write_cellsig_files=function(srat, ensembl_vector, out_location){
  dir.create(out_location)
  mat=srat@assays$RNA@counts
  rownames(mat)=ensembl_vector
  colnames(mat)=paste0(srat@meta.data$annot,":",colnames(mat))
  writeMM(mat, paste0(out_location, "lung.mtx"))
  write.table(rownames(mat), paste0(out_location, "lung_rowNames.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)
  write.table(colnames(mat), paste0(out_location, "lung_columnNames.tsv"), sep = "\t", quote = F, col.names = F, row.names = F)
  
}


write_cellsig_files(merged_lung, ens_ids, '/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_lung/')


merged_lung_nocc=subset(merged_lung, Phase=="G1")

write_cellsig_files(merged_lung_nocc, ens_ids, '/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_lung_nocc/')

merged_lung_nodiv=subset(merged_lung, idents=c("Adult_Dividing_AT2","Adult_Dividing_Basal"),invert=T)

write_cellsig_files(merged_lung_nodiv, ens_ids, '/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_lung_nodiv/')

source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')

lung_tcga='cellsig_lung/out_lung/OutRun_fitExposures.tsv'
fit_lung_tcga=normaliseExposures(lung_tcga)


plotExposures(fit_lung_tcga, show_column_names=F, column_split=lung_meta$narrow_colsplit, column_title_rot=90, cluster_column_slices=F)


lung_meta=data.frame(sample=rownames(mDat_lung),sample_type=mDat_lung$cgc_sample_sample_type, disease_type=mDat_lung$cgc_file_disease_type, histol_diag= mDat_lung$cgc_case_histological_diagnosis)

lung_meta$broad_colsplit=lung_meta$sample_type
lung_meta$broad_colsplit[which(lung_meta$broad_colsplit=="Primary Tumor")]=lung_meta$disease_type[which(lung_meta$broad_colsplit=="Primary Tumor")]

lung_meta$narrow_colsplit=lung_meta$sample_type
lung_meta$narrow_colsplit[which(lung_meta$narrow_colsplit=="Primary Tumor")]=lung_meta$histol_diag[which(lung_meta$narrow_colsplit=="Primary Tumor")]


lung_meta$narrow_colsplit=factor(lung_meta$narrow_colsplit, levels=c("Solid Tissue Normal", 
                                                                     "Lung Adenocarcinoma- Not Otherwise Specified (NOS)",
                                                                     "Lung Adenocarcinoma Mixed Subtype",
                                                                     "Lung Papillary Adenocarcinoma",
                                                                     "Lung Acinar Adenocarcinoma",
                                                                     "Lung Solid Pattern Predominant Adenocarcinoma",
                                                                     "Lung Mucinous Adenocarcinoma",
                                                                     "Lung Micropapillary Adenocarcinoma",
                                                                     "Lung Clear Cell Adenocarcinoma",
                                                                     "Lung Signet Ring Adenocarcinoma",
                                                                     "Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)",
                                                                     "Lung Basaloid Squamous Cell Carcinoma",
                                                                     "Lung Papillary Squamous Cell Carcinoma",
                                                                     "Lung Small Cell Squamous Cell Carcinoma",
                                                                     "Lung Bronchioloalveolar Carcinoma Nonmucinous",
                                                                     "Mucinous (Colloid) Carcinoma", 
                                                                     "Lung Bronchioloalveolar Carcinoma Mucinous",
                                                                     "Recurrent Tumor"))

lung_tcga_nocc='cellsig_lung/out_lung_nocc/OutRun_fitExposures.tsv'
fit_lung_tcga_nocc=normaliseExposures(lung_tcga_nocc)


plotExposures(fit_lung_tcga_nocc, show_column_names=F, column_split=lung_meta$narrow_colsplit, column_title_rot=90, cluster_column_slices=F, row_order=rownames(fit_lung_tcga$exposures))

lung_tcga_nodiv='cellsig_lung/out_lung_nodiv/OutRun_fitExposures.tsv'
fit_lung_tcga_nodiv=normaliseExposures(lung_tcga_nodiv)


plotExposures(fit_lung_tcga_nodiv, show_column_names=F, column_split=lung_meta$narrow_colsplit, column_title_rot=90, cluster_column_slices=F, row_order=rownames(fit_lung_tcga$exposures)[c(1:4,6:9, 11:40)])

out_rasa_final='cellsig_gut/out_rasa_final_tcga/OutRun_fitExposures.tsv'
fit_out_rasa_final=normaliseExposures(out_rasa_final)

rownames(fit_lung_tcga$exposures)[c(1:4,6:9, 11:40)]
plotExposures(fit_out_rasa_final, show_column_names=F)


fetal_bulk_lung='cellsig_lung/out_bulk_lung_nodiv/OutRun_fitExposures.tsv'
fit_fetal_bulk_lung=normaliseExposures(fetal_bulk_lung)


plotExposures(fit_fetal_bulk_lung, row_order=rownames(fit_lung_tcga$exposures)[c(1:4,6:9, 11:40)])


