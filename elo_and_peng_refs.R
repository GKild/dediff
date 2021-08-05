######Peng lung reference#########

source('process_seurat.R')

flung_mtx=readMM('/lustre/scratch117/casm/team274/gk14/Dediff/FetalLungData/meta_data/matrix.mtx')
flung_cells=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/FetalLungData/meta_data/var.csv')
flung_genes=read.csv('/lustre/scratch117/casm/team274/gk14/Dediff/FetalLungData/meta_data/obs.csv')

colnames(flung_mtx)=flung_cells$X
rownames(flung_mtx)=flung_genes$X

rownames(flung_cells)=flung_cells$X

flung_srat=CreateSeuratObject(flung_mtx, meta.data = flung_cells)

flung_srat_processed=process_seurat(flung_srat)
#######Elo Lung reference#########

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

elo_lung_srat@meta.data$annot=paste0("Adult_", elo_lung_srat$annot)
flung_srat_processed@meta.data$annot=paste0("Fetal_", flung_srat_processed$celltype)
elo_lung_mtx=elo_lung_srat@assays$RNA@counts
peng_lung_mtx=flung_srat_processed@assays$RNA@counts


inter=intersect(rownames(elo_lung_mtx), rownames(peng_lung_mtx))

colnames(elo_lung_mtx)=paste0(elo_lung_srat@meta.data$annot,":",colnames(elo_lung_mtx))
colnames(peng_lung_mtx)=paste0(flung_srat_processed@meta.data$annot,":", colnames(peng_lung_mtx)) 


fetal_and_adult_lung=cbind(peng_lung_mtx[inter, ], elo_lung_mtx[inter, ])
ct=c(flung_srat_processed@meta.data$annot,elo_lung_srat@meta.data$annot)

merged_srat=CreateSeuratObject(fetal_and_adult_lung)
Idents(merged_srat)=merged_srat_processed



merged_srat_processed=process_seurat(merged_srat)

merged_srat_processed@meta.data$fetal_or_adult=c(rep("fetal", 52517), rep("adult", 129340))

Idents(merged_srat_processed)=merged_srat_processed@meta.data$annot

DimPlot(merged_srat_processed, group.by = "fetal_or_adult")

rownames(elo_lung_genes)=elo_lung_genes$X

g2=elo_lung_genes[inter,]



merged_srat_processed_g1=subset(merged_srat_processed, Phase=="G1")


#Matts cor
dat = merged_srat_processed@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),merged_srat_processed@meta.data$annot),
                           function(e) rowMeans(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)
#################################
#Gerda's cor
avg_exp=AverageExpression(merged_srat_processed, slot="scale.data")

cordat=cor(avg_exp$RNA)
Heatmap(cordat)
#################################

mtx=merged_srat_processed@assays$RNA@counts
rownames(mtx)=g2$gene_ids.1

mtx_g1=merged_srat_processed_g1@assays$RNA@counts
rownames(mtx_g1)=g2$gene_ids.1

saveRDS(merged_srat_processed, '/lustre/scratch117/casm/team274/gk14/Dediff/elo_and_peng_sc.rds')


writeMM(mtx, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_elo_peng/sc_elo_peng.mtx")
write.table(rownames(mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_elo_peng/sc_elo_peng_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_elo_peng/sc_elo_peng_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


writeMM(mtx_g1, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_elo_peng_nocc/sc_elo_peng.mtx")
write.table(rownames(mtx_g1), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_elo_peng_nocc/sc_elo_peng_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(mtx_g1), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_lung/sc_elo_peng_nocc/sc_elo_peng_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')


lung_cc_fit=normaliseExposures('cellsig_lung/out_elo_peng/OutRun_fitExposures.tsv')
lung_nocc=normaliseExposures('cellsig_lung/out_elo_peng_nocc/OutRun_fitExposures.tsv')

rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
obs=read.csv("onecs/obs.csv")
mDat = colData(rc2)
mDat$tissue = as.character(mDat$xml_tumor_tissue_site)

table(mDat$tissue)

mDat_lung=mDat[which(mDat$tissue=="Lung"),]

lung_meta=data.frame(sample=c(rownames(mDat_lung), "GSM1101684", "GSM1101685", "GSM1101687", "GSM1101693", "GSM1101699",
                              "GSM1101708"),sample_type=c(mDat_lung$cgc_sample_sample_type, rep("fetal_bulk_lung", 6)),
                     disease_type=c(mDat_lung$cgc_file_disease_type,rep("fetal_bulk_lung", 6)),
                     histol_diag= c(mDat_lung$cgc_case_histological_diagnosis, rep("fetal_bulk_lung", 6)))

lung_meta$broad_colsplit=lung_meta$sample_type
lung_meta$broad_colsplit[which(lung_meta$broad_colsplit=="Primary Tumor")]=lung_meta$disease_type[which(lung_meta$broad_colsplit=="Primary Tumor")]

lung_meta$narrow_colsplit=lung_meta$sample_type
lung_meta$narrow_colsplit[which(lung_meta$narrow_colsplit=="Primary Tumor")]=lung_meta$histol_diag[which(lung_meta$narrow_colsplit=="Primary Tumor")]


plotExposures(lung_cc_fit, show_column_names=F,column_split=lung_meta$narrow_colsplit,
              column_title_rot=90, cluster_column_slices=F)
plotExposures(lung_nocc, show_column_names=F,
              column_split=lung_meta$narrow_colsplit,
              column_title_rot=90, cluster_column_slices=F)

mDat_col=readRDS('mDat_gut.rds')
split_cols=c(as.character(mDat_col$gdc_cases.samples.sample_type), rep("fetal_bulk", 8))
split_cols=factor(split_cols, levels = c("Solid Tissue Normal","Primary Tumor","fetal_bulk",
                                         "Recurrent Tumor","Metastatic"))

gut_early_emb_fit=normaliseExposures('cellsig_early_emb/c_genes_out/OutRun_fitExposures.tsv')
plotExposures(gut_early_emb_fit, show_column_names=F, column_split=split_cols)


gut_exp=gut_early_emb_fit$exposures
