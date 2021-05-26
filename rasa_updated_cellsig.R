boardman_meta=read.table('boardman_data/boardman_meta.txt', sep = ",", header = T)
boardman_data=read.delim('boardman_data/boardman_data.tsv', sep = "\t", header = T)

boardman_rel=boardman_data[,c(4,5,grep("count", colnames(boardman_data)))]

colnames(boardman_rel)=gsub("_count", "", colnames(boardman_rel))

length(which(boardman_meta$sample%in%colnames(boardman_rel)))

colnames(boardman_rel)=gsub("Unk.", "", colnames(boardman_rel))

boardman_meta_new=data.frame(sample=colnames(boardman_rel)[3:81])
which(boardman_meta_new$sample%in%boardman_meta$sample)
boardman_meta_new$sample_type="x"

boardman_meta_new$sample_type=boardman_meta$sample_type[match(boardman_meta_new$sample, boardman_meta$sample)]

boardman_meta_new$sample_type[which(is.na(boardman_meta_new$sample_type))]="unknown"

boardman_colsplit=c(boardman_meta_new$sample_type, rep("fetal_bulk",8))
boardman_colsplit[grep("VILLOUS", boardman_colsplit)]="VILLOUS ADENOMA"


mDat_col=readRDS('tcga_meta.rds')
split_cols=c(as.character(mDat_col$gdc_cases.samples.sample_type), rep("fetal_bulk", 8))
split_cols_tcga=factor(split_cols, levels = c("Solid Tissue Normal","Primary Tumor","fetal_bulk",
                                         "Recurrent Tumor","Metastatic"))

row_ord=c( "Adult_Stem cells", "Fetal_Stem cells", "Adult_TA", "Fetal_TA","Adult_BEST4+ colonocyte","Adult_Microfold cell",
           "Adult_Colonocyte","Adult_BEST2+ Goblet cell","Adult_Enteroendocrine","Adult_Tuft", "Adult_Goblet cell",
           "Fetal_Distal progenitor","Fetal_Proximal progenitor", "Fetal_Goblet","Fetal_Paneth","Fetal_Enterocyte",
           "Fetal_Enteroendocrine", "Fetal_BEST4+ epithelial", "Fetal_Colonocyte",  "Fetal_Tuft" , "Intercept")

source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')

rasa_updated_boardman='cellsig_gut/out_rasa_updated_boardman/OutRun_fitExposures.tsv'
fit_boardman=normaliseExposures(rasa_updated_boardman)

plotExposures(fit_boardman, column_split=boardman_colsplit,column_title_rot=90, row_order=row_ord, 
              show_column_names=F,cluster_column_slices=F, show_column_dend=F)

rasa_updated_tcga='cellsig_gut/out_rasa_updated_tcga/OutRun_fitExposures.tsv'
fit_tcga=normaliseExposures(rasa_updated_tcga)

plotExposures(fit_tcga,column_title_rot=90, column_split=split_cols_tcga, row_order=row_ord, 
              show_column_names=F,cluster_column_slices=F, show_column_dend=F)



rasa_updated_g1_boardman='cellsig_gut/out_rasa_updated_g1_boardman/OutRun_fitExposures.tsv'
fit_g1_boardman=normaliseExposures(rasa_updated_g1_boardman)

plotExposures(fit_g1_boardman, column_split=boardman_colsplit,column_title_rot=90, row_order=row_ord, 
              show_column_names=F,cluster_column_slices=F, show_column_dend=F)

rasa_updated_g1_tcga='cellsig_gut/out_rasa_updated_g1_tcga/OutRun_fitExposures.tsv'
fit_g1_tcga=normaliseExposures(rasa_updated_g1_tcga)

plotExposures(fit_g1_tcga,column_title_rot=90, column_split=split_cols_tcga, row_order=row_ord, 
              show_column_names=F,cluster_column_slices=F, show_column_dend=F)


