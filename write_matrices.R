library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(reshape2)
source('logisticRegression.R')



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


paths=paste0("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/TCGA/frag/TCGA_", rownames(mDat_col), ".tsv")

write.table(paths, "colon_bulk_paths.txt", sep = "\t", row.names = F, col.names = F, quote = F)

comm_g=intersect(rownames(epi_seurat_nocc), obs$index)
obs2=obs
rownames(obs2)=obs2$index
obs2=obs2[comm_g,]

epi_mtx=epi_seurat_nocc@assays$RNA@counts
epi_mtx=epi_mtx[rownames(obs2),]
rownames(epi_mtx)=obs2$gene_ids

colnames(epi_mtx)=paste0(epi_seurat_nocc@meta.data$Cluster,":",colnames(epi_mtx))

writeMM(epi_mtx, "cellsig_gut/sc_regev/epi_mtx.mtx")
write.table(rownames(epi_mtx), "cellsig_gut/sc_regev/epi_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(epi_mtx), "cellsig_gut/sc_regev/epi_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


rasa_mtx=lif_nocc_seurat@assays$RNA@counts
rownames(rasa_mtx)=obs$gene_ids

colnames(rasa_mtx)=paste0(lif_nocc_seurat@meta.data$exp_annot,":",colnames(rasa_mtx))

writeMM(rasa_mtx, "cellsig_gut/sc_rasa/rasa.mtx")
write.table(rownames(rasa_mtx), "cellsig_gut/sc_rasa/rasa_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(rasa_mtx), "cellsig_gut/sc_rasa/rasa_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


comb_epi_rasa_rn=intersect(rownames(rasa_mtx), rownames(epi_mtx))

combined_mtx=cbind(rasa_mtx[comb_epi_rasa_rn,], epi_mtx[comb_epi_rasa_rn,])
combined_mtx=as.matrix(combined_mtx)
combined_mtx=Matrix(combined_mtx, sparse = TRUE) 

writeMM(combined_mtx, "cellsig_gut/sc_combined/combined.mtx")
write.table(rownames(combined_mtx), "cellsig_gut/sc_combined/combined_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(combined_mtx), "cellsig_gut/sc_combined/combined_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

source('cellsig_gut/cellSignalAnalysis.R')
tgt='cellsig_gut/out_regev/OutRun_fitExposures.tsv'
tgt2='cellsig_gut/out_rasa/OutRun_fitExposures.tsv'
combined_noCTA='cellsig_gut/out_combined_no_CTA/OutRun_fitExposures.tsv'
combined_CTA='cellsig_gut/out_combined/OutRun_fitExposures.tsv'

fit=normaliseExposures(tgt)
fit2=normaliseExposures(tgt2)
fit_combined_noCTA=normaliseExposures(combined_noCTA)
fit_combined_CTA=normaliseExposures(combined_CTA)

mDat_col$gdc_cases.samples.sample_type=factor(mDat_col$gdc_cases.samples.sample_type, levels = c( "Solid Tissue Normal",
                                                                                                  "Primary Tumor",
                                                                                                  "Recurrent Tumor",
                                                                                                  "Metastatic" ))
plotExposures(fit=fit, show_column_names=F, show_column_dend=F, column_split=mDat_col$gdc_cases.samples.sample_type)
plotExposures(fit=fit2, show_column_names=F, show_column_dend=F, column_split=mDat_col$gdc_cases.samples.sample_type)

plotExposures(fit=fit_combined_noCTA, show_column_names=F, cluster_column_slices=F, column_split=mDat_col$gdc_cases.samples.sample_type)
plotExposures(fit=fit_combined_noCTA, show_column_names=F, cluster_column_slices=F, column_split=mDat_col$gdc_cases.samples.sample_type)

rel_meta=data.frame(stage=mDat_col$xml_stage_event_pathologic_stage, location=mDat_col$tissue,
                    vital_status=mDat_col$cgc_case_vital_status, 
                    disease_type=mDat_col$cgc_file_disease_type,
                    percent_tum_nuc=mDat_col$cgc_slide_percent_tumor_nuclei,
                    percent_tum_cell=mDat_col$cgc_slide_percent_tumor_cells,
                    percent_stromal=mDat_col$cgc_slide_percent_stromal_cells,
                    percent_normal=mDat_col$cgc_slide_percent_normal_cells,
                    percent_necrosis=mDat_col$cgc_slide_percent_necrosis,
                    pathologic_t=mDat_col$cgc_case_pathologic_t,
                    pathologic_n=mDat_col$cgc_case_pathologic_n,
                    gender=mDat_col$xml_gender,
                    age_at_diagnosis=mDat_col$xml_age_at_initial_pathologic_diagnosis, 
                    case_tumor_status=mDat_col$cgc_case_tumor_status,
                    follow_up_tumour_status=mDat_col$cgc_follow_up_tumor_status,
                    follow_up_vital_status=mDat_col$cgc_follow_up_vital_status,
                    case_histological_diagnosis=mDat_col$cgc_case_histological_diagnosis,
                    new_tum_after_treatment=mDat_col$cgc_follow_up_new_tumor_event_after_initial_treatment,
                    neopl_status=mDat_col$xml_person_neoplasm_cancer_status, 
                    has_new_tum_events=mDat_col$xml_has_new_tumor_events_information, 
                    lymphnodes_positive=mDat_col$xml_number_of_lymphnodes_positive_by_he, 
                    pretr_cea_level=mDat_col$xml_preoperative_pretreatment_cea_level, 
                    venous_invasion=mDat_col$xml_venous_invasion, 
                    lymphatic_invasion=mDat_col$xml_lymphatic_invasion, 
                    perineural_invasion=mDat_col$xml_perineural_invasion_present, 
                    kras_mutation=mDat_col$xml_kras_mutation_found, 
                    braf_analysis=mDat_col$xml_braf_gene_analysis_result)
rel_meta=data.frame(icdo=mDat_col$xml_icd_o_3_histology)

plotExposures(fit=fit_combined_noCTA, show_column_names=F, 
              column_split=mDat_col$gdc_cases.samples.sample_type, cluster_column_slices=F, show_column_dend=F)
hmCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
exposureScale=c(0,0.5)
hmColObj = circlize::colorRamp2(seq(exposureScale[1],exposureScale[2],length.out=length(hmCols)),hmCols)
botAnno = HeatmapAnnotation(df = rel_meta,
                            annotation_name_side = 'left')
row_ord_noCTA=c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
                "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
                "Immature Goblet","Secretory TA","fetal_Goblet cell","fetal_Erythocyte","fetal_TA",
                "fetal_T/NK", "fetal_Mesothelium","fetal_BEST4 enterocyte","fetal_B/Plasma","fetal_Enteroendocrine",
                "fetal_Enterocyte","fetal_Neuronal","fetal_Tuft cell","fetal_Mesenchymal","fetal_Endothelial",
                "fetal_crypt","fetal_Myeloid","Intercept")
row_ord_CTA=c("Stem","Cycling TA","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
                "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
                "Immature Goblet","Secretory TA","fetal_Goblet cell","fetal_Erythocyte","fetal_TA",
                "fetal_T/NK", "fetal_Mesothelium","fetal_BEST4 enterocyte","fetal_B/Plasma","fetal_Enteroendocrine",
                "fetal_Enterocyte","fetal_Neuronal","fetal_Tuft cell","fetal_Mesenchymal","fetal_Endothelial",
                "fetal_crypt","fetal_Myeloid","Intercept")
Heatmap(fit_combined_noCTA$exposures, col=hmColObj, cluster_rows=F, show_column_names=F, 
        column_split=mDat_col$gdc_cases.samples.sample_type, cluster_column_slices=F, show_column_dend=F,
        row_order = row_ord_noCTA, bottom_annotation=botAnno)
Heatmap(fit_combined_CTA$exposures, col=hmColObj, cluster_rows=F, show_column_names=F, 
        column_split=mDat_col$gdc_cases.samples.sample_type, cluster_column_slices=F, show_column_dend=F,
        row_order = row_ord_CTA, bottom_annotation=botAnno)


inter_CTA=data.frame(intercept=c(as.data.frame(t(fit_combined_CTA$exposures))$Intercept,
                       as.data.frame(t(fit_combined_noCTA$exposures))$Intercept),
                     data=c(rep("CTA", 715),rep("no CTA", 715)))
inter_stem=data.frame(intercept=c(as.data.frame(t(fit_combined_CTA$exposures))$Stem,
                                 as.data.frame(t(fit_combined_noCTA$exposures))$Stem),
                     data=c(rep("CTA", 715),rep("no CTA", 715)))


ggplot(data = inter_CTA, aes(x=data, y=intercept)) +geom_boxplot()


ggplot(data = inter_stem, aes(x=data, y=intercept)) +geom_boxplot()

inter_tissue=data.frame(intercept=as.data.frame(t(fit_combined_noCTA$exposures))$Intercept, tissue=mDat_col$gdc_cases.samples.sample_type)
inter_tissue_sub=inter_tissue[inter_tissue$tissue%in%c("Solid Tissue Normal", "Primary Tumor"),]

ggplot(data = inter_tissue_sub, aes(x=tissue, y=intercept)) +geom_boxplot()

master_df=as.data.frame(t(fit_combined_noCTA$exposures))
master_df$sample=rownames(master_df)

master_df$sample_type=as.character(mDat_col$gdc_cases.samples.sample_type)



master_df_melted=melt(master_df, id.vars = c("sample", "sample_type"))

master_df_melted=master_df_melted[master_df_melted$sample_type%in%c("Solid Tissue Normal", "Primary Tumor"),]

master_df_melted=master_df_melted[master_df_melted$variable%in%c("Stem",
                                                                 "Intercept", "fetal_crypt",
                                                                 "TA1", "TA2",
                                                                 "Immature Enterocytes 1", 
                                                                 "Immature Enterocytes 2",
                                                                 "Enterocytes"),]

ggplot(master_df_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +scale_fill_manual(values=c("red","lightblue","blue","darkblue","green","grey"))

clinical_2=read.delim('clinical.tsv')
                      