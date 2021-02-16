######################################################################
#don't need need the commented lines once I've extracted CRC metadata# 
######################################################################

# rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
# len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
# obs=read.csv("onecs/obs.csv")
# mDat = colData(rc2)
# 
# mDat$UniqueSampleID = rownames(mDat)
# mDat$age = mDat$cgc_case_age_at_diagnosis*12
# mDat$biopsy = mDat$gdc_cases.samples.sample_type
# mDat$tissue = as.character(mDat$xml_tumor_tissue_site)
# 
# mDat_col=mDat[which(mDat$tissue%in%c("Rectum", "Colon")),]
# mDat_col=mDat_col[,colnames(mDat_col)[which(!colnames(mDat_col)%in%names(which(apply(mDat_col, 2, function(x){sum(is.na(x))})==715)))]]
#saveRDS(mDat_col,'tcga_meta.rds')

mDat_col=readRDS('tcga_meta.rds')
mDat_gtex=read.table('gtex_gut_tissues.txt', sep = "\t", header = T)

split_cols=c(as.character(mDat_col$gdc_cases.samples.sample_type), rep("fetal_bulk", 8), mDat_gtex$tissue)
split_cols=factor(split_cols, levels = c("Solid Tissue Normal","Colon - Transverse","Colon - Sigmoid",
                                         "Small Intestine - Terminal Ileum","Primary Tumor","fetal_bulk",
                                         "Recurrent Tumor","Metastatic"))
combined_noCTA='cellsig_gut/out_with_gtex/OutRun_fitExposures.tsv'
fit_combined_noCTA=normaliseExposures(combined_noCTA)

meta_for_all=data.frame(sample=colnames(fit_combined_noCTA$exposures),
                        sample_type=split_cols,
                        stage=c(mDat_col$xml_stage_event_pathologic_stage, rep(NA,8)))

unique(meta_for_all$stage)

meta_for_all$stage[which(meta_for_all$stage%in%c("Stage I","Stage IA"))]="Stage I"
meta_for_all$stage[which(meta_for_all$stage%in%c("Stage IIA","Stage II" ,"Stage IIB","Stage IIC"))]="Stage II"
meta_for_all$stage[which(meta_for_all$stage%in%c("Stage IIIA","Stage III" ,"Stage IIIB","Stage IIIC"))]="Stage III"
meta_for_all$stage[which(meta_for_all$stage%in%c("Stage IVA","Stage IV" ,"Stage IVB"))]="Stage IV"

row_ord=c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
                "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
                "Immature Goblet","Secretory TA","fetal_crypt","fetal_Goblet cell",
                "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
                "fetal_Enterocyte","NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
                "Macrophages","CD4+ T-cells","Tregs",
                "Follicular","Inflammatory Monocytes","ILCs","Plasma","Intercept")

rel_meta=data.frame(stage=meta_for_all$stage)

hmCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
exposureScale=c(0,0.5)
hmColObj = circlize::colorRamp2(seq(exposureScale[1],exposureScale[2],length.out=length(hmCols)),hmCols)
botAnno = HeatmapAnnotation(df = rel_meta,
                            annotation_name_side = 'left')

plotExposures(fit=fit_combined_noCTA, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=meta_for_all$sample_type, row_order=row_ord)

Heatmap(fit_combined_noCTA$exposures, col=hmColObj, cluster_rows=F, show_column_names=F, 
        column_split=meta_for_all$sample_type, cluster_column_slices=F, show_column_dend=F,
        row_order = row_ord, bottom_annotation=botAnno)

master_df=as.data.frame(t(fit_combined_noCTA$exposures))
master_df$sample=rownames(master_df)

master_df$sample_type=as.character(meta_for_all$sample_type)

master_df_melted=melt(master_df, id.vars = c("sample", "sample_type"))

master_df_melted=master_df_melted[master_df_melted$sample_type%in%c("Solid Tissue Normal", "Primary Tumor", "fetal_bulk"),]



master_df$adult_epi=rowSums(master_df[,c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
  "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
  "Immature Goblet","Secretory TA")])
master_df$fetal_epi=rowSums(master_df[,c("fetal_crypt","fetal_Goblet cell",
                                         "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
                                         "fetal_Enterocyte")])
master_df$immune=rowSums(master_df[,c("NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
                                      "Macrophages","CD4+ T-cells","Tregs",
                                      "Follicular","Inflammatory Monocytes","ILCs")])
master_df$immature_enterocytes=rowSums(master_df[,c("Immature Enterocytes 1","Immature Enterocytes 2")])
master_df$TA=rowSums(master_df[,c("TA 1","TA 2", "Secretory TA")])
master_df$enterocytes=rowSums(master_df[,c("Enterocytes","Best4+ Enterocytes")])
master_df$fetal_enterocytes=rowSums(master_df[,c("fetal_BEST4 enterocyte","fetal_Enterocyte")])

master_rel_ct=master_df[,c("Stem","enterocytes", "M cells",
                           "immature_enterocytes", "TA", "immune", "fetal_enterocytes", "fetal_crypt",
                           "fetal_Goblet cell","fetal_Enteroendocrine", "Intercept", "sample_type", "sample")]
master_rel_ct_melted=melt(master_rel_ct, id.vars = c("sample", "sample_type"))
master_rel_ct_melted=master_rel_ct_melted[master_rel_ct_melted$sample_type%in%c("Solid Tissue Normal", "Primary Tumor", "fetal_bulk"),]
master_rel_ct_melted$variable=factor(master_rel_ct_melted$variable, levels=c("Intercept", "Stem", "fetal_crypt",
                                                                             "enterocytes", "M cells", "immature_enterocytes",
                                                                             "TA", "fetal_enterocytes", "fetal_Goblet cell", 
                                                                             "fetal_Enteroendocrine", "immune"))

ggplot(master_rel_ct_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values=c("grey","red","green",
                             "darkblue","darkblue","darkblue","darkblue",
                             "lightblue","lightblue","lightblue",
                             "orange"))


fetal_v_adult=master_df[,c("sample_type","adult_epi", "fetal_epi", "Intercept", "sample")]
fetal_v_adult_melted=melt(fetal_v_adult, id.vars = c("sample", "sample_type"))
fetal_v_adult_melted=fetal_v_adult_melted[fetal_v_adult_melted$sample_type%in%c("Solid Tissue Normal", "Primary Tumor", "fetal_bulk"),]
ggplot(fetal_v_adult_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +scale_fill_manual(values=c("red","blue","grey"))

fetal_score_df=fetal_v_adult_melted[fetal_v_adult_melted$variable=="fetal_epi",]

ggplot(fetal_score_df, aes(fill=variable, y=value, x=sample_type)) + 
  geom_boxplot() + labs(y="Fetal contribution", x="Sample type")


ggplot(data = fetal_score_df, aes(x=sample_type,y=value)) +
  geom_quasirandom(mapping = aes(x=sample_type,y=value), 
                         dodge.width=.8, shape=19, cex=0.5, alpha=0.4) +
  geom_boxplot(alpha=0.2, outlier.shape=NA) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))+labs(y="Fetal contribution", x="Sample type")

plotExposures(fit=fit_combined_noCTA, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=split_cols, row_order=row_ord, column_title_rot=90)
