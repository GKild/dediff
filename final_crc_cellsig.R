####Final colorectal cellsig plots #######
source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')

tcga_cc=normaliseExposures('/home/jovyan/Dediff/cellsig_gut/out_rasa_final_tcga/OutRun_fitExposures.tsv')
tcga_nocc=normaliseExposures('/home/jovyan/Dediff/cellsig_gut/out_rasa_final_g1_tcga/OutRun_fitExposures.tsv')

boardman_cc=normaliseExposures('/home/jovyan/Dediff/cellsig_gut/out_rasa_final_boardman/OutRun_fitExposures.tsv')
boardman_nocc=normaliseExposures('/home/jovyan/Dediff/cellsig_gut/out_rasa_final_g1_boardman/OutRun_fitExposures.tsv')


plotExposures(tcga_nocc,show_column_names=F)
plotExposures(boardman_nocc,show_column_names=F)



colnames(boardman_nocc$exposures)

boardman_meta=read.table('boardman_data/boardman_meta.txt', sep = ",", header = T)
boardman_meta_new=data.frame(sample=colnames(boardman_nocc$exposures))
which(boardman_meta_new$sample%in%boardman_meta$sample)
boardman_meta_new$sample_type="x"

boardman_meta_new$sample_type=boardman_meta$sample_type[match(boardman_meta_new$sample, boardman_meta$sample)]
boardman_meta_new$sample_type[80:87]="fetal_gut_bulk"
boardman_meta_new$sample_type[which(is.na(boardman_meta_new$sample_type))]="unknown"
boardman_meta_new$sample_type_wide=boardman_meta_new$sample_type
boardman_meta_new$sample_type_wide[grep("VILLOUS", boardman_meta_new$sample_type_wide)]= "VILLOUS ADENOMA"

plotExposures(boardman_nocc,show_column_names=F, column_split=boardman_meta_new$sample_type_wide, column_title_rot=90)

master_df_boardman=as.data.frame(t(boardman_nocc$exposures))
master_df_boardman$sample=rownames(master_df_boardman)

master_df_boardman$sample_type=as.character(boardman_meta_new$sample_type_wide)


master_df_boardman_melted=melt(master_df_boardman, id.vars = c("sample", "sample_type"))

master_df_boardman_melted=master_df_boardman_melted[master_df_boardman_melted$sample_type%in%c("VILLOUS ADENOMA","CANCER","NORMAL EPITH", "fetal_gut_bulk"),]

colnames(master_df_boardman)

master_df_boardman$total_adult_signal=rowSums(master_df_boardman[,grep("Adult", colnames(master_df_boardman))])
master_df_boardman$total_fetal_signal=rowSums(master_df_boardman[,grep("Fetal", colnames(master_df_boardman))])
master_df_boardman=master_df_boardman[master_df_boardman$sample_type%in%c("VILLOUS ADENOMA","CANCER","NORMAL EPITH", "fetal_gut_bulk"),]

pdf('plots/crc_tcga_final_ref_plots/cellsig_2axis_boardman.pdf', height=8, width=8)
ggplot(master_df_boardman, aes(x=total_fetal_signal, y=Intercept, shape=sample_type, color=sample_type)) +
  geom_point()

ggplot(master_df_boardman, aes(x=total_adult_signal, y=total_fetal_signal, shape=sample_type, color=sample_type)) +
  geom_point()

ggplot(master_df_boardman, aes(x=`Fetal_Stem cells`, y=`Adult_Stem cells`, shape=sample_type, color=sample_type)) +
  geom_point()

ggplot(master_df_boardman, aes(x=Intercept, y=`Adult_Stem cells`, shape=sample_type, color=sample_type)) +
  geom_point()
dev.off()

mDat_col=readRDS('mDat_gut.rds')

tcga_meta=data.frame(sample=colnames(tcga_nocc$exposures),
                        sample_type=c(as.character(mDat_col$gdc_cases.samples.sample_type), rep("fetal_gut_bulk", 8)),
                        stage=c(mDat_col$xml_stage_event_pathologic_stage, rep(NA,8)))



plotExposures(tcga_nocc,show_column_names=F, column_split=tcga_meta$sample_type, column_title_rot=90)


master_df_tcga=as.data.frame(t(tcga_nocc$exposures))
master_df_tcga$sample=rownames(master_df_tcga)

master_df_tcga$sample_type=as.character(tcga_meta$sample_type)

master_df_tcga$total_adult_signal=rowSums(master_df_tcga[,grep("Adult", colnames(master_df_tcga))])
master_df_tcga$total_fetal_signal=rowSums(master_df_tcga[,grep("Fetal", colnames(master_df_tcga))])
master_df_tcga=master_df_tcga[master_df_tcga$sample_type%in%c("Primary Tumor", "Solid Tissue Normal", "fetal_gut_bulk"),]

pdf('plots/crc_tcga_final_ref_plots/cellsig_2axis_tcga.pdf', height=8, width=8)
ggplot(master_df_tcga, aes(x=total_fetal_signal, y=Intercept, shape=sample_type, color=sample_type)) +
  geom_point()

ggplot(master_df_tcga, aes(x=total_adult_signal, y=total_fetal_signal, shape=sample_type, color=sample_type)) +
  geom_point()

 ggplot(master_df_tcga, aes(x=`Fetal_Stem cells`, y=`Adult_Stem cells`, shape=sample_type, color=sample_type)) +
  geom_point()

ggplot(master_df_tcga, aes(x=Intercept, y=`Adult_Stem cells`, shape=sample_type, color=sample_type)) +
  geom_point()
dev.off()


df_tum=master_df_tcga[which(master_df_tcga$sample_type=="Primary Tumor"),]

df_tum=df_tum[,c(1:37)]
df_tum$best_lab=colnames(df_tum)[apply(df_tum,1,which.max)]

sort(table(df_tum$best_lab))

df_tum$best_lab=factor(df_tum$best_lab, levels=rev(names(sort(table(df_tum$best_lab)))))

hmCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
exposureScale=c(0,0.5)
hmColObj = circlize::colorRamp2(seq(exposureScale[1],exposureScale[2],length.out=length(hmCols)),hmCols)

meta_for_tum=meta_for_all[which(rownames(df_tum)%in%meta_for_all$sample),]
rownames(meta_for_tum)=meta_for_tum$sample
meta_for_tum=meta_for_tum[rownames(df_tum),]



Heatmap(t(df_tum[,c(rev(names(sort(table(df_tum$best_lab)))))]), cluster_rows = F, show_column_names = F,
        column_split = df_tum$best_lab, column_title_rot = 90, show_column_dend = F, cluster_column_slices = F)



just_tums=mDat_col

rownames(just_tums)=colnames(tcga_nocc$exposures)[1:715]
just_tums=just_tums[rownames(df_tum),]

rel_meta=data.frame(sample=rownames(just_tums), 
                    stage=just_tums$xml_stage_event_pathologic_stage, 
                    location=just_tums$tissue,
                    vital_status=just_tums$cgc_case_vital_status, 
                    disease_type=just_tums$cgc_file_disease_type,
                    percent_tum_nuc=just_tums$cgc_slide_percent_tumor_nuclei,
                    percent_tum_cell=just_tums$cgc_slide_percent_tumor_cells,
                    percent_stromal=just_tums$cgc_slide_percent_stromal_cells,
                    percent_normal=just_tums$cgc_slide_percent_normal_cells,
                    percent_necrosis=just_tums$cgc_slide_percent_necrosis,
                    pathologic_t=just_tums$cgc_case_pathologic_t,
                    pathologic_n=just_tums$cgc_case_pathologic_n,
                    gender=just_tums$xml_gender,
                    age_at_diagnosis=just_tums$xml_age_at_initial_pathologic_diagnosis, 
                    case_tumor_status=just_tums$cgc_case_tumor_status,
                    follow_up_tumour_status=just_tums$cgc_follow_up_tumor_status,
                    follow_up_vital_status=just_tums$cgc_follow_up_vital_status,
                    case_histological_diagnosis=just_tums$cgc_case_histological_diagnosis,
                    new_tum_after_treatment=just_tums$cgc_follow_up_new_tumor_event_after_initial_treatment,
                    neopl_status=just_tums$xml_person_neoplasm_cancer_status, 
                    has_new_tum_events=just_tums$xml_has_new_tumor_events_information, 
                    lymphnodes_positive=just_tums$xml_number_of_lymphnodes_positive_by_he, 
                    pretr_cea_level=just_tums$xml_preoperative_pretreatment_cea_level, 
                    venous_invasion=just_tums$xml_venous_invasion, 
                    lymphatic_invasion=just_tums$xml_lymphatic_invasion, 
                    perineural_invasion=just_tums$xml_perineural_invasion_present, 
                    kras_mutation=just_tums$xml_kras_mutation_found, 
                    braf_analysis=just_tums$xml_braf_gene_analysis_result, 
                    stem_signal=df_tum$`Adult_Stem cells`, 
                    ta_signal=df_tum$Adult_TA, 
                    intercept=df_tum$Intercept)

rel_meta$stage[which(rel_meta$stage%in%c("Stage I","Stage IA"))]="Stage I"
rel_meta$stage[which(rel_meta$stage%in%c("Stage IIA","Stage II" ,"Stage IIB","Stage IIC"))]="Stage II"
rel_meta$stage[which(rel_meta$stage%in%c("Stage IIIA","Stage III" ,"Stage IIIB","Stage IIIC"))]="Stage III"
rel_meta$stage[which(rel_meta$stage%in%c("Stage IVA","Stage IV" ,"Stage IVB"))]="Stage IV"



summary(lm(rel_meta$venous_invasion~df_tum$Adult_TA))


boxplot(rel_meta$venous_invasion, df_tum$Adult_TA)

pdf("plots/crc_tcga_final_ref_plots/high_na_ta_signal.pdf")
ggplot(rel_meta, aes(x=follow_up_tumour_status, y=ta_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=follow_up_vital_status, y=ta_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=new_tum_after_treatment, y=ta_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=perineural_invasion, y=ta_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=kras_mutation, y=ta_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=braf_analysis, y=ta_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=pretr_cea_level, y=ta_signal)) +geom_point()
dev.off()

pdf("plots/crc_tcga_final_ref_plots/high_na_stem_signal.pdf")
ggplot(rel_meta, aes(x=follow_up_tumour_status, y=stem_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=follow_up_vital_status, y=stem_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=new_tum_after_treatment, y=stem_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=perineural_invasion, y=stem_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=kras_mutation, y=stem_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=braf_analysis, y=stem_signal)) +geom_boxplot()
ggplot(rel_meta, aes(x=pretr_cea_level, y=stem_signal)) +geom_point()
dev.off()

high_na_cols=c("follow_up_tumour_status", "follow_up_vital_status", "new_tum_after_treatment",
               "pretr_cea_level", "perineural_invasion", "kras_mutation",  "braf_analysis")


rel_meta_low_nas=rel_meta[,colnames(rel_meta)[which(!colnames(rel_meta)%in%high_na_cols)]]

cat_cols=c("stage", "location", "vital_status", "disease_type",
           "pathologic_t", "pathologic_n", "gender",
           "case_tumor_status", "case_histological_diagnosis",
           "neopl_status", "has_new_tum_events", 
           "venous_invasion", "lymphatic_invasion")

rel_meta_low_nas[,cat_cols] <- lapply(rel_meta_low_nas[,cat_cols], factor)

sum(complete.cases(rel_meta_low_nas))

lm1=lm(ta_signal ~ stage + location + vital_status +disease_type + percent_tum_nuc + percent_tum_cell + percent_stromal +
         percent_normal + percent_necrosis +pathologic_t +pathologic_n + gender + age_at_diagnosis +case_tumor_status +
         case_histological_diagnosis + neopl_status + has_new_tum_events + lymphnodes_positive + venous_invasion +lymphatic_invasion, data = rel_meta_low_nas)

summary(lm1)

lm2=lm(stem_signal ~ stage + location + vital_status +disease_type + percent_tum_nuc + percent_tum_cell + percent_stromal +
             percent_normal + percent_necrosis +pathologic_t +pathologic_n + gender + age_at_diagnosis +case_tumor_status +
             case_histological_diagnosis + neopl_status + has_new_tum_events + lymphnodes_positive + venous_invasion +lymphatic_invasion, data = rel_meta_low_nas)

summary(lm2)


lm3=lm(intercept ~ stage + location + vital_status +disease_type + percent_tum_nuc + percent_tum_cell + percent_stromal +
         percent_normal + percent_necrosis +pathologic_t +pathologic_n + gender + age_at_diagnosis +case_tumor_status +
         case_histological_diagnosis + neopl_status + has_new_tum_events + lymphnodes_positive + venous_invasion +lymphatic_invasion, data = rel_meta_low_nas)

summary(lm3)


