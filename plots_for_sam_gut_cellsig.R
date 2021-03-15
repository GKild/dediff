library(dplyr)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)
source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')
mDat_col=readRDS('tcga_meta.rds')
split_cols=c(as.character(mDat_col$gdc_cases.samples.sample_type), rep("fetal_bulk", 8))
split_cols=factor(split_cols, levels = c("Solid Tissue Normal","Primary Tumor","fetal_bulk",
                                         "Recurrent Tumor","Metastatic"))
combined_noCTA='cellsig_gut/out_combined/OutRun_fitExposures.tsv'
fit_combined_noCTA=normaliseExposures(combined_noCTA)

row_ord=c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
          "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
          "Immature Goblet","Secretory TA","fetal_crypt","fetal_Goblet cell",
          "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
          "fetal_Enterocyte","NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
          "Macrophages","CD4+ T-cells","Tregs",
          "Follicular","Inflammatory Monocytes","ILCs","Plasma","Intercept")

plotExposures(fit=fit_combined_noCTA, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=split_cols, row_order=row_ord)

master_df=as.data.frame(t(fit_combined_noCTA$exposures))
master_df$sample=rownames(master_df)

master_df$sample_type=as.character(split_cols)

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

master_df$other_adult_epi=rowSums(master_df[,c("enterocytes","M cells","immature_enterocytes",
                                                                 "TA")])
master_df$other_fetal_epi=rowSums(master_df[,c("fetal_enterocytes","fetal_Goblet cell",
                                                                 "fetal_Enteroendocrine")])



master_rel_ct=master_df[,c("Stem","other_adult_epi", "immune", "fetal_crypt",
                           "other_fetal_epi", "Intercept", "sample_type", "sample")]
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
fetal_v_adult$fetal_to_adult=fetal_v_adult$fetal_epi/(fetal_v_adult$adult_epi+fetal_v_adult$fetal_epi)
fetal_v_adult_melted=melt(fetal_v_adult, id.vars = c("sample", "sample_type"))
fetal_v_adult_melted=fetal_v_adult_melted[fetal_v_adult_melted$sample_type%in%c("Solid Tissue Normal", "Primary Tumor", "fetal_bulk"),]
ggplot(fetal_v_adult_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +scale_fill_manual(values=c("red","blue","grey"))

fetal_score_df=fetal_v_adult_melted[fetal_v_adult_melted$variable=="fetal_to_adult",]
fetal_score_df$sample_type[which(fetal_score_df$sample_type=="fetal_bulk")]="Fetal Normal bulk"
fetal_score_df$sample_type[which(fetal_score_df$sample_type=="Solid Tissue Normal")]="Adult Normal bulk"
fetal_score_df$sample_type=factor(fetal_score_df$sample_type, levels=c("Fetal Normal bulk", "Adult Normal bulk", "Primary Tumor"))
ggplot(fetal_score_df, aes(fill=variable, y=value, x=sample_type)) + 
  geom_boxplot() + labs(y="Fetal contribution", x="Sample type")


ggplot(data = fetal_score_df, aes(x=sample_type,y=value)) +
  geom_boxplot(fill=NA,outlier.shape=NA, colour="grey18") + theme_bw() +
  geom_quasirandom(mapping = aes(x=sample_type,y=value), 
                   dodge.width=.8, shape=19, cex=0.5, alpha=0.4) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))+labs(y="Fetal contribution", x="Sample type")

wilcox.test(fetal_score_df$value[which(fetal_score_df$sample_type=="Solid Tissue Normal")],
            fetal_score_df$value[which(fetal_score_df$sample_type=="Primary Tumor"
                                       t.test(fetal_score_df$value[which(fetal_score_df$sample_type=="Solid Tissue Normal")],
                                              fetal_score_df$value[which(fetal_score_df$sample_type=="Primary Tumor")])
                                       
                                       
                                       fetal_cont_plot=function()
                                         