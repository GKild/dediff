boardman_meta=read.table('boardman_meta.txt', sep = ",", header = T)
boardman_data=read.delim('boardman_data.tsv', sep = "\t", header = T)

boardman_rel=boardman_data[,c(4,5,grep("count", colnames(boardman_data)))]

colnames(boardman_rel)=gsub("_count", "", colnames(boardman_rel))

length(which(boardman_meta$sample%in%colnames(boardman_rel)))

colnames(boardman_rel)=gsub("Unk.", "", colnames(boardman_rel))

boardman_meta_new=data.frame(sample=colnames(boardman_rel)[3:81])
which(boardman_meta_new$sample%in%boardman_meta$sample)
boardman_meta_new$sample_type="x"

boardman_meta_new$sample_type=boardman_meta$sample_type[match(boardman_meta_new$sample, boardman_meta$sample)]

boardman_meta_new$sample_type[which(is.na(boardman_meta_new$sample_type))]="unknown"

colnames(boardman_rel)[1]="geneLengths"
colnames(boardman_rel)[2]="geneName"
  
for (x in colnames(boardman_rel)[3:81]) {
    write.table(boardman_rel[1:64253,c("geneName","geneLengths",x)], paste0('cellsig_gut/boardman_data/',x,'.tsv'), sep="\t", row.names = F, col.names = T, quote = F)
  }

boardman_paths=paste0('/home/jovyan/Dediff/cellsig_gut/boardman_data/', colnames(boardman_rel)[3:81], '.tsv')

write.table(boardman_paths, "/home/jovyan/Dediff/cellsig_gut/boardman_paths.txt",sep = "\t", row.names = F, col.names = F, quote=F)


source('/home/jovyan/Dediff/matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')

baordman='cellsig_gut/out_boardman/OutRun_fitExposures.tsv'
fit_boardman=normaliseExposures(baordman)
row_ord=c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
          "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
          "Immature Goblet","Secretory TA","fetal_crypt","fetal_Goblet cell",
          "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
          "fetal_Enterocyte","NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
          "Macrophages","CD4+ T-cells","Tregs",
          "Follicular","Inflammatory Monocytes","ILCs","Plasma","Intercept")


plotExposures(fit_boardman, show_column_names=F, column_split=boardman_meta_new$sample_type, row_order=row_ord, column_title_rot=90, cluster_column_slices=F)


master_df_boardman=as.data.frame(t(fit_boardman$exposures))
master_df_boardman$sample=rownames(master_df_boardman)

master_df_boardman$sample_type=as.character(boardman_meta_new$sample_type)

master_df_boardman$sample_type[grep("VILLOUS", master_df_boardman$sample_type)]="VILLOUS ADENOMA"

master_df_boardman_melted=melt(master_df_boardman, id.vars = c("sample", "sample_type"))

master_df_boardman_melted=master_df_boardman_melted[master_df_boardman_melted$sample_type%in%c("VILLOUS ADENOMA","CANCER","NORMAL EPITH"),]



master_df_boardman$adult_epi=rowSums(master_df_boardman[,c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
                                                           "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
                                                           "Immature Goblet","Secretory TA")])
master_df_boardman$fetal_epi=rowSums(master_df_boardman[,c("fetal_crypt","fetal_Goblet cell",
                                                           "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
                                                           "fetal_Enterocyte")])
master_df_boardman$immune=rowSums(master_df_boardman[,c("NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
                                                        "Macrophages","CD4+ T-cells","Tregs",
                                                        "Follicular","Inflammatory Monocytes","ILCs")])
master_df_boardman$immature_enterocytes=rowSums(master_df_boardman[,c("Immature Enterocytes 1","Immature Enterocytes 2")])
master_df_boardman$TA=rowSums(master_df_boardman[,c("TA 1","TA 2", "Secretory TA")])
master_df_boardman$enterocytes=rowSums(master_df_boardman[,c("Enterocytes","Best4+ Enterocytes")])
master_df_boardman$fetal_enterocytes=rowSums(master_df_boardman[,c("fetal_BEST4 enterocyte","fetal_Enterocyte")])

boardman_master_rel_ct=master_df_boardman[,c("Stem","enterocytes", "M cells",
                                             "immature_enterocytes", "TA", "immune", "fetal_enterocytes", "fetal_crypt",
                                             "fetal_Goblet cell","fetal_Enteroendocrine", "Intercept", "sample_type", "sample")]
boardman_master_rel_ct_melted=melt(boardman_master_rel_ct, id.vars = c("sample", "sample_type"))
boardman_master_rel_ct_melted=boardman_master_rel_ct_melted[boardman_master_rel_ct_melted$sample_type%in%c(c("VILLOUS ADENOMA","CANCER","NORMAL EPITH")),]
boardman_master_rel_ct_melted$sample=factor(boardman_master_rel_ct_melted$sample,
                                            levels = master_df_boardman$sample[order(master_df_boardman$Stem, decreasing = F)])
master_df_boardman$sample[order(master_df_boardman$Stem, decreasing = F)]


boardman_master_rel_ct_melted$variable=factor(boardman_master_rel_ct_melted$variable, levels=c("Intercept",
                                                                                               "enterocytes", "M cells", "immature_enterocytes",
                                                                                               "TA", "fetal_enterocytes", "fetal_Goblet cell", 
                                                                                               "fetal_Enteroendocrine", "immune", "fetal_crypt", "Stem"))

ggplot(boardman_master_rel_ct_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_fill_manual(values=c("grey",
                             "darkblue","darkblue","darkblue","darkblue",
                             "lightblue","lightblue","lightblue",
                             "orange", "green", "red"))
ggplot(boardman_master_rel_ct_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  #strip.background = element_blank(),
  #panel.border = element_rect(colour = "black")) +
  scale_fill_manual(values=c("grey","red","green",
                             "darkblue","darkblue","darkblue","darkblue",
                             "lightblue","lightblue","lightblue",
                             "orange"))


boardman_fetal_v_adult=master_df_boardman[,c("sample_type","adult_epi", "fetal_epi", "Intercept", "sample")]
boardman_fetal_v_adult_melted=melt(boardman_fetal_v_adult, id.vars = c("sample", "sample_type"))
boardman_fetal_v_adult_melted=boardman_fetal_v_adult_melted[boardman_fetal_v_adult_melted$sample_type%in%c(c("VILLOUS ADENOMA","CANCER","NORMAL EPITH")),]
ggplot(boardman_fetal_v_adult_melted, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity") + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +scale_fill_manual(values=c("red","blue","grey"))

boardman_fetal_score_df=boardman_fetal_v_adult_melted[boardman_fetal_v_adult_melted$variable=="fetal_epi",]
boardman_fetal_score_df$sample_type=factor(boardman_fetal_score_df$sample_type, levels=c("NORMAL EPITH", "VILLOUS ADENOMA", "CANCER"))

ggplot(boardman_fetal_score_df, aes(fill=variable, y=value, x=sample_type)) + 
  geom_boxplot() + labs(y="Fetal contribution", x="Sample type")


ggplot(data = boardman_fetal_score_df, aes(x=sample_type,y=value)) +
  geom_quasirandom(mapping = aes(x=sample_type,y=value), 
                   dodge.width=.8, shape=19, cex=0.5, alpha=0.4) +
  geom_boxplot(alpha=0.1, outlier.shape=NA) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))+labs(y="Fetal contribution", x="Sample type")






boardman_master_rel_ct_melted$sample=factor(boardman_master_rel_ct_melted$sample,
                                            levels = master_df_boardman$sample[order(master_df_boardman$Stem, decreasing = F)])
master_df_boardman$sample[order(master_df_boardman$Stem, decreasing = F)]

board_cellcomp=boardman_master_rel_ct_melted
board_cellcomp$variable=as.character(board_cellcomp$variable)
board_cellcomp$variable[which(board_cellcomp$variable=="fetal_crypt")]="Fetal Stem"
board_cellcomp$variable[which(board_cellcomp$variable=="Stem")]="Adult Stem"
board_cellcomp$variable[which(board_cellcomp$variable%in%c("enterocytes","M cells","immature_enterocytes",
                                                           "TA"))]="Other Adult Epithelium"
board_cellcomp$variable[which(board_cellcomp$variable%in%c("fetal_enterocytes","fetal_Goblet cell",
                                                           "fetal_Enteroendocrine"))]="Other Fetal Epithelium"
board_cellcomp$variable[which(board_cellcomp$variable=="immune")]="Adult Immune"

board_cellcomp$variable=factor(board_cellcomp$variable, levels=c("Intercept", "Adult Immune",
                                                                 "Other Adult Epithelium", 
                                                                 "Other Fetal Epithelium",
                                                                 "Fetal Stem", 
                                                                 "Adult Stem"))
board_cellcomp$sample_type=factor(board_cellcomp$sample_type, levels=c("NORMAL EPITH", "VILLOUS ADENOMA", "CANCER"))

ggplot(board_cellcomp, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity", width = 1) + facet_grid(. ~ sample_type, scales="free_x", space="free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_fill_manual(values=c("grey", "#DDCC77",
                             "#4477AA",
                             "#88CCEE",
                             "#332288", "#AA4499"))

tcga_cellcomp=master_rel_ct_melted

tcga_cellcomp$sample=factor(tcga_cellcomp$sample, 
                            levels = master_df$sample[order(master_df$Stem, decreasing = F)])

tcga_cellcomp$variable=as.character(tcga_cellcomp$variable)
tcga_cellcomp$variable[which(tcga_cellcomp$variable=="fetal_crypt")]="Fetal Stem"
tcga_cellcomp$variable[which(tcga_cellcomp$variable=="Stem")]="Adult Stem"
tcga_cellcomp$variable[which(tcga_cellcomp$variable%in%c("enterocytes","M cells","immature_enterocytes",
                                                         "TA"))]="Other Adult Epithelium"
tcga_cellcomp$variable[which(tcga_cellcomp$variable%in%c("fetal_enterocytes","fetal_Goblet cell",
                                                         "fetal_Enteroendocrine"))]="Other Fetal Epithelium"
tcga_cellcomp$variable[which(tcga_cellcomp$variable=="immune")]="Adult Immune"

tcga_cellcomp$variable=factor(tcga_cellcomp$variable, levels=c("Intercept", "Adult Immune",
                                                               "Other Adult Epithelium", 
                                                               "Other Fetal Epithelium",
                                                               "Fetal Stem", 
                                                               "Adult Stem"))

g= ggplot(tcga_cellcomp, aes(fill=variable, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity", width = 1) + facet_grid(. ~ sample_type, scales="free_x", space="free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) +
  scale_fill_manual(values=c("grey", "#DDCC77",
                             "#4477AA",
                             "#88CCEE",
                             "#332288", "#AA4499"))
library(grid)
gt = ggplot_gtable(ggplot_build(g))
gt$widths[5] = 9*gt$widths[5]
gt$widths[9] = 4*gt$widths[9]
grid.draw(gt)