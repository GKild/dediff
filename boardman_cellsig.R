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
