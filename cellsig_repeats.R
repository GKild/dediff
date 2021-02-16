test1='cellsig_gut/out_test1/OutRun_fitExposures.tsv'
fit_test1=normaliseExposures(test1)

test2='cellsig_gut/out_test2/OutRun_fitExposures.tsv'
fit_test2=normaliseExposures(test2)

test3='cellsig_gut/out_test3/OutRun_fitExposures.tsv'
fit_test3=normaliseExposures(test3)

row_ord=c("Stem","Enterocyte Progenitors", "Enteroendocrine","Tuft","Goblet","TA 1","TA 2","M cells",
          "Enterocytes","Best4+ Enterocytes","Immature Enterocytes 1","Immature Enterocytes 2",
          "Immature Goblet","Secretory TA","fetal_crypt","fetal_Goblet cell",
          "fetal_BEST4 enterocyte","fetal_Enteroendocrine",
          "fetal_Enterocyte","NKs","MT-hi","GC","CD8+ T-cells","Mast cells","DCs",
          "Macrophages","CD4+ T-cells","Tregs",
          "Follicular","Inflammatory Monocytes","ILCs","Plasma","Intercept")


pdf("cellsig_repeats.pdf", height = 10, width=10)
print(plotExposures(fit=fit_combined_noCTA, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=meta_for_all$sample_type, row_order=row_ord))


print(plotExposures(fit=fit_test1, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=meta_for_all$sample_type, row_order=row_ord))

print(plotExposures(fit=fit_test2, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=meta_for_all$sample_type, row_order=row_ord))

print(plotExposures(fit=fit_test3, show_column_names=F, cluster_column_slices=F, show_column_dend=F,
              column_split=meta_for_all$sample_type, row_order=row_ord))

dev.off()


# non-stem tumours

df_tum=master_df[which(master_df$sample_type=="Primary Tumor"),]

df_tum=df_tum[,c(1:33)]
df_tum$best_lab=colnames(df_tum)[apply(df_tum,1,which.max)]

sort(table(df_tum$best_lab))

df_tum$best_lab=factor(df_tum$best_lab, levels=rev(names(sort(table(df_tum$best_lab)))))

hmCols = c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
exposureScale=c(0,0.5)
hmColObj = circlize::colorRamp2(seq(exposureScale[1],exposureScale[2],length.out=length(hmCols)),hmCols)

meta_for_tum=meta_for_all[which(rownames(df_tum)%in%meta_for_all$sample),]
rownames(meta_for_tum)=meta_for_tum$sample
meta_for_tum=meta_for_tum[rownames(df_tum),]

rel_meta=data.frame(stage=meta_for_tum$stage)
botAnno = HeatmapAnnotation(df = rel_meta,
                            annotation_name_side = 'left')

Heatmap(t(df_tum[,c(rev(names(sort(table(df_tum$best_lab)))))]), cluster_rows = F, show_column_names = F,
        column_split = df_tum$best_lab, col=hmColObj, show_column_dend = F, cluster_column_slices = F,
        column_gap = unit(c(2,rep(5,11)), "mm"), row_order= levels(df_tum$best_lab), column_title_rot = 90,
        width = unit(6 ,"in"), height = unit(4, "in"), bottom_annotation=botAnno)
