r_mtx=readMM('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/matrix.mtx')
r_genes=read.csv('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/obs.csv')
r_meta=read.csv('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/var.csv')
rownames(r_meta)=r_meta$X

rownames(r_mtx)=r_genes$gene_ids
colnames(r_mtx)=r_meta$X


which(r_meta$Age_group%in%c("Adult", "First trim", "Second trim"))
which(r_meta$Region%in%c("APD", "LargeInt", "REC"))

r_meta_age=r_meta[which(r_meta$Age_group%in%c("Adult", "First trim", "Second trim")),]
r_meta_age_region=r_meta_age[which(r_meta_age$Region%in%c("APD", "LargeInt", "REC")),]


r_rel_mtx=r_mtx[,r_meta_age_region$X]


rasa_srat=CreateSeuratObject(r_rel_mtx, meta.data = r_meta_age_region)

saveRDS(rasa_srat, '/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/rasa_rel_srat.rds')

rasa_srat=readRDS('/lustre/scratch117/casm/team274/gk14/rasa_final_object_soupx/rasa_rel_srat.rds')

rasa_srat@meta.data$annotation_narrow=rasa_srat@meta.data$Integrated_05

rasa_srat@meta.data$annotation_narrow[which(rasa_srat@meta.data$annotation_narrow%in%c("D cells (SST+)","EC cells (NPW+)", "EC cells (TAC1+)",
                                                                           "EECs","I cells (CCK+)", "L cells (PYY+)","M/X cells (MLN/GHRL+)","N cells (NTS+)", "Progenitor (NEUROG3+)"))]="Enteroendocrine"


rasa_srat@meta.data$annotation_wide=rasa_srat@meta.data$category
rasa_srat@meta.data$fetal_or_adult=rasa_srat@meta.data$Diagnosis
rasa_srat@meta.data$fetal_or_adult[which(rasa_srat@meta.data$fetal_or_adult=="fetal")]="Fetal_"
rasa_srat@meta.data$fetal_or_adult[which(rasa_srat@meta.data$fetal_or_adult=="Healthy adult")]="Adult_"
rasa_srat@meta.data$annotation_wide=paste0(rasa_srat@meta.data$fetal_or_adult, rasa_srat@meta.data$annotation_wide)

rasa_srat@meta.data$annotation_narrow=paste0(rasa_srat@meta.data$fetal_or_adult,rasa_srat@meta.data$annotation_narrow)


rasa_srat@meta.data$final_annot=rasa_srat@meta.data$annotation_wide
rasa_srat@meta.data$final_annot[which(rasa_srat@meta.data$final_annot%in%c("Adult_Epithelial", "Fetal_Epithelial"))]=rasa_srat@meta.data$annotation_narrow[which(rasa_srat@meta.data$final_annot%in%c("Adult_Epithelial", "Fetal_Epithelial"))]
 
rasa_srat_processed=process_seurat(rasa_srat)


#Matts cor
dat = rasa_srat_processed@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),rasa_srat_processed@meta.data$final_annot),
                           function(e) rowMeans(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)
#################################
#Gerda's cor
Idents(rasa_srat_processed)=rasa_srat_processed@meta.data$final_annot
avg_exp=AverageExpression(rasa_srat_processed, slot="scale.data")

cordat=cor(avg_exp$RNA)
Heatmap(cordat)
#################################


rasa_srat_processed=subset(rasa_srat_processed, Integrated_05%in%c("CLDN10+ cells"), invert=T)
DimPlot(rasa_srat_processed, group.by = "Phase")



rasa_srat_processed_g1=subset(rasa_srat_processed, Phase=="G1")

obs=read.csv("onecs/obs.csv")
rasa_mtx=rasa_srat_processed@assays$RNA@counts
rownames(rasa_mtx)=obs$gene_ids

colnames(rasa_mtx)=paste0(rasa_srat_processed@meta.data$final_annot,":",colnames(rasa_mtx))

writeMM(rasa_mtx, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_rasa_final/rasa.mtx")
write.table(rownames(rasa_mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_rasa_final/rasa_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(rasa_mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_rasa_final/rasa_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


rasa_mtx_g1=rasa_srat_processed_g1@assays$RNA@counts
rownames(rasa_mtx_g1)=obs$gene_ids

colnames(rasa_mtx_g1)=paste0(rasa_srat_processed_g1@meta.data$final_annot,":",colnames(rasa_mtx_g1))

writeMM(rasa_mtx_g1, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_rasa_final_g1/rasag1.mtx")
write.table(rownames(rasa_mtx_g1), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_rasa_final_g1/rasag1_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(rasa_mtx_g1), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_rasa_final_g1/rasag1_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)






