#Tims early embryo
early_emb=readRDS('/lustre/scratch117/casm/team274/tc16/Embryo/embryo_integrated_allembryos_filtered.Rdata')

DimPlot(early_emb)

table(Idents(early_emb))

early_emb = CellCycleScoring(early_emb, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


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






early_emb@meta.data$final_annot=as.character(Idents(early_emb))


rasa_plus_early_emb=merge(rasa_srat_processed, early_emb)

scale_dat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  return(srat)
}

rasa_plus_early_emb=scale_dat(rasa_plus_early_emb)
Idents(rasa_plus_early_emb)=rasa_plus_early_emb@meta.data$final_annot
avg_exp=AverageExpression(rasa_plus_early_emb, slot="scale.data")

cordat=cor(avg_exp$RNA)
Heatmap(cordat)

#Matts cor
dat = rasa_plus_early_emb@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),rasa_plus_early_emb@meta.data$final_annot),
                           function(e) rowMeans(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)
#################################
#Gerda's cor

dim(rasa_plus_early_emb)

dim(rasa_srat)

obs=read.csv("onecs/obs.csv")
obs

length(which(rownames(early_emb)%in%rownames(rasa_srat_processed)))


##### make an embryo mtx for cellsig
rownames(obs)=obs$index

comm_genes=intersect(obs$index, rownames(early_emb))

obs2=obs[comm_genes,]

early_emb_mtx=early_emb@assays$RNA@counts[comm_genes,]

rownames(early_emb_mtx)=obs2$gene_ids

colnames(early_emb_mtx)=paste0(early_emb@meta.data$final_annot,":",colnames(early_emb_mtx))

writeMM(early_emb_mtx, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_early_emb/early_emb.mtx")
write.table(rownames(early_emb_mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_early_emb/early_emb_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(early_emb_mtx), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_early_emb/early_emb_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


rasa_mtx=rasa_srat_processed@assays$RNA@counts[comm_genes,]
rownames(rasa_mtx)=obs2$gene_ids

colnames(rasa_mtx)=paste0(rasa_srat_processed@meta.data$final_annot,":",colnames(rasa_mtx))

rasa_and_emb=cbind(early_emb_mtx, rasa_mtx)

writeMM(rasa_and_emb, "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_gut_early_emb/early_emb.mtx")
write.table(rownames(rasa_and_emb), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_gut_early_emb/early_emb_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(rasa_and_emb), "/lustre/scratch117/casm/team274/gk14/Dediff/cellsig_gut/sc_gut_early_emb/early_emb_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')

rasa_emb='cellsig_gut/out_early_emb_rasa/OutRun_fitExposures.tsv'
fit_rasa_emb=normaliseExposures(rasa_emb)


plotExposures(fit_rasa_emb, show_column_names=F)


early_emb='cellsig_gut/out_early_emb/OutRun_fitExposures.tsv'
fit_early_emb=normaliseExposures(early_emb)

mDat_col=readRDS('mDat_gut.rds')
split_cols=c(as.character(mDat_col$gdc_cases.samples.sample_type), rep("fetal_bulk", 8))
split_cols=factor(split_cols, levels = c("Solid Tissue Normal","Primary Tumor","fetal_bulk",
                                         "Recurrent Tumor","Metastatic"))
plotExposures(fit_early_emb, show_column_names=F, column_split=split_cols)


bla=fit_early_emb$exposures

bla=t(bla)

bla=as.data.frame(bla)

bla$epi_by_intercept=(bla$Epiblast-bla$Intercept)
bla$hypo_by_intercept=(bla$Hypoblast-bla$Intercept)



bla=t(bla)

fit_2=fit_early_emb
fit_2$exposures=bla

plotExposures(fit_2, show_column_names=F, column_split=split_cols)

bla=fit_early_emb$exposures
exposures= bla[1:4,]

exposures = t(t(exposures)/colSums(exposures))

fit_2$exposures=exposures


plotExposures(fit_2, show_column_names=F, column_split=split_cols)
