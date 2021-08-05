library(loomR)
gut_loom= connect("/lustre/scratch119/casm/team274sb/gk14/allgutcombined.loom", skip.validate =T)
genes=gut_loom[['row_attrs/Gene']][]
cells=gut_loom[['col_attrs/CellID']][]





unspliced_counts_per_cell = gut_loom$map(FUN = rowSums, MARGIN = 2, chunk.size = 10000, dataset.use = "layers/unspliced", 
                         display.progress = TRUE)
spliced_counts_per_cell = gut_loom$map(FUN = rowSums, MARGIN = 2, chunk.size = 10000, dataset.use = "layers/spliced", 
                                         display.progress = TRUE)

ambiguous_counts_per_cell = gut_loom$map(FUN = rowSums, MARGIN = 2, chunk.size = 10000, dataset.use = "layers/ambiguous", 
                                       display.progress = TRUE)


spliced_unspliced_df=data.frame(cells, unspliced_counts_per_cell, spliced_counts_per_cell)

spliced_unspliced_df$prop_unspliced=spliced_unspliced_df$unspliced_counts_per_cell/(spliced_unspliced_df$unspliced_counts_per_cell+spliced_unspliced_df$spliced_counts_per_cell)


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


colnames(rasa_srat)[1:10]
spliced_unspliced_df$cells[1:10]

rasa_meta=rasa_srat@meta.data

spliced_unspliced_df$sample=sapply(strsplit(spliced_unspliced_df$cells, ":"), "[[", 1)
spliced_unspliced_df$barcode=sapply(strsplit(spliced_unspliced_df$cells, ":"), "[[", 2)

spliced_unspliced_df$barcode=str_replace(spliced_unspliced_df$barcode, "x", "")

spliced_unspliced_df$sample=str_replace(spliced_unspliced_df$sample, "cellranger302_count_", "")
spliced_unspliced_df$sample=str_replace(spliced_unspliced_df$sample, "_GRCh38-3_0_0", "")
spliced_unspliced_df$sample=str_replace(spliced_unspliced_df$sample, "cellranger310_count_", "")

spliced_unspliced_df$sample[grep("^[0-9]{5}", spliced_unspliced_df$sample)]=sapply(str_split(spliced_unspliced_df$sample[grep("^[0-9]{5}", spliced_unspliced_df$sample)], "_", n=2), "[", 2)


spliced_unspliced_df$new_barcode=paste0(spliced_unspliced_df$barcode, "-1-", spliced_unspliced_df$sample)

rownames(spliced_unspliced_df)=spliced_unspliced_df$new_barcode

ratio_rel=spliced_unspliced_df[rownames(rasa_meta), ]
ratio_rel$cell_type=rasa_meta$final_annot
ratio_rel$Age_group=rasa_meta$Age_group
ratio_rel$chemistry=rasa_meta$X10X
ratio_rel$sample_name=rasa_meta$Sample.name
ratio_rel$nFeature=rasa_meta$nFeature_RNA

ratio_rel$Age_group[which(ratio_rel$Age_group%in%c("First trim","Second trim"))]="Fetal"


ratio_rel$cell_type=factor(ratio_rel$cell_type, levels = c(
  "Adult_B cells", "Fetal_B cells", 
  "Adult_T cells", "Fetal_T cells",
  "Adult_Plasma cells", "Fetal_Red blood cells", 
  "Adult_Mesenchymal", "Fetal_Mesenchymal", 
  "Adult_Myeloid", "Fetal_Myeloid", 
  "Adult_Neuronal", "Fetal_Neuronal", 
  "Adult_Enteroendocrine", "Fetal_Enteroendocrine",
  "Adult_Endothelial", "Fetal_Endothelial",
  "Adult_BEST2+ Goblet cell", "Adult_Goblet cell", "Fetal_Goblet cell", 
  "Adult_Colonocyte", "Fetal_Colonocyte", "Fetal_Enterocyte",
  "Adult_BEST4+ epithelial", "Fetal_BEST4+ epithelial", 
  "Adult_Stem cells", "Fetal_Stem cells", "Fetal_Proximal progenitor", "Fetal_Distal progenitor",
  "Adult_TA", "Fetal_TA", 
  "Adult_Tuft", "Fetal_Tuft",
  "Adult_Microfold cell", "Fetal_Microfold cell",
  "Adult_Paneth", "Fetal_Paneth", "Fetal_CLDN10+ cells"))

ratio_rel$ct2=as.character(ratio_rel$cell_type)
ratio_rel$ct2=sapply(strsplit(ratio_rel$ct2, "_"), "[", 2)

ratio_rel_3prime=ratio_rel[ratio_rel$chemistry=="3'",]
ratio_rel_5prime=ratio_rel[ratio_rel$chemistry=="5'",]


ggplot(ratio_rel, aes(x=sample_name, y=prop_unspliced, color=chemistry)) + facet_wrap(. ~ct2, ncol = 4) +
geom_boxplot() +theme(axis.text.x = element_text(angle = 90))

ggplot(ratio_rel_5prime, aes(x=cell_type, y=prop_unspliced, color=Age_group)) +
  geom_boxplot() +theme(axis.text.x = element_text(angle = 90))

ggplot(ratio_rel, aes(x=cell_type, y=prop_unspliced, color=Age_group)) +
  geom_boxplot() +theme(axis.text.x = element_text(angle = 90))


ratio_rel_3prime[which(ratio_rel_3prime$sample_name=="BRC2043"),]


ggplot(ratio_rel, aes(x=sample_name, y=nFeature, color=Age_group)) + facet_wrap(. ~ct2, ncol = 4) +
  geom_boxplot() +theme(axis.text.x = element_text(angle = 90))


ggplot(ratio_rel, aes(x=cell_type, y=nFeature, color=Age_group)) +
  geom_boxplot() +theme(axis.text.x = element_text(angle = 90))
