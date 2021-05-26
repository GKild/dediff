library(Matrix)
library(Seurat)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes


#Fetal single cell from Alison Simmons' group
simmons_fetal=readRDS('simmons_fetal_data/epithelium.RDS')

process_seurat =function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

simmons_fetal_processed=process_seurat(simmons_fetal)


DimPlot(simmons_fetal_processed, group.by = "Location")
simmons_fetal_colon=subset(simmons_fetal_processed, BroadLocation=="Colon")

DimPlot(simmons_fetal_colon, group.by = "old.ident")

DimPlot(simmons_fetal_colon,cells.highlight = rownames(simmons_fetal_processed@meta.data)[which(simmons_fetal_processed@meta.data$old.ident%in%c("Distal Stem Cells"))])


simmons_fetal_colon@meta.data$old.ident


table(simmons_fetal_colon@meta.data[which(simmons_fetal_colon@meta.data$old.ident%in%c("Distal Stem Cells")),]$Phase)

#Fetal and adult gut epithelium data from Rasa - updated!
r_a_mtx=readMM('rasa_stuff/rasa_adult_updated/matrix.mtx')
r_a_genes=read.csv('rasa_stuff/rasa_adult_updated/obs.csv')
r_a_meta=read.csv('rasa_stuff/rasa_adult_updated/var.csv')
rownames(r_a_meta)=r_a_meta$X

rownames(r_a_mtx)=r_a_genes$gene_ids
colnames(r_a_mtx)=r_a_meta$X

r_f_mtx=readMM('rasa_stuff/rasa_fetal_updated/matrix.mtx')
r_f_genes=read.csv('rasa_stuff/rasa_fetal_updated/obs.csv')
r_f_meta=read.csv('rasa_stuff/rasa_fetal_updated/var.csv')
rownames(r_f_meta)=r_f_meta$X

rownames(r_f_mtx)=r_f_genes$gene_ids
colnames(r_f_mtx)=r_f_meta$X

#create seurat objects from dGc matrices
r_a_srat=CreateSeuratObject(r_a_mtx, meta.data = r_a_meta)
r_f_srat=CreateSeuratObject(r_f_mtx, meta.data = r_f_meta)

#only want healthy adults + large intestine from adult data
r_a_srat=subset(r_a_srat, Diagnosis=="Healthy adult")
r_a_srat=subset(r_a_srat, Region%in%c("LargeInt", "REC"))
# only want large int from fetal data
r_f_srat=subset(r_f_srat, Region%in%c("LargeInt"))

# update fetal metadata with Rasa's column 

check1=read.csv('rasa_stuff/rasa_fetal_updated/fetal_epi_anno (1).csv')
rownames(check1)=check1$X
check1=check1[colnames(r_f_srat),]

r_f_srat@meta.data$annotation_2=check1$annotation_2
r_a_srat@meta.data$annotation_2=r_a_srat@meta.data$annotation
r_a_srat@meta.data$annotation_2[which(r_a_srat@meta.data$annotation_2%in%c("D cells (SST+)", "EC cells (TAC1+)",
                                                                       "EECs", "L cells (PYY+)","N cells (NTS+)", "Progenitor (NEUROG3+)"))]="Enteroendocrine"

r_a_srat@meta.data$annotation_2=paste0("Adult_", r_a_srat@meta.data$annotation_2)
r_f_srat@meta.data$annotation_2=paste0("Fetal_", r_f_srat@meta.data$annotation_2)


saveRDS(r_f_srat, "rasa_stuff/rasa_fetal_updated/filtered_fetal_updated.rds")
saveRDS(r_a_srat, "rasa_stuff/rasa_adult_updated/filtered_adult_updated.rds")

r_f_srat=readRDS("rasa_stuff/rasa_fetal_updated/filtered_fetal_updated.rds")
r_a_srat=readRDS("rasa_stuff/rasa_adult_updated/filtered_adult_updated.rds")


rasa_updated_srat=merge(r_a_srat, r_f_srat)

#add adult immune cells into the reference --pending, waiting on data from Rasa


rasa_updated_srat_processed=process_seurat(rasa_updated_srat)

DimPlot(rasa_updated_srat_processed, group.by = "annotation_2")

DimPlot(rasa_updated_srat_processed, cells.highlight = rownames(rasa_updated_srat_processed@meta.data)[which(rasa_updated_srat_processed@meta.data$annotation=="Fetal_Stem cells")])
DimPlot(rasa_updated_srat_processed, cells.highlight = rownames(rasa_updated_srat_processed@meta.data)[which(rasa_updated_srat_processed@meta.data$annotation=="Adult_Stem cells")])

DimPlot(rasa_updated_srat_processed, group.by = "batch") + theme(legend.position="none")

Idents(rasa_updated_srat_processed)=rasa_updated_srat_processed@meta.data$annotation

DotPlot(rasa_updated_srat_processed, features = c("LGR5", "OLFM4", "SI", "ANPEP", "CLCA1", "ZG16", "BEST4", "OTOP2"))
DimPlot(rasa_updated_srat_processed, group.by = "Phase")

rasa_updated_srat_processed_g1=subset(rasa_updated_srat_processed, Phase=="G1")

DotPlot(rasa_updated_srat_processed_g1, features = c("LGR5", "OLFM4", "SI", "ANPEP", "CLCA1", "ZG16", "BEST4", "OTOP2"))

obs=read.csv("onecs/obs.csv")
rasa_mtx=rasa_updated_srat_processed@assays$RNA@counts
rownames(rasa_mtx)=obs$gene_ids

colnames(rasa_mtx)=paste0(rasa_updated_srat_processed@meta.data$annotation,":",colnames(rasa_mtx))

writeMM(rasa_mtx, "cellsig_gut/sc_rasa_updated/rasa.mtx")
write.table(rownames(rasa_mtx), "cellsig_gut/sc_rasa_updated/rasa_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(rasa_mtx), "cellsig_gut/sc_rasa_updated/rasa_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)


rasa_mtx_g1=rasa_updated_srat_processed_g1@assays$RNA@counts
rownames(rasa_mtx_g1)=obs$gene_ids

colnames(rasa_mtx_g1)=paste0(rasa_updated_srat_processed_g1@meta.data$annotation,":",colnames(rasa_mtx_g1))

writeMM(rasa_mtx_g1, "cellsig_gut/sc_rasa_updated/rasag1.mtx")
write.table(rownames(rasa_mtx_g1), "cellsig_gut/sc_rasa_updated/rasag1_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(rasa_mtx_g1), "cellsig_gut/sc_rasa_updated/rasag1_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
