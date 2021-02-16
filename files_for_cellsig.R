library(Matrix)
library(Seurat)
#Regev data - epithelial and immune matrices separately
epi_mtx=readMM('regev_adult_gut/gene_sorted-Epi.matrix.mtx')
epi_barcodes=read.table('regev_adult_gut/Epi.barcodes2.tsv')
epi_genes=read.table('regev_adult_gut/Epi.genes.tsv')

imm_mtx=readMM('regev_adult_gut/gene_sorted-Imm.matrix.mtx')
imm_barcodes=read.table('regev_adult_gut/Imm.barcodes2.tsv')
imm_genes=read.table('regev_adult_gut/Imm.genes.tsv')


meta=read.table('regev_adult_gut/all.meta2.txt', sep = "\t", header = T)
rownames(meta)=meta$NAME
rownames(epi_mtx)=epi_genes$V1
colnames(epi_mtx)=epi_barcodes$V1
epi_meta=meta[colnames(epi_mtx), ]
healthy_epi_meta=epi_meta[epi_meta$Health=="Healthy",]
healthy_epi_mtx=epi_mtx[,rownames(healthy_epi_meta)]
epi_seurat=CreateSeuratObject(healthy_epi_mtx, meta.data = healthy_epi_meta)

rownames(imm_mtx)=imm_genes$V1
colnames(imm_mtx)=imm_barcodes$V1
imm_meta=meta[colnames(imm_mtx), ]
healthy_imm_meta=imm_meta[imm_meta$Health=="Healthy",]
healthy_imm_mtx=imm_mtx[,rownames(healthy_imm_meta)]
imm_seurat=CreateSeuratObject(healthy_imm_mtx, meta.data = healthy_imm_meta)

srat_comb=merge(epi_seurat,imm_seurat)

s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes
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

srat_comb=process_seurat(srat_comb)

ct=rownames(srat_comb@meta.data)[which(!srat_comb@meta.data$Cluster%in%c("Cycling TA", "Cycling T", "Cycling B","Cycling Monocytes"))]
srat_comb_nocc=subset(srat_comb, cells=ct)
srat_comb_nocc=subset(srat_comb_nocc, Phase=="G1")

#combine some CTs into groups because we don't care about this level of detail

srat_comb_nocc@meta.data$new_ident=srat_comb_nocc@meta.data$Cluster

srat_comb_nocc@meta.data$new_ident[grep("CD4", srat_comb_nocc@meta.data$new_ident)]="CD4+ T-cells"
srat_comb_nocc@meta.data$new_ident[grep("CD8", srat_comb_nocc@meta.data$new_ident)]="CD8+ T-cells"
srat_comb_nocc@meta.data$new_ident[grep("Mast", srat_comb_nocc@meta.data$new_ident)]="Mast cells"
srat_comb_nocc@meta.data$new_ident[grep("DC", srat_comb_nocc@meta.data$new_ident)]="DCs"


# only epithelial data for Rasa's large intestine fetal data

li_seurat=readRDS('large_intestine_rasa.rds')
li_seurat_fetal=subset(li_seurat, Age%in%c("F78", "F66", "F67", "F72", "F73"))
lif_nocc_seurat=subset(li_seurat_fetal, Phase=="G1")

fetal_epi=subset(lif_nocc_seurat, exp_annot%in%c("fetal_Goblet cell","fetal_Enterocyte",
                                                "fetal_BEST4 enterocyte","fetal_crypt",
                                                "fetal_Enteroendocrine"))
# prep matrices for combination, combine, save in the format needed for cell signal analysis

obs=read.csv("onecs/obs.csv")

comm_g=intersect(rownames(srat_comb_nocc), obs$index)
obs2=obs
rownames(obs2)=obs2$index
obs2=obs2[comm_g,]

adult_mtx=srat_comb_nocc@assays$RNA@counts
adult_mtx=adult_mtx[rownames(obs2),]
rownames(adult_mtx)=obs2$gene_ids
colnames(adult_mtx)=paste0(srat_comb_nocc@meta.data$new_ident,":",colnames(adult_mtx))

fetal_mtx=fetal_epi@assays$RNA@counts
rownames(fetal_mtx)=obs$gene_ids
colnames(fetal_mtx)=paste0(fetal_epi@meta.data$exp_annot,":",colnames(fetal_mtx))

comb_epi_rasa_rn=intersect(rownames(fetal_mtx), rownames(adult_mtx))

combined_mtx=cbind(fetal_mtx[comb_epi_rasa_rn,], adult_mtx[comb_epi_rasa_rn,])
combined_mtx=as.matrix(combined_mtx)
combined_mtx=Matrix(combined_mtx, sparse = TRUE) 

writeMM(combined_mtx, "cellsig_gut/sc_combined/combined.mtx")
write.table(rownames(combined_mtx), "cellsig_gut/sc_combined/combined_rowNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(colnames(combined_mtx), "cellsig_gut/sc_combined/combined_columnNames.tsv", sep = "\t", quote = F, col.names = F, row.names = F)

#sort out fetal bulk into format for cell signal analysis

load('Foetal_RawCounts_Zilbauer.RData')

zilb_counts=counts
zilb_meta=foetal.info

zilb_meta$sample.id[which(zilb_meta$sample.type=="purified")]

zilb_fetal_counts=zilb_counts[,zilb_meta$sample.id[which(zilb_meta$sample.type=="purified")]]

exmp=read.table(
  'matts-bits-of-code-master/cellSignalAnalysis/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv',
  sep = "\t", header = T)
rownames(exmp)=exmp$geneName
exp1=intersect(rownames(exmp), rownames(zilb_counts))
exmp=exmp[exp1,]
zilb_fetal_counts=zilb_fetal_counts[exp1,]

zilb_fetal_counts$geneLengths=exmp$geneLengths
zilb_fetal_counts$geneName=rownames(zilb_fetal_counts)

for (x in zilb_meta$sample.id[which(zilb_meta$sample.type=="purified")]) {
  a = zilb_fetal_counts[,c("geneName","geneLengths",x)]
  write.table(a, paste0('cellsig_gut/fetal_bulk/',x,'.tsv'), sep="\t", row.names = F, col.names = T, quote = F)
}

gtex_meta=read.delim('/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/GTEX/Metadata.tsv')

rel_ids=gtex_meta$UniqueSampleID[which(gtex_meta$Tissue%in%c("Colon - Transverse","Colon - Sigmoid","Small Intestine - Terminal Ileum"))]
rel_gtex_meta=gtex_meta[which(gtex_meta$Tissue%in%c("Colon - Transverse","Colon - Sigmoid","Small Intestine - Terminal Ileum")),]
gtex_paths=paste0('/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/GTEX/frag/',rel_ids, '.tsv')
gtex_tissues=data.frame(path=gtex_paths, tissue=rel_gtex_meta$Tissue)

write.table(gtex_paths, 'cellsig_gut/gtex_paths.txt', sep = "\t", quote=F, row.names = F, col.names = F)
write.table(gtex_tissues, 'gtex_gut_tissues.txt', sep = "\t", row.names = F, col.names = T, quote=F)
