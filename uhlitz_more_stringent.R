process_seurat2 =function(srat){
  srat=subset(srat, subset = pMT < 0.5)
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 10)
  srat = FindNeighbors(srat, dims=1:10)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:10, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}


uhlitz_epi_tumor_mt_filt=process_seurat2(uhlitz_epi_tumor)

row_ord=c("Ent. Prog. 1",  "Ent. Prog. 2",  "Enterocytes 1", "Enterocytes 2", "Enterocytes 3", "Goblet",
          "Imma. Gob. 1" , "Imma. Gob. 2","Tuft", "Stem" ,"Stem/TA 1","Stem/TA 2","Stem/TA 3","TC1","TC2","TC3","TC4") 
DimPlot(uhlitz_epi_tumor_mt_filt, group.by = "cell_type_epi_custom", label=T)

ps_uhlitz_tum_mt=predictSimilarity(ref_train, uhlitz_epi_tumor_mt_filt@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor_mt_filt@meta.data$cell_type_epi_custom)
ps_uhlitz_tum_mt_merged=predictSimilarity(ref_train_merged, uhlitz_epi_tumor_mt_filt@assays$RNA@counts[comm_genes,], uhlitz_epi_tumor_mt_filt@meta.data$cell_type_epi_custom)

similarityHeatmap(ps_uhlitz_tum_mt, column_order=col_ord, row_order=row_ord)


large_int_rasa=readRDS('rasa_stuff/large_intestine_rasa.rds')



uhlitz_epi_new=process_seurat2(uhlitz_epi)
uhlitz_epi_new@meta.data$new_ident=paste0(uhlitz_epi_new@meta.data$sample_origin, "_", uhlitz_epi_new@meta.data$cell_type_epi_simple)

uhlitz_epi_new_nocc=subset(uhlitz_epi_new, Phase=="G1")
epi_new_ps=predictSimilarity(ref_train, uhlitz_epi_new@assays$RNA@counts[comm_genes,], uhlitz_epi_new@meta.data$new_ident)
epi_new_nocc_ps=predictSimilarity(ref_train, uhlitz_epi_new_nocc@assays$RNA@counts[comm_genes,], uhlitz_epi_new_nocc@meta.data$new_ident)

similarityHeatmap(epi_new_ps, column_order=col_ord)
similarityHeatmap(epi_new_nocc_ps, column_order=col_ord)

fetal_crypt=subset(lif_nocc_seurat, exp_annot%in%c("fetal_crypt"))

new_int=intersect(comm_genes, rownames(fetal_crypt))

bound_mtx=cbind(uhlitz_epi_new@assays$RNA@counts[new_int, ], fetal_crypt@assays$RNA@counts[new_int, ])

celltypes=c(uhlitz_epi_new@meta.data$new_ident, fetal_crypt@meta.data$exp_annot)

new_ps=predictSimilarity(ref_train, bound_mtx, celltypes, logits = F)
new_ps2=predictSimilarity(ref_train_merged, bound_mtx, celltypes, logits = F)

similarityHeatmap(new_ps, column_order=col_ord)
similarityHeatmap(new_ps2, column_order=col_ord_merged)

new_ps=as.data.frame(new_ps)

new_ps$celltypes=celltypes

new_ps=new_ps[new_ps$celltypes%in%c("Normal_Stem", "Tumor_TC1", "Tumor_TC2", "Tumor_TC3", "Tumor_TC4", "fetal_crypt", "fetal_B/Plasma"),]

new_ps=new_ps[, c("Stem", "fetal_crypt", "celltypes")]

new_ps_melted=melt(new_ps, id.vars = "celltypes")

new_ps_melted$variable=as.character(new_ps_melted$variable)
new_ps_melted$variable[which(new_ps_melted$variable=="fetal_crypt")]="Fetal Stem"
new_ps_melted$variable[which(new_ps_melted$variable=="Stem")]="Adult Stem"


new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Normal_Stem")]="Normal Stem"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC1")]="Tumor Cluster 1"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC2")]="Tumor Cluster 2"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC3")]="Tumor Cluster 3"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC4")]="Tumor Cluster 4"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="fetal_crypt")]="Fetal Stem"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="fetal_B/Plasma")]="Fetal B/Plasma"

new_ps_melted$celltypes=factor(new_ps_melted$celltypes, levels=c("Normal Stem", "Fetal Stem", "Tumor Cluster 1","Tumor Cluster 2",
                                                                 "Tumor Cluster 3", "Tumor Cluster 4", "Fetal B/Plasma"))

ggplot(data = new_ps_melted, aes(x=celltypes, y=value, fill=variable)) +geom_boxplot()

ggplot(data = new_ps_melted, aes(x=celltypes,y=value, fill=variable)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) + theme_bw() + scale_fill_manual(values=c("#114477", "#117744"))+
  geom_quasirandom(mapping = aes(x=celltypes,y=value, group=variable), 
                   dodge.width=.8, shape=19, cex=0.5, alpha=0.1) +
  labs(x="Cell type/cluster", y="Probability of match against reference",fill="Reference cell type")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))



###############
#new ref with SCPs 

adr_all=readRDS('/lustre/scratch117/casm/team274/gk14/cellxgene/adr_all.rds')
scps_mtx=subset(adr_all, Annotation=="SCPs")


ref_subset=subset(ref_comb, celltype%in%c("Stem" ,  "fetal_crypt", "Enterocytes" ))


int1=intersect(rownames(scps_mtx), rownames(ref_subset))
int2=intersect(int1, rownames(uhlitz_epi_new))

mtx_with_ctrls=cbind(uhlitz_epi_new@assays$RNA@counts[int2, ],
                     scps_mtx@assays$RNA@counts[int2, ],
                     ref_subset@assays$RNA@counts[int2,])
celltypes_with_ctrls=c(uhlitz_epi_new@meta.data$new_ident, scps_mtx@meta.data$Annotation, ref_subset@meta.data$celltype)

dim(mtx_with_ctrls)
length(celltypes_with_ctrls)

#####train a new model for this

ref_train_ctrls=trainModel(ref_comb@assays$RNA@counts[int2,], ref_comb@meta.data$celltype, workers=NULL)

ps_ctrls=predictSimilarity(ref_train_ctrls, mtx_with_ctrls, logits = F)

ps_ctrls2=predictSimilarity(ref_train_ctrls, mtx_with_ctrls, celltypes_with_ctrls, logits = T)
similarityHeatmap(ps_ctrls2, column_order=col_ord)
new_ps=as.data.frame(ps_ctrls)

new_ps$celltypes=celltypes_with_ctrls

new_ps=new_ps[new_ps$celltypes%in%c("Stem", "Tumor_TC1", "Tumor_TC2", "Tumor_TC3", "Tumor_TC4", "fetal_crypt", "Enterocytes", "SCPs"),]

new_ps=new_ps[, c("Stem", "fetal_crypt", "celltypes")]

new_ps_melted=melt(new_ps, id.vars = "celltypes")

new_ps_melted$variable=as.character(new_ps_melted$variable)
new_ps_melted$variable[which(new_ps_melted$variable=="fetal_crypt")]="Fetal Stem"
new_ps_melted$variable[which(new_ps_melted$variable=="Stem")]="Adult Stem"


new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC1")]="Tumor Cluster 1"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC2")]="Tumor Cluster 2"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC3")]="Tumor Cluster 3"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC4")]="Tumor Cluster 4"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="fetal_crypt")]="Fetal Stem"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Enterocytes")]="Adult Enterocytes"

new_ps_melted$celltypes=factor(new_ps_melted$celltypes, levels=c("Stem", "Fetal Stem", "Tumor Cluster 1","Tumor Cluster 2",
                                                                 "Tumor Cluster 3", "Tumor Cluster 4", "SCPs", "Adult Enterocytes"))

ggplot(data = new_ps_melted, aes(x=celltypes, y=value, fill=variable)) +geom_boxplot()

ggplot(data = new_ps_melted, aes(x=celltypes,y=value, fill=variable)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) + theme_bw() + scale_fill_manual(values=c("#114477", "#117744"))+
  geom_quasirandom(mapping = aes(x=celltypes,y=value, group=variable), 
                   dodge.width=.8, shape=19, cex=0.5, alpha=0.1) +
  labs(x="Cell type/cluster", y="Probability of match against reference",fill="Reference cell type")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))

new_ps_melted$celltypes=as.character(new_ps_melted$celltypes)
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC1")]="Tumor"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC2")]="Tumor"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC3")]="Tumor"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Tumor_TC4")]="Tumor"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="fetal_crypt")]="Fetal Stem"
new_ps_melted$celltypes[which(new_ps_melted$celltypes=="Enterocytes")]="Adult Enterocytes"

new_ps_melted$celltypes=factor(new_ps_melted$celltypes, levels=c("Stem", "Fetal Stem", "Tumor", "SCPs", "Adult Enterocytes"))

ggplot(data = new_ps_melted, aes(x=celltypes, y=value, fill=variable)) +geom_boxplot()

ggplot(data = new_ps_melted, aes(x=celltypes,y=value, fill=variable)) +
  geom_boxplot(outlier.shape=NA, alpha=0.7) + theme_bw() + scale_fill_manual(values=c("#114477", "#117744"))+
  geom_quasirandom(mapping = aes(x=celltypes,y=value, group=variable), 
                   dodge.width=.8, shape=19, cex=0.5, alpha=0.1) +
  labs(x="Cell type/cluster", y="Probability of match against reference",fill="Reference cell type")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text=element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1))
