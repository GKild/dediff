paths=c('/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Xanthogranuloma/GOSH033/cellranger302_count_33315_CG_SB_NB8715412_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Xanthogranuloma/GOSH033/cellranger302_count_33315_CG_SB_NB8715413_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Xanthogranuloma/GOSH033/cellranger302_count_33315_CG_SB_NB8715414_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Ewings/GOSH032/cellranger302_count_33315_CG_SB_NB8715409_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Ewings/GOSH032/cellranger302_count_33315_CG_SB_NB8715410_GRCh38-1_2_0/filtered_feature_bc_matrix/',
        '/lustre/scratch117/casm/team274/my4/oldScratch/Data/Disease/Ewings/GOSH032/cellranger302_count_33315_CG_SB_NB8715409_GRCh38-1_2_0/filtered_feature_bc_matrix/')
names(paths)=c("GOSH33","GOSH33","GOSH33","GOSH32","GOSH32","GOSH32")


obj=Read10X(paths)

obj_srat=CreateSeuratObject(obj)

process_10x <- function(srat){
  srat = NormalizeData(srat)
  srat =FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  #srat = CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  srat = ScaleData(srat)
  srat = RunPCA(srat, npcs = 50)
  srat = FindNeighbors(srat, dims=1:50)
  srat = FindClusters(srat, resolution = 1)
  srat = RunUMAP(srat, dims=1:50, min.dist = 0.5, n.neighbors = 50)
  return(srat)
}

obj_processed=process_10x(obj_srat)


p1 <- DimPlot(obj_processed, group.by = "orig.ident")
p2 <- FeaturePlot(obj_processed, features=c("CD99","PAX7","CD68","PTPRC"))
plot_grid(p1,p2)
