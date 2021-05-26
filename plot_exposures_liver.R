mDat_liver=readRDS('mDat_liver.rds')

paths=paste0("/lustre/scratch117/casm/team274/my4/oldScratch/ProjectsExtras/CellDeconvolution/Results/bulkRNAseq/TCGA/frag/TCGA_", rownames(mDat_liver), ".tsv")

write.table(paths, "cellsig_liver/liver_tcga_paths.txt", sep = "\t", row.names = F, col.names = F, quote = F)
source('matts-bits-of-code-master/cellSignalAnalysisV2/cellSignalAnalysis.R')

cycling_tcga='cellsig_liver/out_tcga_cycling/OutRun_fitExposures.tsv'
fit_cycling_tcga=normaliseExposures(cycling_tcga)


plotExposures(fit_cycling_tcga, show_column_names=F, column_split=mDat_liver$cgc_sample_sample_type)

noncycling_tcga='cellsig_liver/out_tcga_noncycling/OutRun_fitExposures.tsv'
fit_noncycling_tcga=normaliseExposures(noncycling_tcga)


plotExposures(fit_noncycling_tcga, show_column_names=F, column_split=mDat_liver$cgc_sample_sample_type)

