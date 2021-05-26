rc2 = readRDS('/lustre/scratch119/realdata/mdt1/team274/my4/bulkRNAseq/rse_gene_TCGA.RDS')
len=read.table('matts-bits-of-code-master/cellSignalAnalysisV2/bulkData/SangerProjectTARGET_ALL_phase1/SangerProjectTARGET_ALL_phase1_SRR452447.tsv', sep = "\t", header = T)
obs=read.csv("onecs/obs.csv")
mDat = colData(rc2)


mDat$tissue = as.character(mDat$xml_tumor_tissue_site)

table(mDat$tissue)

mDat_liver=mDat[which(mDat$tissue=="Liver"),]

mDat_liver$xml_tumor_type

mDat_liver=mDat_liver[,colnames(mDat_liver)[which(!colnames(mDat_liver)%in%names(which(apply(mDat_liver, 2, function(x){sum(is.na(x))})==ncol(mDat_liver))))]]

saveRDS(mDat_liver, "mDat_liver.rds")


stjudes_hepatob=read.delim('/lustre/scratch119/realdata/mdt1/team274/ek12/Fetal_tissue_ness/Bulk_leukemia/StJudes/salmon_out/stJudes_bulk_txiConvGlen.csv')
stjudes_meta=read.delim('/lustre/scratch119/realdata/mdt1/team274/ek12/Fetal_tissue_ness/Bulk_leukemia/StJudes/SAMPLE_INFO.txt')
stjudes_meta2=read.delim('/lustre/scratch119/realdata/mdt1/team274/ek12/Fetal_tissue_ness/Bulk_leukemia/StJudes/SAMPLE_INFO_updated.txt')

st_judes_hepat_meta=stjudes_meta[which(stjudes_meta$attr_diagnosis=="Hepatoblastoma"),]
