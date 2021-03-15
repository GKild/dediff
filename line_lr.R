line_counts=read.delim('crc_line/GSE149224_RSH.all.counts.txt', sep = " ", header = T)
line_meta=read.csv('crc_line/GSE149224_meta.information.csv')

line_mtx=as.matrix(line_counts)
line_mtx=Matrix(line_mtx, sparse = T)

ref_comb=merge(srat_comb_nocc, fetal_epi)
ref_comb@meta.data$celltype=c(srat_comb_nocc@meta.data$new_ident, fetal_epi@meta.data$exp_annot)

shared_rows=intersect(rownames(ref_comb), rownames(line_mtx))

ref_shared=ref_comb@assays$RNA@counts[shared_rows,]
line_mtx_shared=line_mtx[shared_rows,]

source('logisticRegression.R')

ref_tr=trainModel(ref_shared, ref_comb@meta.data$celltype, workers=NULL)

line_pred=predictSimilarity(ref_tr, line_mtx_shared, line_meta$dose)

similarityHeatmap(line_pred)
