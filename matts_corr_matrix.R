dat = rasa_updated_srat_processed_g1@assays$RNA@counts
#Assuming “clusters” is a vector of the same length as dat indicating which cell is of which type
dat = do.call(cbind,lapply(split(seq(ncol(dat)),rasa_updated_srat_processed_g1@meta.data$annotation_2),
                           function(e) rowSums(dat[,e,drop=FALSE])))
dat = t(t(dat)/colSums(dat))
out = cor(dat)
Heatmap(out)


Idents(rasa_updated_srat_processed_g1)=rasa_updated_srat_processed_g1@meta.data$annotation_2
rasa_updated_srat_processed_g1 = ScaleData(rasa_updated_srat_processed_g1, features=rownames(rasa_updated_srat_processed_g1))

avg_exp=AverageExpression(rasa_updated_srat_processed_g1, slot="scale.data")

cordat=cor(avg_exp$RNA)

Heatmap(cordat)


