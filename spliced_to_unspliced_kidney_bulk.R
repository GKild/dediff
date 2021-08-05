###kidney_velo###
library(velocyto.R)
library(loomR)
loom_paths=list.files('/lustre/scratch117/casm/team274/gk14/bulk_velo/velo_kidney_looms', full.names = T)

loom_lists=lapply(loom_paths, function(x){
  ldat = connect(x, skip.validate =T)
  gene.names = ldat[["row_attrs/Gene"]][]
  ambi=ldat[["layers/ambiguous"]][,]
  spliced=ldat[["layers/spliced"]][,]
  unspliced=ldat[["layers/unspliced"]][,]
  
  df=data.frame(ambi,spliced, unspliced)
  rownames(df)=make.unique(gene.names)
  return(df)
  
})


ldat_100_perc=connect("/lustre/scratch117/casm/team274/gk14/bulk_velo/subsample_test/PR36159g_full_bam/PR36159g_1EG38.loom", skip.validate =T)

ambi_100_perc=ldat_100_perc[["layers/ambiguous"]][,]
spliced_100_perc=ldat_100_perc[["layers/spliced"]][,]
unspliced_100_perc=ldat_100_perc[["layers/unspliced"]][,]

df_100_perc=data.frame(ambi_100_perc,spliced_100_perc, unspliced_100_perc)

sum(df_100_perc$ambi_100_perc)
sum(loom_lists$PR36159g_ROWZJ.loom$ambi)
sum(loom_lists$`5279STDY7449361_DQSBZ.loom`$ambi)/sum(loom_lists$PR36159g_ROWZJ.loom$ambi)

bla=sapply(loom_lists, function(x){
  return(colSums(x))
})

bla_melt=melt(bla)


ggplot(bla_melt, aes(x = Var2, y = value, fill=Var1))+
  geom_bar(position = "fill",stat="identity")+ theme(axis.text.x = element_text(angle = 90))



sum(df_100_perc$ambi_100_perc)
sum(df_100_perc$spliced_100_perc)
sum(df_100_perc$unspliced_100_perc)

bla_melt

df1=data.frame(Var1=c("ambi", "spliced", "unspliced"),
           Var2=c("PR36159g_GRCH37.5","PR36159g_GRCH37.5","PR36159g_GRCH37.5"),
           value=c(5143249, 21071172, 57871043))


bla2=rbind(bla_melt, df1)

ggplot(bla2, aes(x = Var2, y = value, fill=Var1))+
  geom_bar(position = "fill",stat="identity")+ theme(axis.text.x = element_text(angle = 90))


exposures = t(t(exposures)/colSums(exposures))
