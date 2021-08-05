library(velocyto.R)
ldat_1_perc = connect("/lustre/scratch117/casm/team274/gk14/bulk_velo/PR36159g_1perc_sub/PR36159g_3m_reads_5S35M.loom", skip.validate =T)
ldat_10_perc= connect("/lustre/scratch117/casm/team274/gk14/bulk_velo/PR36159g_10perc_sub/PR36159g_30m_reads_Q8EUX.loom", skip.validate =T)
ldat_20_perc=connect("/lustre/scratch117/casm/team274/gk14/bulk_velo/PR36159g_20perc_sub/PR36159g_60m_reads_Z9YSD.loom", skip.validate =T)
ldat_50_perc=connect("/lustre/scratch117/casm/team274/gk14/bulk_velo/PR36159g_50perc_sub/PR36159g_150m_reads_9BTCC.loom", skip.validate =T)
ldat_100_perc=connect("/lustre/scratch117/casm/team274/gk14/bulk_velo/PR36159g_full_bam/PR36159g_1EG38.loom", skip.validate =T)


gene.names = ldat_1_perc[["row_attrs/Gene"]][]

ambi_1_perc=ldat_1_perc[["layers/ambiguous"]][,]
spliced_1_perc=ldat_1_perc[["layers/spliced"]][,]
unspliced_1_perc=ldat_1_perc[["layers/unspliced"]][,]

df_1_perc=data.frame(ambi_1_perc,spliced_1_perc, unspliced_1_perc)


ambi_10_perc=ldat_10_perc[["layers/ambiguous"]][,]
spliced_10_perc=ldat_10_perc[["layers/spliced"]][,]
unspliced_10_perc=ldat_10_perc[["layers/unspliced"]][,]

df_10_perc=data.frame(ambi_10_perc,spliced_10_perc, unspliced_10_perc)

ambi_20_perc=ldat_20_perc[["layers/ambiguous"]][,]
spliced_20_perc=ldat_20_perc[["layers/spliced"]][,]
unspliced_20_perc=ldat_20_perc[["layers/unspliced"]][,]

df_20_perc=data.frame(ambi_20_perc,spliced_20_perc, unspliced_20_perc)

ambi_50_perc=ldat_50_perc[["layers/ambiguous"]][,]
spliced_50_perc=ldat_50_perc[["layers/spliced"]][,]
unspliced_50_perc=ldat_50_perc[["layers/unspliced"]][,]

df_50_perc=data.frame(ambi_50_perc,spliced_50_perc, unspliced_50_perc)

ambi_100_perc=ldat_100_perc[["layers/ambiguous"]][,]
spliced_100_perc=ldat_100_perc[["layers/spliced"]][,]
unspliced_100_perc=ldat_100_perc[["layers/unspliced"]][,]

df_100_perc=data.frame(ambi_100_perc,spliced_100_perc, unspliced_100_perc)


colnames(df_1_perc)=c("ambiguous", "spliced", "unspliced")
colnames(df_10_perc)=c("ambiguous", "spliced", "unspliced")
colnames(df_20_perc)=c("ambiguous", "spliced", "unspliced")
colnames(df_50_perc)=c("ambiguous", "spliced", "unspliced")
colnames(df_100_perc)=c("ambiguous", "spliced", "unspliced")



df_1_perc$subsample="1_percent"
df_10_perc$subsample="10_percent"
df_20_perc$subsample="20_percent"
df_50_perc$subsample="50_percent"
df_100_perc$subsample="100_percent"


final_df=rbind(df_1_perc, df_10_perc, df_20_perc, df_50_perc, df_100_perc)

final_df_melt=melt(final_df)
ggplot(final_df_melt, aes(x = subsample, y = value, fill=variable))+
  geom_col(position="fill", width = 0.7)

plot(df_100_perc$unspliced+df_100_perc$spliced, df_100_perc$spliced/(df_100_perc$unspliced+df_100_perc$spliced), log="x", cex=0.01)



plot(hist(df_100_perc$spliced/(df_100_perc$unspliced+df_100_perc$spliced), n=100))
