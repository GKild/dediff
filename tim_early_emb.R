#Tims early embryo
early_emb=readRDS('/lustre/scratch117/casm/team274/tc16/Embryo/embryo_integrated_allembryos_filtered.Rdata')

DimPlot(early_emb, group.by = "Phase")

table(Idents(early_emb))
