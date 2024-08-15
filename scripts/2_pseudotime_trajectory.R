# Narrow down to the cluster 0, 1, 2

# load data and add covariates and estimates
adnp.sub = readRDS('Data_store/adnp_all_seu_light.rds')
adnp.male = adnp.sub[, adnp.sub$background!='1599' & adnp.sub$orig.ident != 'WTC11-Homo2']

load('Data_store/adnp_male_umap_coord_meta_data.RData')
adnp.male@reductions$umap = umap.coord.male
adnp.male$seurat_clusters = meta.data.male$seurat_clusters

adnp.cds = readRDS('Data_store/adnp_cds_with_pt.rds')
adnp.cds.new <- order_cells(adnp.cds)
saveRDS(adnp.cds.new, file = 'Data_store/adnp_cds_with_pt_cluster0_root.rds')

adnp.male$pseudotime = adnp.cds@principal_graph_aux$UMAP$pseudotime

umap.select = which(adnp.male@reductions$umap@cell.embeddings[, 1] < 10 & 
                      adnp.male@reductions$umap@cell.embeddings[, 2] < 5 & 
                      adnp.male@reductions$umap@cell.embeddings[, 1] - 4*adnp.male@reductions$umap@cell.embeddings[, 2] > -8 & 
                      adnp.male$pseudotime < Inf)

adnp.male$pseudotime.new = adnp.cds.new@principal_graph_aux$UMAP$pseudotime
umap.select.new = which(adnp.male@reductions$umap@cell.embeddings[, 1] < 10 & 
                          adnp.male@reductions$umap@cell.embeddings[, 2] < 5 & 
                          adnp.male@reductions$umap@cell.embeddings[, 1] - 4*adnp.male@reductions$umap@cell.embeddings[, 2] > -8 & 
                          adnp.male$pseudotime.new < Inf)

adnp.male.sub = adnp.male[, umap.select]
adnp.male.sub.new = adnp.male[, umap.select.new]

adnp.graph.test <- graph_test(adnp.cds, neighbor_graph="knn", cores=8)
adnp.graph.test <- adnp.graph.test[order(adnp.graph.test$morans_I, decreasing = T), ]
adnp.pt.id <- row.names(subset(adnp.graph.test, morans_I > 0.2 & q_value < 0.01))
pdf(file = 'pseudotime_focused_root_cluster0.pdf', width = 8, height = 6)
FeaturePlot(adnp.male.sub.new, features = 'pseudotime.new')
dev.off()

adnp.cds.new = readRDS('Data_store/adnp_cds_with_pt_cluster0_root.rds')
plot_cells(adnp.cds.new, color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


adnp.graph.test <- graph_test(adnp.cds.new, neighbor_graph = 'knn', cores = 8)
adnp.graph.test <- adnp.graph.test[order(adnp.graph.test$morans_I, decreasing = T), ]
adnp.pt.id <- row.names(subset(adnp.graph.test, morans_I > 0.2 & q_value < 0.01))

pdf(file = 'pseudotime_focused_root_cluster0.pdf', width = 8, height = 6)
FeaturePlot(adnp.male.sub.new, features = 'pseudotime.new')
dev.off()


### Start with cluster 1 
adnp.cds.new <- readRDS('Data_store/adnp_cds_with_pt_cluster1_root.rds')
adnp.graph.test.cluster1 = graph_test(adnp.cds.new, neighbor_graph = 'knn', cores = 8)
adnp.graph.test.cluster1 <- adnp.graph.test.cluster1[order(adnp.graph.test.cluster1$morans_I, decreasing = T), ]
adnp.pt.id <- row.names(subset(adnp.graph.test.cluster1, morans_I > 0.5 & q_value < 0.01))
save(adnp.graph.test.cluster1, adnp.pt.id, file = 'Data_store/adnp_graph_test_cluster1.RData')

adnp.male$pseudotime = adnp.cds.new@principal_graph_aux$UMAP$pseudotime
umap.select.new = which(adnp.male@reductions$umap@cell.embeddings[, 1] < 10 & 
                          adnp.male@reductions$umap@cell.embeddings[, 2] < 5 & 
                          adnp.male@reductions$umap@cell.embeddings[, 1] - 4*adnp.male@reductions$umap@cell.embeddings[, 2] > -8 & 
                          adnp.male$pseudotime < Inf)

adnp.male.sub.new = adnp.male[, umap.select.new]

cor.pt.cluster1 = apply(adnp.male.sub.new@assays$RNA$data, 1, cor, y = adnp.male.sub.new$pseudotime)
save(cor.pt.cluster1, file = 'Data_store/cor_pt_genes_cluster1.RData')
FeaturePlot(adnp.male.sub.new, features = names(which(abs(cor.pt.cluster1) > 0.5)) )


load('Data_store/adnp_graph_test_cluster1.RData')
# row.names(subset(adnp.graph.test.cluster1, morans_I > 0.4 & q_value < 0.01))
#[1] "PHOX2B"     "MALAT1"     "AC105389.3" "GRID2"      "BNC2"       "PTPRG"      "UNC5D"      "TAC1"       "PLD5"       "FOXP2"      "SNTG1"      "FAM155A"    "ROBO1" 

# add new label, focus on only the mature cells
adnp.male.sub.new$cluster.new = 'progenitor'
adnp.male.sub.new$cluster.new[adnp.male.sub.new$pseudotime.new > 10 & adnp.male.sub.new$seurat_clusters == 2] = 'mature.WT'
adnp.male.sub.new$cluster.new[adnp.male.sub.new$pseudotime.new > 10 & adnp.male.sub.new$seurat_clusters == 1] = 'mature.Mut'

temp = FindMarkers(adnp.male.sub.new, ident.1 = 'mature.WT',ident.2 = 'mature.Mut', group.by = 'cluster.new')
temp.up = temp[temp$avg_log2FC > 1 & temp$p_val_adj < 0.01, ]
temp.down = temp[temp$avg_log2FC < -1 & temp$p_val_adj < 0.01, ]
write.csv(temp.up, file = 'Data_store/mature_compare_up.csv')
write.csv(temp.down, file = 'Data_store/mature_compare_down.csv')

write.table(rownames(temp.up), file = 'Data_store/up_de_mature_wilcox.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(temp.down), file = 'Data_store/down_de_mature_wilcox.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

df.pt = NULL
for(g in adnp.pt.id[1:10]){
  df.temp = data.frame(pt = adnp.male$pseudotime, 
                       exprs = adnp.male@assays$RNA@layers$data[match(g, rownames(adnp.male)), ], 
                       gene.name = g)
  df.temp = df.temp[df.temp$pt < Inf, ]
  df.pt = rbind(df.pt, df.temp)
}
df.pt$gene.name = factor(df.pt$gene.name, levels = adnp.pt.id[1:10])
ggplot(df.pt, aes(pt, exprs, color = pseudotime)) + geom_point() + 
  geom_spline(span = 2) + 
  facet_wrap(~ gene.name, ncol = 2) + theme_minimal(base_size = 14) + scale_color_viridis(option="inferno")
