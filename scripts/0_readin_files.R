# Read in files
library(Seurat)
library(Matrix)
library(ggplot2)
# file.names = list.files('your_path_to_data/')
# Use shared genes 
gene.names = NULL
counts = NULL
meta.data = NULL
for(f in file.names){
  temp = readRDS(paste0('your_path_to_data/', f))
  temp = temp[, temp$percent.mt < 10]
  temp = temp[, temp$doubletFinder == 'singlets']
  if(is.null(gene.names)){
    gene.names = rownames(temp)
  }else{
    gene.names = intersect(gene.names, rownames(temp))
  }
  if(is.null(counts)){
    counts = temp@assays$RNA$counts
  }else{
    counts = cbind(counts[gene.names, ], temp@assays$RNA$counts[gene.names, ])
  }
  if(is.null(meta.data)){
    meta.data = temp@meta.data[, c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'doubletFinder')]
  }else{
    meta.data = rbind(meta.data, temp@meta.data[, c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'doubletFinder')])
  }
  cat(f, '\n')
}
rm(temp)

rownames(meta.data) = paste0(meta.data$orig.ident, '-', rownames(meta.data))
colnames(counts) = rownames(meta.data)

adnp.all.seu = CreateSeuratObject(counts = counts, meta.data = meta.data, min.cells = 3, min.features = 200)
background = c('1599', '1599', '1599', 'PGP1', 'PGP1', 'PGP1', 'WTC11', 'WTC11', 'WTC11')
genotype = c('Het', 'Homo', 'WT', 'Het', 'Homo', 'WT', 'Het', 'Homo', 'WT')
adnp.all.seu$background = factor(background[as.numeric(adnp.all.seu$orig.ident)], levels = c('1599', 'PGP1', 'WTC11'))
adnp.all.seu$genotype = factor(genotype[as.numeric(adnp.all.seu$orig.ident)], levels = c('WT', 'Het', 'Homo'))

VlnPlot(adnp.all.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

adnp.all.seu <- NormalizeData(adnp.all.seu, normalization.method = "LogNormalize", scale.factor = 10000)

adnp.all.seu <- FindVariableFeatures(adnp.all.seu, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adnp.all.seu), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(adnp.all.seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(adnp.all.seu)
adnp.all.seu <- ScaleData(adnp.all.seu, features = all.genes)

# Add genotype and background in meta.data
background = c('1599', '1599', '1599', 'PGP1', 'PGP1', 'PGP1', 'WTC11', 'WTC11', 'WTC11')
genotype = c('Het', 'Homo', 'WT', 'Het', 'Homo', 'WT', 'Het', 'Homo', 'WT')
adnp.all.seu$background = factor(background[as.numeric(adnp.all.seu$orig.ident)], levels = c('1599', 'PGP1', 'WTC11'))
adnp.all.seu$genotype = factor(genotype[as.numeric(adnp.all.seu$orig.ident)], levels = c('WT', 'Het', 'Homo'))
saveRDS(adnp.all.seu, file = 'Data_store/adnp_all_seu.rds')

adnp.all.seu <- RunPCA(adnp.all.seu, features = VariableFeatures(object = adnp.all.seu))
DimPlot(adnp.all.seu, reduction = "pca", group.by = 'orig.ident')
DimPlot(adnp.all.seu, reduction = "pca", group.by = 'genotype')


ElbowPlot(adnp.all.seu, ndims = 50)

adnp.all.seu <- RunUMAP(adnp.all.seu, dims = 1:30)
adnp.all.seu <- RunTSNE(adnp.all.seu, dims = 1:30)
umap.coord = adnp.all.seu@reductions$umap
tsne.coord = adnp.all.seu@reductions$tsne
save(umap.coord, file = '9celline_umap_coord.RData')
save(tsne.coord, file = '9cellline_tsne_coord.RData')


DimPlot(adnp.all.seu, reduction = 'tsne', group.by = 'orig.ident')
DimPlot(adnp.all.seu, reduction = 'tsne', group.by = 'background')
DimPlot(adnp.all.seu, reduction = 'tsne', group.by = 'genotype')

DimPlot(adnp.all.seu, reduction = 'umap', group.by = 'orig.ident')
pdf('plots/UMAP_9_cellines_background.pdf', width = 10, height = 4)
DimPlot(adnp.all.seu, reduction = 'umap', group.by = 'background') + facet_wrap( ~ background) + scale_color_brewer(palette = 'Set2')
dev.off()

pdf('plots/UMAP_9_celllines_genotypes.pdf', width = 10, height = 10)
DimPlot(adnp.all.seu, reduction = 'umap', group.by = 'genotype', split.by = 'background') + facet_grid(background ~ genotype) + scale_color_brewer(palette = 'Set1')
dev.off()

exprs.perc = colMeans(adnp.all.seu@assays$RNA$counts > 0)
exprs.perc.for.hist = data.frame(exprs.percent = exprs.perc, genotype = adnp.all.seu$genotype, background = adnp.all.seu$background)

pdf('plots/UMAP_9_celllines_hist_exprs_percent.pdf', width = 10, height = 10)
ggplot(exprs.perc.for.hist, aes(exprs.percent)) + geom_histogram() + facet_grid(genotype ~ background) + theme_minimal(base_size = 15)
dev.off()

# We are going to move WTC11-Homo from further analysis right now. 

adnp.8.seu = adnp.all.seu[, adnp.all.seu$orig.ident != 'WTC11-Homo2']
rm(adnp.all.seu)



