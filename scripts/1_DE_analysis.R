# de_analysis background
library(Seurat)

# Work with adnp.8.seu first
# Need pre-filtering for genes no matter methods we are going to use.
adnp.all.seu = readRDS('Data_store/adnp_all_seu_light.rds')
adnp.sub = adnp.all.seu[, adnp.all.seu$orig.ident != 'WTC11-Homo2']
# Expressed in at 5% of the cells
# selected.genes = rownames(adnp.sub)[rowMeans(adnp.sub@assays$RNA$counts[, adnp.sub$genotype != 'Homo'] > 0) > 0.05]
# 10517 genes left

# Start with the Seurat FindMarkers WT vs. Het for 3 backgrounds separately
background = c('1599', 'PGP1', 'WTC11')
genotype = c('WT', 'Het', 'Homo')
#### MAST
# list.de.mast = list()
# for(s in background){
#  temp = adnp.sub[, adnp.sub$background == s & adnp.sub$genotype != 'Homo']
#  de.mast = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Het', test.use = 'MAST', features = selected.genes, group.by = 'genotype', logfc.threshold = 0)
#  list.de.mast[[s]] = de.mast
# }
# save(list.de.mast, file = 'Data_store/list_de_seurat_mast.RData')

#for(s in background[1:2]){
#  temp = adnp.sub[, adnp.sub$background == s]
#  temp$gt = factor('WT', levels = c('WT', 'Mut'))
#  temp$gt[temp$genotype != 'WT'] = 'Mut'
#  
#  de.mast = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Mut', test.use = 'MAST', features = selected.genes, group.by = 'gt', logfc.threshold = 0) 
#  save(de.mast, file = paste0('Data_store/de_mast_with_Homo_', s, '.RData'))
#}

temp = adnp.sub[, adnp.sub$background != '1599']
selected.genes = rownames(temp)[rowMeans(temp@assays$RNA$counts > 0) > 0.05]

de.mast.male = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Het', test.use = 'MAST', features = selected.genes, group.by = 'genotype', logfc.threshold = 0)

temp$gt = factor('WT', levels = c('WT', 'Mut'))
temp$gt[temp$genotype != 'WT'] = 'Mut'
de.mast.male.Homo = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Mut', test.use = 'MAST', features = selected.genes, group.by = 'gt', logfc.threshold = 0) 
save(de.mast.male, de.mast.male.Homo, file = 'Data_store/de_mast_male_comb.RData')


########## DE for each cell line ##################
# DE analysis
# Start with each cell line separately
# We tested two situations: WT vs. Het; and WT vs. (Het and Homo)
# I have removed WTC11-Homo2 b/c the expression percentage is not good at all. 
# I believe there are something wrong with this dataset
library(Seurat)

# Work with adnp.8.seu first
# Need pre-filtering for genes no matter methods we are going to use.
adnp.all.seu = readRDS('Data_store/adnp_all_seu_light.rds')
adnp.sub = adnp.all.seu[, adnp.all.seu$orig.ident != 'WTC11-Homo2']

hist(rowMeans(adnp.sub@assays$RNA$counts[, adnp.sub$genotype != 'Homo'] > 0), breaks = 100)
# Expressed in at 5% of the cells
selected.genes = rownames(adnp.sub)[rowMeans(adnp.sub@assays$RNA$counts[, adnp.sub$genotype != 'Homo'] > 0) > 0.05]
# 10517 genes left

# Start with the Seurat FindMarkers WT vs. Het for 3 backgrounds separately
background = c('1599', 'PGP1', 'WTC11')
genotype = c('WT', 'Het', 'Homo')
#list.de.mast = list()
list.de.wilcox = list()
for(s in background){
  temp = adnp.sub[, adnp.sub$background == s & adnp.sub$genotype != 'Homo']
  #de.mast = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Het', test.use = 'MAST', group.by = 'genotype', logfc.threshold = 0)
  #list.de.mast[[s]] = de.mast
  de.wilcox = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Het', features = selected.genes, group.by = 'genotype', logfc.threshold = 0)
  list.de.wilcox[[s]] = de.wilcox
}

df.log2FC = data.frame(genes = selected.genes, F1599_log2FC = list.de.wilcox$`1599`[selected.genes, 'avg_log2FC'], 
           PGP1_log2FC = list.de.wilcox$PGP1[selected.genes, 'avg_log2FC'], 
           WTC11_log2FC = list.de.wilcox$WTC11[selected.genes, 'avg_log2FC'])

pdf('plots/scatter_logFC_wilcox.pdf', width = 8, height = 8)
plot(df.log2FC[, -1], pch = '.')
dev.off()

gene.names = intersect(rownames(list.de.wilcox$`1599`), rownames(list.de.wilcox$PGP1))
plot(list.de.wilcox$`1599`[gene.names, ]$avg_log2FC, list.de.wilcox$PGP1[gene.names, ]$avg_log2FC)
cor(list.de.wilcox$`1599`[gene.names, ]$avg_log2FC, list.de.wilcox$PGP1[gene.names, ]$avg_log2FC)
# [1] -0.07273881
gene.names = intersect(rownames(list.de.wilcox$`1599`), rownames(list.de.wilcox$WTC11))
plot(list.de.wilcox$`1599`[gene.names, ]$avg_log2FC, list.de.wilcox$WTC11[gene.names, ]$avg_log2FC)
cor(list.de.wilcox$`1599`[gene.names, ]$avg_log2FC, list.de.wilcox$WTC11[gene.names, ]$avg_log2FC)
# [1] 0.06057116
gene.names = intersect(rownames(list.de.wilcox$PGP1), rownames(list.de.wilcox$WTC11))
plot(list.de.wilcox$PGP1[gene.names, ]$avg_log2FC, list.de.wilcox$WTC11[gene.names, ]$avg_log2FC)
cor(list.de.wilcox$PGP1[gene.names, ]$avg_log2FC, list.de.wilcox$WTC11[gene.names, ]$avg_log2FC)
# [1] 0.8242791


# Read in bulk data DE results
de.bulk = read.csv('hiPSCs_GABA_ADNP.csv')
# 18813 genes
# length(intersect(de.bulk$Symbol, selected.genes))
# [1] 9614

df.log2FC$bulk_logFC = NA
common.genes = intersect(de.bulk$Symbol, selected.genes)
df.log2FC$bulk_logFC[match(common.genes, selected.genes)] = de.bulk$logFC[match(common.genes, de.bulk$Symbol)]
bulk_logFC = df.log2FC$bulk_logFC
a = df.log2FC$F1599_log2FC
cor(a[!is.na(a+bulk_logFC)], bulk_logFC[!is.na(a+bulk_logFC)])
# [1] -0.1445868
a = df.log2FC$PGP1_log2FC
cor(a[!is.na(a+bulk_logFC)], bulk_logFC[!is.na(a+bulk_logFC)])
# [1] 0.03433324
a = df.log2FC$WTC11_log2FC
cor(a[!is.na(a+bulk_logFC)], bulk_logFC[!is.na(a+bulk_logFC)])
# [1] 0.007702086

# The pi1 analysis 
bulk.de.genes = de.bulk$Symbol[de.bulk$adj.P.Val < 0.05]
bulk.de.genes = bulk.de.genes[bulk.de.genes!= '--']
# 1397 genes
pi1.discover.bulk = unlist(lapply(list.de.wilcox, function(x){
  1-propTrueNull(x[intersect(bulk.de.genes, rownames(x)), 'p_val'])
}))
# 1599      PGP1     WTC11 
# 0.8608303 0.8527850 0.8760821
pi1.discover.sc = unlist(lapply(list.de.wilcox, function(x){
  deg = rownames(x)[x$p_val_adj < 0.05]
  1-propTrueNull(de.bulk[match(intersect(deg, de.bulk$Symbol), de.bulk$Symbol), 'P.Value'])
}))

# 1599      PGP1     WTC11 
# 0.2695400 0.2641266 0.2693396 
#### MAST
list.de.mast = list()
for(s in background){
  temp = adnp.sub[, adnp.sub$background == s & adnp.sub$genotype != 'Homo']
  de.mast = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Het', test.use = 'MAST', features = selected.genes, group.by = 'genotype', logfc.threshold = 0)
  list.de.mast[[s]] = de.mast
}

load('Data_store/list_de_seurat_mast.RData')
# list.de.mast
df.log2FC = data.frame(genes = selected.genes, F1599_log2FC = list.de.mast$`1599`[selected.genes, 'avg_log2FC'], 
                       PGP1_log2FC = list.de.mast$PGP1[selected.genes, 'avg_log2FC'], 
                       WTC11_log2FC = list.de.mast$WTC11[selected.genes, 'avg_log2FC'])

pdf('plots/scatter_logFC_mast.pdf', width = 8, height = 8)
plot(df.log2FC[, -1], pch = '.')
dev.off()

gene.names = intersect(rownames(list.de.mast$`1599`), rownames(list.de.mast$PGP1))
plot(list.de.mast$`1599`[gene.names, ]$avg_log2FC, list.de.mast$PGP1[gene.names, ]$avg_log2FC)
cor(list.de.mast$`1599`[gene.names, ]$avg_log2FC, list.de.mast$PGP1[gene.names, ]$avg_log2FC)
# [1] -0.07273881
gene.names = intersect(rownames(list.de.mast$`1599`), rownames(list.de.mast$WTC11))
plot(list.de.mast$`1599`[gene.names, ]$avg_log2FC, list.de.mast$WTC11[gene.names, ]$avg_log2FC)
cor(list.de.mast$`1599`[gene.names, ]$avg_log2FC, list.de.mast$WTC11[gene.names, ]$avg_log2FC)
# [1] 0.06057116
gene.names = intersect(rownames(list.de.mast$PGP1), rownames(list.de.mast$WTC11))
plot(list.de.mast$PGP1[gene.names, ]$avg_log2FC, list.de.mast$WTC11[gene.names, ]$avg_log2FC)
cor(list.de.mast$PGP1[gene.names, ]$avg_log2FC, list.de.mast$WTC11[gene.names, ]$avg_log2FC)
# [1] 0.8242791


df.log2FC$bulk_logFC = NA
common.genes = intersect(de.bulk$Symbol, selected.genes)
df.log2FC$bulk_logFC[match(common.genes, selected.genes)] = de.bulk$logFC[match(common.genes, de.bulk$Symbol)]
bulk_logFC = df.log2FC$bulk_logFC
a = df.log2FC$F1599_log2FC
cor(a[!is.na(a+bulk_logFC)], bulk_logFC[!is.na(a+bulk_logFC)])
# [1] -0.1445868
a = df.log2FC$PGP1_log2FC
cor(a[!is.na(a+bulk_logFC)], bulk_logFC[!is.na(a+bulk_logFC)])
# [1] 0.03433324
a = df.log2FC$WTC11_log2FC
cor(a[!is.na(a+bulk_logFC)], bulk_logFC[!is.na(a+bulk_logFC)])
# [1] 0.007702086

# The pi1 analysis 
bulk.de.genes = de.bulk$Symbol[de.bulk$adj.P.Val < 0.05]
bulk.de.genes = bulk.de.genes[bulk.de.genes!= '--']
# 1397 genes
pi1.discover.bulk = unlist(lapply(list.de.mast, function(x){
  1-propTrueNull(x[intersect(bulk.de.genes, rownames(x)), 'p_val'])
}))
#1599      PGP1     WTC11 
#0.9922047 0.9867811 0.9995104 

pi1.discover.sc = unlist(lapply(list.de.mast, function(x){
  deg = rownames(x)[x$p_val_adj < 0.05]
  1-propTrueNull(de.bulk[match(intersect(deg, de.bulk$Symbol), de.bulk$Symbol), 'P.Value'])
}))
# 1599      PGP1     WTC11 
# 0.2639662 0.2604336 0.2618699 

####### Combine the logFC from all DE analysis together

df.log2FC = data.frame(genes = selected.genes, 
                       F1599_MAST = list.de.mast$`1599`[selected.genes, 'avg_log2FC'], 
                       PGP1_MAST = list.de.mast$PGP1[selected.genes, 'avg_log2FC'], 
                       WTC11_MAST = list.de.mast$WTC11[selected.genes, 'avg_log2FC'], 
                       F1599_Wilcox = list.de.wilcox$`1599`[selected.genes, 'avg_log2FC'], 
                       PGP1_Wilcox = list.de.wilcox$PGP1[selected.genes, 'avg_log2FC'], 
                       WTC11_Wilcox = list.de.wilcox$WTC11[selected.genes, 'avg_log2FC'])

df.log2FC$bulk_logFC = NA
common.genes = intersect(de.bulk$Symbol, selected.genes)
df.log2FC$bulk_logFC[match(common.genes, selected.genes)] = de.bulk$logFC[match(common.genes, de.bulk$Symbol)]

pdf('plots/scatter_logFC_mast.pdf', width = 8, height = 8)
plot(df.log2FC[, -1], pch = '.')
dev.off()

pairwise.cor = function(mat){
  n = ncol(mat)
  cname = colnames(mat)
  cor.mat = matrix(NA, nrow = n, ncol = n)
  diag(cor.mat) = 1
  colnames(cor.mat) <- rownames(cor.mat) <- cname
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      a = mat[, i]; b = mat[, j]
      a = a[!is.na(rowSums(mat[, c(i, j)]))]; b = b[!is.na(rowSums(mat[, c(i, j)]))]
      cor.mat[i, j] <- cor.mat[j, i] <- cor(a, b)
    }
  }
  return(cor.mat)
}
pairwise.cor(as.matrix(df.log2FC[, -1]))
#              F1599_MAST   PGP1_MAST  WTC11_MAST F1599_Wilcox PGP1_Wilcox WTC11_Wilcox   bulk_logFC
#F1599_MAST            1 -0.07273881 0.060571156   1.00000000 -0.07273881  0.060571156 -0.144586847
#PGP1_MAST    -0.07273881          1 0.824279120  -0.07273881  1.00000000  0.824279120  0.034333243
#WTC11_MAST    0.06057116  0.82427912          1   0.06057116  0.82427912  1.000000000  0.007702086
#F1599_Wilcox  1.00000000 -0.07273881 0.060571156           1 -0.07273881  0.060571156 -0.144586847
#PGP1_Wilcox  -0.07273881  1.00000000 0.824279120  -0.07273881          1  0.824279120  0.034333243
#WTC11_Wilcox  0.06057116  0.82427912 1.000000000   0.06057116  0.82427912           1  0.007702086
#bulk_logFC   -0.14458685  0.03433324 0.007702086  -0.14458685  0.03433324  0.007702086           1






list.de.full = c(lapply(list.de.wilcox, function(x){
  rownames(x)[x$p_val_adj < 0.01 & abs(x$avg_log2FC) > 1]
}), 

lapply(list.de.mast, function(x){
  rownames(x)[x$p_val_adj < 0.01 & abs(x$avg_log2FC) > 1]
})
)

names(list.de.full) = c(paste0('Wilcox_', background), paste0('MAST_', background))
list.de.full$bulk = de.bulk$Symbol[de.bulk$adj.P.Val < 0.01 & abs(de.bulk$logFC) > 1]
list.de.full$bulk = list.de.full$bulk[list.de.full$bulk != '--']

# Wilcox_1599  Wilcox_PGP1 Wilcox_WTC11    MAST_1599    MAST_PGP1   MAST_WTC11    bulk 
#    557          267          320          558          267          321          300

library(UpSetR)
pdf('plots/UpSet_Wilcox.pdf', width = 8, height = 6)
upset(fromList(list.de.full[c(1:3)]), order.by = 'freq')
dev.off()

pdf('plots/UpSet_MAST.pdf', width = 8, height = 6)
upset(fromList(list.de.full[c(4:6)]), order.by = 'freq')
dev.off()

pdf('plots/UpSet_Wilcox_with_bulk.pdf', width = 8, height = 6)
upset(fromList(list.de.full[c(1:3, 7)]), order.by = 'freq')
dev.off()

pdf('plots/UpSet_MAST_with_bulk.pdf', width = 8, height = 6)
upset(fromList(list.de.full[c(4:7)]), order.by = 'freq')
dev.off()


# Down-regulated only for single-cell data
list.down.wilcox = lapply(list.de.wilcox, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC > 1]
})
list.down.mast = lapply(list.de.mast, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC > 1]
})
pdf('plots/UpSet_Wilcox_down.pdf', width = 8, height = 6)
upset(fromList(list.down.wilcox), order.by = 'freq')
dev.off()

pdf('plots/UpSet_MAST_down.pdf', width = 8, height = 6)
upset(fromList(list.down.mast), order.by = 'freq')
dev.off()
# [1] "AL161804.1" "LINC01829" 

# Up-regulated only for single-cell data
list.up.wilcox = lapply(list.de.wilcox, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC < -1]
})
list.up.mast = lapply(list.de.mast, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC < -1]
})
pdf('plots/UpSet_Wilcox_up.pdf', width = 8, height = 6)
upset(fromList(list.up.wilcox), order.by = 'freq')
dev.off()

pdf('plots/UpSet_MAST_up.pdf', width = 8, height = 6)
upset(fromList(list.up.mast), order.by = 'freq')
dev.off()
# [1] "SHISA6"     "KHDRBS2"    "HSPA1A"     "VCAN"       "KHDRBS2-OT" "SLIT2"      "GNGT1"      "BRINP2"     "XDH"        "TRPC4"      "LINC02254" 
common.down = intersect(intersect(list.down.wilcox$`1599`, list.down.wilcox$PGP1), list.down.wilcox$WTC11)
common.up = intersect(intersect(list.up.wilcox$`1599`, list.up.wilcox$PGP1), list.up.wilcox$WTC11)
# We are interested in 1. the intersect between three backgrounds. 
a = VlnPlot(adnp.8.seu[, adnp.8.seu$background == '1599'], features = common.down, group.by = 'genotype') 
b = VlnPlot(adnp.8.seu[, adnp.8.seu$background == 'PGP1'], features = common.down, group.by = 'genotype') 
c = VlnPlot(adnp.8.seu[, adnp.8.seu$background == 'WTC11'], features = common.down, group.by = 'genotype') 

library(cowplot)
pdf('plots/Vln_plot_down_2_genes.pdf', width = 6, height = 9)
plot_grid(a, b, c, align = 'h', nrow = 3)
dev.off()

pdf('plots/Vln_plot_up_11_genes_1599.pdf', width = 9, height = 6)
VlnPlot(adnp.8.seu[, adnp.8.seu$background == '1599'], features = common.up, group.by = 'genotype') 
dev.off()
pdf('plots/Vln_plot_up_11_genes_PGP1.pdf', width = 9, height = 6)
VlnPlot(adnp.8.seu[, adnp.8.seu$background == 'PGP1'], features = common.up, group.by = 'genotype') 
dev.off()
pdf('plots/Vln_plot_up_11_genes_WTC11.pdf', width = 9, height = 6)
VlnPlot(adnp.8.seu[, adnp.8.seu$background == 'WTC11'], features = common.up, group.by = 'genotype') 
dev.off()

# Combine Het and Homo get together for 1599 and PGP1

# list.de.with.homo.wilcox = list()
for(s in background[1:2]){
  temp = adnp.sub[, adnp.sub$background == s]
  temp$gt = factor('WT', levels = c('WT', 'Mut'))
  temp$gt[temp$genotype != 'WT'] = 'Mut'
  
  de.wilcox = FindMarkers(object = temp, ident.1 = 'WT', ident.2 = 'Mut', features = selected.genes, group.by = 'gt', logfc.threshold = 0) 
  save(de.wilcox, file = paste0('Data_store/de_wilcox_with_Homo_', s, '.RData'))
}
load('Data_store/de_wilcox_with_Homo_1599.RData')
list.de.wilcox$Homo_1599 = de.wilcox
load('Data_store/de_wilcox_with_Homo_PGP1.RData')
list.de.wilcox$Homo_PGP1 = de.wilcox
save(list.de.wilcox, file = 'Data_store/list_de_seurat_wilcox_with_Homo.RData')
df.log2FC = data.frame(genes = selected.genes, 
                       F1599_log2FC = list.de.wilcox$`1599`[selected.genes, 'avg_log2FC'], 
                       F1599_log2FC_Homo = list.de.wilcox$Homo_1599[selected.genes, 'avg_log2FC'], 
                       PGP1_log2FC = list.de.wilcox$PGP1[selected.genes, 'avg_log2FC'], 
                       PGP1_log2FC_Homo = list.de.wilcox$Homo_PGP1[selected.genes, 'avg_log2FC'], 
                       WTC11_log2FC = list.de.wilcox$WTC11[selected.genes, 'avg_log2FC'])

pdf('plots/scatter_logFC_wilcox_with_homo.pdf', width = 8, height = 8)
plot(df.log2FC[, -1], pch = '.')
dev.off()
cor.wilcox = pairwise.cor(df.log2FC[, -1])

#                   F1599_log2FC F1599_log2FC_Homo PGP1_log2FC PGP1_log2FC_Homo WTC11_log2FC
#F1599_log2FC        1.00000000       0.925989312 -0.07273881      -0.05134056  0.060571156
#F1599_log2FC_Homo   0.92598931       1.000000000 -0.06450400      -0.02599503  0.004535334
#PGP1_log2FC        -0.07273881      -0.064503997  1.00000000       0.97129431  0.824279120
#PGP1_log2FC_Homo   -0.05134056      -0.025995034  0.97129431       1.00000000  0.858106523
#WTC11_log2FC        0.06057116       0.004535334  0.82427912       0.85810652  1.000000000

list.de.full$Wilcox_1599_Homo = rownames(list.de.wilcox$Homo_1599)[list.de.wilcox$Homo_1599$p_val_adj < 0.01 & list.de.wilcox$Homo_1599$avg_log2FC < -1]
list.de.full$Wilcox_PGP1_Homo = rownames(list.de.wilcox$Homo_PGP1)[list.de.wilcox$Homo_PGP1$p_val_adj < 0.01 & list.de.wilcox$Homo_PGP1$avg_log2FC < -1]

upset(fromList(list.de.full[grep('Wilcox', names(list.de.full))]), order.by = 'freq')
# 7 common genes 

list.down.wilcox$Homo_1599 = rownames(list.de.wilcox$Homo_1599)[list.de.wilcox$Homo_1599$p_val_adj < 0.01 & list.de.wilcox$Homo_1599$avg_log2FC > 1]
list.down.wilcox$Homo_PGP1 = rownames(list.de.wilcox$Homo_PGP1)[list.de.wilcox$Homo_PGP1$p_val_adj < 0.01 & list.de.wilcox$Homo_PGP1$avg_log2FC > 1]

list.up.wilcox$Homo_1599 = rownames(list.de.wilcox$Homo_1599)[list.de.wilcox$Homo_1599$p_val_adj < 0.01 & list.de.wilcox$Homo_1599$avg_log2FC < -1]
list.up.wilcox$Homo_PGP1 = rownames(list.de.wilcox$Homo_PGP1)[list.de.wilcox$Homo_PGP1$p_val_adj < 0.01 & list.de.wilcox$Homo_PGP1$avg_log2FC < -1]


list.up.wilcox = lapply(list.de.wilcox, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC < -1]
})
list.up.mast = lapply(list.de.mast, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC < -1]
})
pdf('plots/UpSet_Wilcox_up_Homo.pdf', width = 8, height = 6)
upset(fromList(list.up.wilcox), order.by = 'freq')
dev.off()

pdf('plots/UpSet_Wilcox_down_Homo.pdf', width = 8, height = 6)
upset(fromList(list.down.wilcox), order.by = 'freq')
dev.off()
# [1] "SHISA6"     "KHDRBS2"    "HSPA1A"     "VCAN"       "KHDRBS2-OT" "SLIT2"      "GNGT1"      "BRINP2"     "XDH"        "TRPC4"      "LINC02254" 
names( which(table(unlist(list.down.wilcox)) == 5) )
names( which(table(unlist(list.up.wilcox)) == 5) )
# [1] "HSPA1A"     "KHDRBS2"    "KHDRBS2-OT" "SHISA6"     "SLIT2"      "TRPC4"      "VCAN" 

common.down = names( which(table(unlist(list.down.wilcox)) == 5) )
common.up = names( which(table(unlist(list.up.wilcox)) == 5) )


pdf('plots/Vln_plot_up_7_genes_1599.pdf', width = 8, height = 9)
VlnPlot(adnp.8.seu[, adnp.8.seu$background == '1599'], features = common.up, group.by = 'genotype') 
dev.off()
pdf('plots/Vln_plot_up_7_genes_PGP1.pdf', width = 8, height = 9)
VlnPlot(adnp.8.seu[, adnp.8.seu$background == 'PGP1'], features = common.up, group.by = 'genotype') 
dev.off()
pdf('plots/Vln_plot_up_7_genes_WTC11.pdf', width = 8, height = 9)
VlnPlot(adnp.8.seu[, adnp.8.seu$background == 'WTC11'], features = common.up, group.by = 'genotype') 
dev.off()

FeaturePlot(adnp.8.seu, features = common.up)

save(selected.genes, list.de.full, list.de.mast, list.de.wilcox, list.down.mast, 
     list.up.mast, list.up.wilcox,file = 'Data_store/list_de_ensemble.RData')


### load and make figures ###### 
#
df.log2FC.mast = data.frame(do.call('cbind', lapply(list.de.mast, function(x){
  temp = rep(NA, length(selected.genes))
  names(temp) = selected.genes
  temp[rownames(x)] = x$avg_log2FC
  return(temp)
})) ) 
names(df.log2FC.mast) = c('F1599', 'PGP1', 'WTC11', 'Homo_PGP1', 'Homo_1599')
rownames(df.log2FC.mast) = selected.genes

pdf('plots/scatter_logFC_mast_with_Homo.pdf', width = 8, height = 8)
plot(df.log2FC.mast,  pch = '.', main = 'MAST log2FC')
dev.off()

df.log2FC.wilcox = data.frame(do.call('cbind', lapply(list.de.wilcox, function(x){
  temp = rep(NA, length(selected.genes))
  names(temp) = selected.genes
  temp[rownames(x)] = x$avg_log2FC
  return(temp)
})) ) 
names(df.log2FC.wilcox) = c('F1599', 'PGP1', 'WTC11', 'Homo_1599', 'Homo_PGP1')
rownames(df.log2FC.wilcox) = selected.genes
df.log2FC.wilcox = df.log2FC.wilcox[, c('F1599', 'PGP1', 'WTC11', 'Homo_PGP1', 'Homo_1599')]

pdf('plots/scatter_logFC_wilcox_with_Homo.pdf', width = 8, height = 8)
plot(df.log2FC.wilcox,  pch = '.', main = 'Wilcox log2FC')
dev.off()

list.up.wilcox = lapply(list.de.wilcox, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC < -1]
})
list.up.mast = lapply(list.de.mast, function(x){
  rownames(x)[x$p_val_adj < 0.01 & x$avg_log2FC < -1]
})

# We have the same results for MAST and Wilcox
common.down = names(which(table(unlist(list.down.wilcox)) == 5))
common.up = names(which(table(unlist(list.up.wilcox)) == 5))


df.de.genes = do.call('cbind', lapply(list.de.full, function(x){
  if(length(x) < 588){
    x = c(x, rep('', 588 - length(x)))
  }
  return(x)
}))
write.csv(df.de.genes, quote = FALSE, file = 'Data_store/DE_genes_all.csv')

# With only the male samples, and ignore the batch effect
adnp.sub = readRDS('Data_store/adnp_all_seu_light.rds')
adnp.male = adnp.sub[, adnp.sub$background!='1599' & adnp.sub$orig.ident != 'WTC11-Homo2']

# Get the DE genes from this dataset 
selected.genes = rownames(adnp.male)[rowMeans(adnp.male@assays$RNA$counts > 0) > 0.05]

# Differential expression
adnp.male$gt = factor('WT', levels = c('WT', 'Mut'))
adnp.male$gt[adnp.male$genotype != 'WT'] = 'Mut'
de.mast.male.Homo = FindMarkers(object = adnp.male, ident.1 = 'WT', ident.2 = 'Mut', test.use = 'MAST', features = selected.genes, group.by = 'gt', logfc.threshold = 0) 
de.mast.male = FindMarkers(object = adnp.male, ident.1 = 'WT', ident.2 = 'Het', test.use = 'MAST', features = selected.genes, group.by = 'genotype', logfc.threshold = 0)
de.wilcox.male.Homo = FindMarkers(object = adnp.male, ident.1 = 'WT', ident.2 = 'Mut', features = selected.genes, group.by = 'gt', logfc.threshold = 0) 
de.wilcox.male = FindMarkers(object = adnp.male, ident.1 = 'WT', ident.2 = 'Het', features = selected.genes, group.by = 'genotype', logfc.threshold = 0)

write.table(selected.genes, file = 'male_genes_background_for_de.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
save(de.mast.male, de.mast.male.Homo, de.wilcox.male, de.wilcox.male.Homo, file = 'Data_store/de_mast_male_comb.RData')

down = rownames(de.mast.male)[de.mast.male$p_val_adj < 0.05 & de.mast.male$avg_log2FC > 0.5]
up = rownames(de.mast.male)[de.mast.male$p_val_adj < 0.05 & de.mast.male$avg_log2FC < -0.5]

write.table(down, file = 'Data_store/down_de_male_mast.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(up, file = 'Data_store/up_de_male_mast.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

down = rownames(de.mast.male.Homo)[de.mast.male.Homo$p_val_adj < 0.05 & de.mast.male.Homo$avg_log2FC > 0.5]
up = rownames(de.mast.male.Homo)[de.mast.male.Homo$p_val_adj < 0.05 & de.mast.male.Homo$avg_log2FC < -0.5]
write.table(down, file = 'Data_store/down_de_male_mast_homo.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(up, file = 'Data_store/up_de_male_mast_homo.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)


down = rownames(de.wilcox.male)[de.wilcox.male$p_val_adj < 0.05 & de.wilcox.male$avg_log2FC > 0.5]
up = rownames(de.wilcox.male)[de.wilcox.male$p_val_adj < 0.05 & de.wilcox.male$avg_log2FC < -0.5]
write.table(down, file = 'Data_store/down_de_male_wilcox.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(up, file = 'Data_store/up_de_male_wilcox.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

down = rownames(de.wilcox.male.Homo)[de.wilcox.male.Homo$p_val_adj < 0.05 & de.wilcox.male.Homo$avg_log2FC > 0.5]
up = rownames(de.wilcox.male.Homo)[de.wilcox.male.Homo$p_val_adj < 0.05 & de.wilcox.male.Homo$avg_log2FC < -0.5]
write.table(down, file = 'Data_store/down_de_male_wilcox_homo.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(up, file = 'Data_store/up_de_male_wilcox_homo.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

