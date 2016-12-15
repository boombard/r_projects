library(Seurat)
source('../r_utilities/utilities.R')

rm(list=ls())
load('preprocessed.RData')

membrane_genes <- membrane_db$Gene_ID
cellFilter <- pData(sc_qc)$tissue == "Blister"

cellFilter <- logical(length=dim(pData(sc_qc))[1]) == F

set.seed(0)
# Use Seurat to find highly variable genes
pbmc <- new("seurat", raw.data = exprs(sc_qc[endog_genes, cellFilter]))
# pbmc <- Setup(pbmc, min.cells = 3, min.genes = 1600, do.logNormalize = F, total.expr = 1e4, project = "house_dustmite")
pbmc <- Setup(pbmc, do.logNormalize = T, project = "house_dustmite")

pbmc <- MeanVarPlot(pbmc,
                    fxn.x = expMean,
                    fxn.y = logVarDivMean,
                    x.low.cutoff = 0.015,
                    y.cutoff = 0.5,
                    do.contour = F)
print(paste0("Variable genes detected: ", length(pbmc@var.genes)))

pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes)
# VizPCA(pbmc, 1:2)
# PCAPlot(pbmc, 1, 2)

PCElbowPlot(pbmc, num.pc = 20)

pc_use = 15

# PCHeatmap(pbmc, pc.use = 1:15, cells.use = 500, do.balanced = T, label.columns = F, use.full = F)
# pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = F)
pbmc <- FindClusters(pbmc, pc.use = 1:pc_use, resolution = 1., print.output = 0, save.SNN = T)
pbmc <- RunTSNE(pbmc, dims.use = 1:pc_use, do.fast = T, perplexity = 30)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = T, min.pct = 0.25, thresh.use = 0.25)

pbmc.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20

# pdf("plots/seurat_cluster_heatmap.pdf", width = 10, height = 15)
dev.new()
DoHeatmap(pbmc, genes.use = top20$gene, order.by.ident = T, slim.col.label = T, remove.key = T)
# dev.off()

dev.new()
TSNEPlot(pbmc)

dev.new()
plotLabelledClusters(pbmc@tsne.rot$tSNE_1,
                     pbmc@tsne.rot$tSNE_2,
                     pbmc@ident,
                     colour = pData(sc_qc)$cell_type,
                     shape = pData(sc_qc)$tissue)

# tsne_seurat <- pbmc@tsne.rot
# tsne_seurat$ident <- pbmc@ident
# centers <- tsne_seurat %>% group_by(ident) %>% summarise(x = median(tSNE_1), 
#                                                          y = median(tSNE_2))
# centers <- data.frame(centers)
# 
# pData(sc_qc)[rownames(pData(sc_qc)) %in% attributes(pbmc@ident[pbmc@ident == 1])$names, "tissue"]
# 
# p <- ggplot(data = pbmc@tsne.rot, aes(x=tSNE_1, y=tSNE_2)) +
#      geom_point(aes(colour = colour_category, shape = shape_category))
# p2 <- p + geom_point(data = centers, aes(x = x, y = y), size = 0, alpha = 0) +
#       geom_text(data = centers, aes(x = x, y = y, label=ident), size = 8, fontface = "bold", alpha = 0.6)
# print(p2)

# pdf("plots/tsne_all_cluster.pdf")

# dev.off()

cd_genes <- c("CD4", "CD3G", "CD8A", "CD8B")
cd45 <- c("PTPRC")
central_memory <- c("CCR7", "SELL")
sp_dp_genes <- c("KIT", "IL7R", "PTGDR2", "CCR7", "SELL")
cytokine <- c("IFNG", "IL2", "IL10", "GZMB")
b_cell <- c("CD19")
gene_names[grep("IL13", gene_names)]

FeaturePlot(pbmc, cd45,
            cols.use = c("green", "blue"))

# seuratDEHeatmap(pbmc, pbmc.markers,
#                 filename = "plots/cluster_de_heatmap_blood.pdf")


log_expression <- log10(t(sc_qc[endog_genes, ]@assayData$tpm + 1))

var_membrane_genes <- pbmc@var.genes[
  pbmc@var.genes %in% membrane_genes]
length(var_membrane_genes)
  
var_membrane_filter <- colnames(log_expression) %in% var_membrane_genes
length(var_membrane_filter)

# Get the frequency of the expression within the cells
membrane_expression <- pbmc@raw.data[var_membrane_genes, ]
expression_distribution <- rowSums(membrane_expression > 0)
frequent_membrane_genes <- names(expression_distribution > 20)

heatmap(t(log_expression[, var_membrane_filter]))
source('../r_utilities/utilities.R')

cells_ident <- factor(pbmc@ident)
table(cells_ident)
cells_use <- pbmc@cell.names[order(cells_ident)]
colsep_use <- cumsum(table(cells_ident))
col_lab <- rep("", length(cells_use))
col_lab[round(colsep_use - table(cells_ident) / 2) + 1] = levels(cells_ident)
cex_col <- 0.2 + 1/log10(length(unique(cells_ident)))

# pdf("plots/membrane_gene_cluster_heatmap.pdf", width = 10, height = 25)
dev.new()
heatmap.2(
  minmax(pbmc@scale.data[frequent_membrane_genes, ],
         min = -2.5, max = 2.5),
  Rowv = T,
  Colv = NA,
  trace = "none",
  colsep = colsep_use,
  labCol = col_lab,
  cexCol = cex_col,
  dendrogram = "row",
  col = pyCols,
  keysize = 0.5,
  key.title = NA)
# dev.off()
