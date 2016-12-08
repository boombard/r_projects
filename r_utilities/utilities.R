
# Replace ensemble gene IDs with symbols

logVarMinus <- function (x) {
  return(mean(var(as.numeric(x)) - 1) + 1)
}

logMeanMinus <- function (x) {
  return(log(mean(exp(as.numeric(x)) - 1) + 1))
}

plotGeneDispersion <- function (logExpression) {
  gene_dispersion <- apply(logExpression, 1, disperse)
  gene_exprs_mean <- apply(logExpression, 1, mean)
  plot(gene_exprs_mean, gene_dispersion)
}

replace_ensembl_gene <- function(df, ensembl_map, axis = 1, drop_duplicates = F) {
  ensembl_names <- character()
  if (axis == 1) {
    ensembl_names <- rownames(df)
  }
  else if (axis == 2) {
    ensembl_names <- colnames(df)
  }
  symbols <- ensembl_map[match(ensembl_names, ensembl_map$Ensembl.Gene.ID),
                        'Associated.Gene.Name']
  new_names <- ifelse(is.na(symbols), ensembl_names, as.character(symbols)) 
  duplicate_index <- duplicated(new_names)
 
  if (axis == 1) {
    if (drop_duplicates) {
      df <- df[!duplicate_index, ]
      rownames(df) <- new_names[!duplicate_index]
    }
    else {
      rownames(df) <- new_names
    }
  }
  else if (axis == 2) {
    if (drop_duplicates) {
      df <- df[, !duplicate_index]
      colnames(df) <- new_names[!duplicate_index]
    }
    else {
      colnames(df) <- new_names
    }
  }
  
  return(df)
}

runSeurat <- function(raw_data, disp_cutoff = 0.5, pc_use = 20, 
                      interactive = FALSE, colour_category = NULL,
                      shape_category = NULL) {
  set.seed(0)
  # Use Seurat to find highly variable genes
  pbmc <- new("seurat", raw.data = raw_data)
  # pbmc <- Setup(pbmc, min.cells = 3, min.genes = 1600, do.logNormalize = F, total.expr = 1e4, project = "house_dustmite")
  pbmc <- Setup(pbmc, do.logNormalize = T, project = "house_dustmite")
  
  pbmc <- MeanVarPlot(pbmc,
                      fxn.x = expMean,
                      fxn.y = logVarDivMean,
                      x.low.cutoff = 0.015,
                      y.cutoff = disp_cutoff,
                      do.contour = F)
  print(paste0("Variable genes detected: ", length(pbmc@var.genes)))
  
  if (interactive) {
    disp_cutoff <- readline(prompt = "Y cutoff: ")
    disp_cutoff <- as.numeric(disp_cutoff)
  }
  
  pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes)
  VizPCA(pbmc, 1:2)
  PCAPlot(pbmc, 1, 2)
  PCElbowPlot(pbmc, num.pc = pc_use)
  
  if (interactive) {
    pc_use <- readline(prompt = "Number of PCs to use: ")
    pc_use <- as.integer(pc_use)
  }
  
  # PCHeatmap(pbmc, pc.use = 1:15, cells.use = 500, do.balanced = T, label.columns = F, use.full = F)
  # pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = F)
  pbmc <- FindClusters(pbmc, pc.use = 1:pc_use, resolution = 0.8, print.output = 0, save.SNN = T)
  pbmc <- RunTSNE(pbmc, dims.use = 1:pc_use, do.fast = T, perplexity = 15)
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = T, min.pct = 0.25, thresh.use = 0.25)
  
  pbmc.markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
  
  # pdf("plots/seurat_cluster_heatmap.pdf", width = 10, height = 15)
  # DoHeatmap(pbmc, genes.use = top20$gene, order.by.ident = T, slim.col.label = T, remove.key = T)
  # dev.off()
  
  TSNEPlot(pbmc)
  
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
  
  
  return(list("cellsObject" = pbmc, "clusterMarkers" = pmbc.markers))
}

plotLabelledClusters <- function(xData, yData,
                                 clusterData, colour = NULL, shape = NULL) {
  plotData <- data.frame(xData, yData, clusterData)
  centers <- plotData %>% group_by(clusterData) %>% 
    summarise(x = median(xData),  y = median(yData))
  centers <- data.frame(centers)
  
  if (is.null(colour)) {
    colour <- clusterData
  }
  print(length(colour))
  print(length(shape))
  print(length(plotData))
  
  # pData()[rownames(pData(sc_qc)) %in% attributes(pbmc@ident[pbmc@ident == 1])$names, "tissue"]
  
  # pdf("plots/variable_genes_tsne.pdf")
  p <- ggplot(data = data.frame(xData, yData), aes(x = xData, y = yData)) +
       geom_point(aes(colour = colour, shape = shape))
  p2 <- p + 
        geom_point(data = centers, aes(x = x, y = y), size = 0, alpha = 0) +
        geom_text(data = centers, aes(x = x, y = y, label=clusterData), 
                  size = 8, fontface = "bold", alpha = 0.6)
  print(p2)
}

seuratDEHeatmap <- function(cellsObject, markers, filename = NULL) {
    
  markers %>% group_by(cluster) %>% top_n(20, avg_diff) -> top20
  
  if (!is.null(filename)) {
    pdf(filename, width = 10, height = 15)
  }
  DoHeatmap(cellsObject, genes.use = top20$gene, order.by.ident = T, 
            slim.col.label = T, remove.key = T)
   
  if (!is.null(filename)) {
    dev.off()
  }
}

plotHeatMap <- function(data, ..., ordering = "NA") {
  pallete <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
  
  heatmap.2(
    data,
    density.info = "none",
    trace = "none",
    col = topo.colors(75),
    dendrogram = "row",
    colsep = ordering
  )
}
