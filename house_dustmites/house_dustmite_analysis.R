library(scater, quietly = TRUE)
library(dplyr)
library(SC3)
library(data.table)
library(gProfileR)
library(ggplot2)
library(gplots)
library(GGally)
library(M3Drop)
library(DESeq2)
library(limma)
library(Seurat)
library(gProfileR)
options(stringsAsFactors = FALSE)

tpm_data <- data.frame(
  fread('/Users/ge2/data/house_dustmites/results/tpm_t.csv'), row.names = 1)
qc <- data.frame(
  fread('/Users/ge2/data/house_dustmites/qc/salmon_qc.csv'), row.names = 1)
metadata <- data.frame(
  fread('/Users/ge2/data/house_dustmites/technical_metadata.csv'), 
  row.names = 1)
mt_genes <- read.csv(
  '/Users/ge2/data/annotations/human_mt_genes.tsv',
  header = TRUE,  sep = '\t')
symbol_map <- read.table(
  "/Users/ge2/data/annotations/human_annotation_map.tsv", 
  header = TRUE, sep = '\t')
membrane_db <- read.table(
  "/Users/ge2/data/annotations/plasma_membrane_rvt.csv", header = T, sep = ",")

source('../r_utilities/utilities.R')

tpm_data <- replace_ensembl_gene(tpm_data, symbol_map, axis = 2, drop_duplicates = T)
gene_names <- colnames(tpm_data)
membrane_genes <- membrane_db$Gene_ID

# Merge the qc and metadata files
annotation_df <- merge(metadata, qc, by = 0, all.x = TRUE, sort = TRUE)
rownames(annotation_df) <- annotation_df[,"sample_name"]

# Change the rownames of the tpm
sample_names <- annotation_df[
  match(rownames(tpm_data), annotation_df[,'Row.names']), "sample_name"]
rownames(tpm_data) <- sample_names

# Filter out QC samples not in the TPM table
annotation_df <- annotation_df[rownames(annotation_df) %in% rownames(tpm_data), ]

# Sort the rows of qc and tpm to ensure equality
annotation_df <- annotation_df[ order(row.names(annotation_df)), ]
tpm_data <- tpm_data[ order(row.names(tpm_data)), ]

sample_types <- unique(annotation_df$sample_type)

# Control Wells
control_labels <- c("Bulk control", "Empty", "RNA")
annotation_df$is_control <- annotation_df$sample_type %in% control_labels

control_annotation <- annotation_df[annotation_df$is_control, ]
cell_annotation <- annotation_df[!annotation_df$is_control, ]

control_tpm <- tpm_data[annotation_df$is_control, ]
cell_tpm <- tpm_data[!annotation_df$is_control, ]

# Remove samples from other experiments
alien_samples <- c("ARS004 Blood ILC cells", 
                   "ARS005 Blood CD4 T cells", 
                   "ARS005 Blood ILC cells")
alien_indices = cell_annotation$sample_type %in% alien_samples
cell_annotation <- cell_annotation[!alien_indices, ]
cell_tpm <- cell_tpm[!alien_indices, ]

pheno_data <- new("AnnotatedDataFrame", cell_annotation)
rownames(pheno_data) <- pheno_data$sample_name

sceset <- scater::newSCESet(
  tpmData = t(cell_tpm),
  phenoData = pheno_data
)

# Get the ERCC and MT values
ercc <- featureNames(sceset)[grepl("ERCC.", featureNames(sceset))]
mt_list <- mt_genes$Ensembl.Gene.ID

# Define minimum expression limit
is_exprs(sceset) <- exprs(sceset) > 1

sceset <- scater::calculateQCMetrics(
  sceset,
  feature_controls = list(ERCC = ercc, MT = mt_list)
)

jpeg(filename = "~/r_projects/house_dustmites/quality_pairplots.jpeg", quality = 100, width = 1000, height = 1000)
pdf(file = "~/r_projects/house_dustmites/quality_pairplots.pdf", width = 10, height = 10)
ggscatmat(
  sceset@phenoData@data,
  columns = c("num_processed", "total_features", 
              "percent_mapped", "pct_tpm_feature_controls_MT"),
  color = "sample_type"
)
title("Test")
dev.off()

hist(sceset$total_features, breaks = 100)
abline(v = 1600, col = "red")

# Define filters
filter_total_features = sceset$total_features > 1600
filter_num_mapped = sceset$num_mapped > 150000
filter_pct_mt <- sceset$pct_tpm_feature_controls_MT < 20
filter_pct_mapped <- sceset$percent_mapped > 50

sceset$use <- (
  filter_total_features &
  filter_num_mapped &
  filter_pct_mt &
  filter_pct_mapped
)

sceset$quality_color <- ifelse(sceset$use, "accept", "reject")

# Default filter

sceset$use_default <- (
  !sceset$filter_on_total_features &
  !sceset$filter_on_total_counts &
  !sceset$filter_on_pct_tpm_feature_controls_ERCC &
  !sceset$filter_on_pct_tpm_feature_controls_MT
)

# Automatic PCA filter on quality features (Doesn't really work for TPM)
sceset <- scater::plotPCA(
  sceset, 
  size_by = "total_features",
  shape_by = "use",
  pca_data_input = "pdata",
  detect_outliers = TRUE,
  return_SCESet  = TRUE,
  exprs_values = "tpm",
  selected_variables = c("pct_exprs_top_100_features", 
                         "total_features", 
                         "pct_exprs_feature_controls",
                         "n_detected_feature_controls",
                         "log10_exprs_endogenous_features",
                         "log10_counts_feature_controls")
)

# Remove detected outliers
sceset$use <- sceset$use & !sceset$outlier

pdf("~/r_projects/house_dustmites/quality_filter_pairplots.pdf", 
    height = 10,
    width = 10)
ggscatmat(
  sceset@phenoData@data,
  columns = c("num_processed", "total_features", 
              "percent_mapped", "pct_tpm_feature_controls_MT"),
  color = "quality_color"
)
dev.off()

# Around 60%
sum(sceset$use) / length(sceset$use) * 100

# Gene Filtering
scater::plotQC(sceset, type = "highest-expression", exprs_values = "tpm")
scater::plotQC(sc_qc, type = "exprs-freq-vs-mean")

# Remove undetectable genes
filter_genes <- apply(tpm(sceset[, pData(sceset)$use]), 1, 
                      function (x) length(x[x > 1]) > 2)
sum(filter_genes)
fData(sceset)$use <- filter_genes

dim(sceset[fData(sceset)$use, pData(sceset)$use])

quality_summary <- aggregate(pData(sceset)$use,
                             by = list("sample_type"=pData(sceset)[, "sample_type"]),
                             FUN=mean)
quality_summary <- merge(
  quality_summary,
  aggregate(pData(sceset)$use,
            by = list("sample_type"=pData(sceset)[, "sample_type"]),
            FUN=length),
  by = "sample_type")
colnames(quality_summary) <- c("sample_type", "pct_pass", "total")
quality_summary
write.table(quality_summary, file = "/Users/ge2/Desktop/quality_pass.csv", sep = ",")

# Apply the final QC filters
sc_qc <- sceset[fData(sceset)$use, pData(sceset)$use]
endog_genes <- !fData(sc_qc)$is_feature_control

# First, log Transform
log_expression <- log10(t(sc_qc[endog_genes, ]@assayData$tpm + 1))
plotGeneDispersion(t(log_expression))

cell_type <- phenoData(sc_qc)$cell_type
tissue <- phenoData(sc_qc)$tissue

seuratFilter <- sc_qc$tissue == "Blister"
seuratFilter <- logical(length = length(sc_qc$tissue)) == F

results <- runSeurat(
  exprs(sc_qc[endog_genes, seuratFilter]), 
  interactive = T, colour_category = cell_type,
  shape_category = tissue, disp_cutoff = 0.5, pc_use = 15
)

pdf("plots/tsne_all_cluster.pdf")
plotLabelledClusters(results$cellsObject@tsne.rot$tSNE_1,
                     results$cellsObject@tsne.rot$tSNE_2,
                     results$cellsObject@ident,
                     colour = pData(sc_qc[endog_genes, seuratFilter])$cell_type,
                     shape = pData(sc_qc[endog_genes, seuratFilter])$tissue)
dev.off()

cd_genes <- c("CD4", "CD3G", "CD8A", "CD8B")
cd45 <- c("PTPRC")
central_memory <- c("CCR7", "SELL")
sp_dp_genes <- c("KIT", "IL7R", "PTGDR2", "CCR7", "SELL")
cytokine <- c("IFNG", "IL2", "IL10", "GZMB")
b_cell <- c("CD19")
gene_names[grep("IL13", gene_names)]

FeaturePlot(results$cellsObject, b_cell,
            cols.use = c("green", "blue"))

seuratDEHeatmap(results$cellsObject, results$clusterMarkers,
                filename = "plots/cluster_de_heatmap_blood.pdf")

var_membrane_genes <- results$cellsObject@var.genes[
  results$cellsObject@var.genes %in% membrane_genes]
length(var_membrane_genes)
  
var_membrane_filter <- colnames(log_expression) %in% var_membrane_genes
length(var_membrane_filter)

# Get the frequency of the expression within the cells
membrane_expression <- results$cellsObject@raw.data[var_membrane_genes, ]
expression_distribution <- rowSums(membrane_expression > 0)
frequent_membrane_genes <- names(expression_distribution > 20)

heatmap(t(log_expression[, var_membrane_filter]))
source('../r_utilities/utilities.R')

cells_ident <- factor(results$cellsObject@ident)
table(cells_ident)
cells_use <- results$cellsObject@cell.names[order(cells_ident)]
colsep_use <- cumsum(table(cells_ident))
col_lab <- rep("", length(cells_use))
col_lab[round(colsep_use - table(cells_ident) / 2) + 1] = levels(cells_ident)
cex_col <- 0.2 + 1/log10(length(unique(cells_ident)))

pdf("plots/membrane_gene_cluster_heatmap.pdf", width = 10, height = 25)
heatmap.2(
  minmax(results$cellsObject@scale.data[frequent_membrane_genes, ],
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
dev.off()

# PCA analysis
sc_pca <- prcomp(log_expression)
plot(sc_pca, type = "l")
pct_explained <- sc_pca$sdev / sum(sc_pca$sdev) * 100
plot(cumsum(pct_explained), type = "l")

# PCA models the data poorly

# Let's take the 1000 most highly expressed genes
top_gene_n <- 500

total_feature_exprs <- apply(log_expression, 2, mean)
total_exprs <- sum(log_expression)

gene_order <- order(total_feature_exprs, decreasing = TRUE)
percent_expression <- 100 * sum(total_feature_exprs[gene_order][1:top_gene_n]) / total_exprs
top_genes <- colnames(log_expression)[gene_order][1:top_gene_n]

pdf("~/r_projects/house_dustmites/pca.pdf")
scater::plotPCA(sc_qc[endog_genes, ],
               colour_by = "cell_type",
               size_by = "total_features",
               shape_by = "tissue",
               exprs_values = "tpm")
dev.off()

unique(pData(sc_qc)$cell_type)

pdf("~/r_projects/house_dustmites/tsne.pdf")
# , pData(sc_qc)$cell_type == "SP_and_DP" & pData(sc_qc)$tissue == "Blood"
scater::plotTSNE(sc_qc[top_genes, ],
                 ntop = 2000,
                 perplexity = 20,
                 colour_by = "cell_type",
                 shape_by = "tissue",
                 size_by = "total_features",
                 exprs_values = "tpm")
dev.off()

# DE analysis
m3_normalised_data <- M3Drop::M3DropCleanData(
  sceset@assayData$tpm,
  is.counts = FALSE, 
  min_detected_genes = 1600)

m3_model <- M3Drop::M3DropDropoutModels(m3_normalised_data$data)
DE_genes <- M3Drop::M3DropDifferentialExpression(
  sceset@assayData$tpm,
  mt_method = "fdr",
  mt_threshold = 0.01
)

# DESeq2
browseVignettes("DESeq2")

# Does not detect any DE genes

pdf("~/r_projects/house_dustmites/confounding_factor_analysis.pdf")
scater::plotQC(sc_qc[endog_genes, ],
               type = "expl",
               exprs_values = "tpm",
               variables = c("total_features",
                             "num_processed",
                             "percent_mapped",
                             "pct_feature_controls_MT"))
dev.off()

# Replace the ENSEMBL codes
log_expression <- replace_ensembl_gene(log_expression, symbol_map, axis = 1)

# Linear models for each gene
cd3_filter <- pData(sc_qc)$cell_type == "CD3"
tcell_filter <- pData(sc_qc)$cell_type == "T-Cell"

cd3_pdata <- pData(sc_qc)[cd3_filter, ]
tcell_pdata <- pData(sc_qc)[tcell_filter, ]

# DE between blood CD3 and blister CD3
tissue_factor <- as.factor(cd3_pdata$tissue)
design <- model.matrix(~ 0 + tissue_factor)
colnames(design) <- gsub("tissue_factor", "", colnames(design))

dim(log_expression[cd3_filter, ])
colnames(design)
contrast_matrix <- makeContrasts(Blood - Blister, levels = design)

fit <- lmFit(t(log_expression[cd3_filter, ]), design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef = 1, adjust = "BH")

results <- decideTests(fit2)
vennDiagram(results)
# colnames(design) <- unique(sc_qc$tissue)
length(results)
gprofiler(results$ID[1:100])
