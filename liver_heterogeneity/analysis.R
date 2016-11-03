library(scater)
library(SC3)
library(data.table)
library(gProfileR)
source('../r_utilities/utilities.R')
options(stringsAsFactors = FALSE)

# Load the data
tpm <- data.frame(fread('/Users/ge2/data/liver_heterogeneity/results/quality_tpm.csv'), row.names = 1)
qc <- data.frame(read.csv('/Users/ge2/data/liver_heterogeneity/qc/salmon_qc.csv'), row.names = 1)
metadata <- data.frame(read.csv('/Users/ge2/data/liver_heterogeneity/qc/salmon_qc.csv'), row.names = 1)
metrics <- data.frame(read.csv('/Users/ge2/data/liver_heterogeneity/results/metrics.csv'), row.names = 1)
symbol_map <- read.table("/Users/ge2/data/annotations/human_annotation_map.tsv", 
                         header = TRUE, sep = '\t')
cycle_genes <- read.table("/Users/ge2/data/annotations/cell_cycle_genes.tsv", 
                          header = FALSE, sep = '\t')
cycle_gene_symbols <- as.character(cycle_genes[, 3])

tpm <- replace_ensembl_gene(tpm, symbol_map, axis = 1)
colnames(tpm)

# Filter out QC samples not in the TPM table
qc <- qc[ rownames(qc) %in% rownames(tpm), ]
metrics <- metrics[ rownames(metrics) %in% rownames(tpm), ]

# Sort the rows of qc and tpm to ensure equality
qc <- qc[ order(row.names(qc)), ]
tpm <- tpm[ order(row.names(tpm)), ]
metrics <- metrics[ order(row.names(metrics)), ]

# Filter out genes that aren't expressed in any cell
filter_genes <- apply(tpm, 2, function (x) length(x[x > 1]) >= 2);
tpm <- tpm[, filter_genes]

# Create the scater object
qc_data <- new("AnnotatedDataFrame", data=qc)
metrics_data <- new("AnnotatedDataFrame", data=metrics)

sceset <- scater::newSCESet(
  tpmData = t(tpm),
  phenoData = metrics_data
)

# Plot PCA
scater::plotPCA(sceset, exprs_values="tpm", size_by="num_processed")
 
# Confounding factors for PCA
scater::plotQC(sceset, type='find-pcs', variable='read_no', exprs_values='tpm')
scater::plotQC(sceset, type='find-pcs', variable='mapping_rate', exprs_values='tpm')
scater::plotQC(sceset, type='find-pcs', variable='mt_pct', exprs_values='tpm')
scater::plotQC(sceset, type='find-pcs', variable='expressed_genes', exprs_values='tpm')

# Seems as though the confounding factors can't readily be explained by principal
# components.

scater::plotQC(sceset, type='expl', exprs_values='tpm',
               variables=c("mapping_rate", "read_no", "mt_pct", "expressed_genes"))

# Normal linear regression doesn't indicate confounding factors either

# Look at the SC3 clustering
is_exprs(sceset) <- exprs(sceset) > 1
sceset <- calculateQCMetrics(sceset)
sceset <- sc3(sceset, ks=2:4, exprs_values = 'tpm')

sc3_plot_consensus(sceset, k = 3)

# Look at the DE

sc3_plot_de_genes(sceset, k = 3)
sc3_summarise_results(sceset, k = 3)

de_genes <- sceset@sc3$biology$`3`$de.genes

profile <- gprofiler(rownames(de_genes)[1:100], organism = "hsapiens", 
                     ordered_query = T, correction_method = "fdr", hier_filtering = "moderate")

# Try out seurat

library(seurat)
liver <- new('seurat', raw.data = t(tpm))
liver <- Setup(liver, min.cells = 3, min.genes = 200, do.logNormalize = T, total.expr = 1e4, project = 'liver')
liver <- MeanVarPlot(liver, fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                    y.cutoff = 0.5, do.contour = F)