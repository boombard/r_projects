library(scater, quietly = TRUE)
library(dplyr)
library(data.table)
library(ggplot2)
library(gplots)
library(GGally)
options(stringsAsFactors = FALSE)

tpm_data <- data.frame(
  fread('data/tpm_t.csv'), row.names = 1)
qc <- data.frame(
  fread('data/salmon_qc.csv'), row.names = 1)
metadata <- data.frame(
  fread('data/technical_metadata.csv'), 
  row.names = 1)
mt_genes <- read.csv(
  '../annotations/human_mt_genes.tsv',
  header = TRUE,  sep = '\t')
symbol_map <- read.table(
  "../annotations/human_annotation_map.tsv", 
  header = TRUE, sep = '\t')
membrane_db <- read.table(
  "../annotations/plasma_membrane_rvt.csv", header = T, sep = ",")

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

# pdf(file = "~/r_projects/house_dustmites/quality_pairplots.pdf", width = 10, height = 10)
# ggscatmat(
#   sceset@phenoData@data,
#   columns = c("num_processed", "total_features", 
#               "percent_mapped", "pct_tpm_feature_controls_MT"),
#   color = "sample_type"
# )
# dev.off()

# hist(sceset$total_features, breaks = 100)
# abline(v = 1600, col = "red")

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

# pdf("~/r_projects/house_dustmites/quality_filter_pairplots.pdf", 
#     height = 10,
#     width = 10)
# ggscatmat(
#   sceset@phenoData@data,
#   columns = c("num_processed", "total_features", 
#               "percent_mapped", "pct_tpm_feature_controls_MT"),
#   color = "quality_color"
# )
# dev.off()

# Around 60%
sum(sceset$use) / length(sceset$use) * 100

# Gene Filtering
# scater::plotQC(sceset, type = "highest-expression", exprs_values = "tpm")
# scater::plotQC(sc_qc, type = "exprs-freq-vs-mean")

# Remove undetectable genes
filter_genes <- apply(tpm(sceset[, pData(sceset)$use]), 1, 
                      function (x) length(x[x > 1]) > 2)
sum(filter_genes)
fData(sceset)$use <- filter_genes

dim(sceset[fData(sceset)$use, pData(sceset)$use])

# quality_summary <- aggregate(pData(sceset)$use,
#                              by = list("sample_type"=pData(sceset)[, "sample_type"]),
#                              FUN=mean)
# quality_summary <- merge(
#   quality_summary,
#   aggregate(pData(sceset)$use,
#             by = list("sample_type"=pData(sceset)[, "sample_type"]),
#             FUN=length),
#   by = "sample_type")
# colnames(quality_summary) <- c("sample_type", "pct_pass", "total")
# quality_summary
# write.table(quality_summary, file = "/Users/ge2/Desktop/quality_pass.csv", sep = ",")

# Apply the final QC filters
sc_qc <- sceset[fData(sceset)$use, pData(sceset)$use]
endog_genes <- !fData(sc_qc)$is_feature_control

save(sc_qc, endog_genes, symbol_map, membrane_db, file = "preprocessed.RData")
