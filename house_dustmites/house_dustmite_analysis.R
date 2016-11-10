library(scater, quietly = TRUE)
library(SC3)
library(data.table)
library(gProfileR)
library(ggplot2)
library(GGally)
source('../r_utilities/utilities.R')
options(stringsAsFactors = FALSE)

tpm <- data.frame(fread('/Users/ge2/data/house_dustmites/results/tpm_t.csv'),
                  row.names = 1)
qc <- data.frame(fread('/Users/ge2/data/house_dustmites/qc/salmon_qc.csv'),
                 row.names = 1)
metadata <- data.frame(fread('/Users/ge2/data/house_dustmites/technical_metadata.csv'),
                       row.names = 1)
mt_genes <- read.csv('/Users/ge2/data/annotations/human_mt_genes.tsv', 
                      header = TRUE, sep = '\t')
symbol_map <- read.table("/Users/ge2/data/annotations/human_annotation_map.tsv", 
                         header = TRUE, sep = '\t')

tpm <- replace_ensembl_gene(tpm, symbol_map, axis = 1)
rownames(tpm)

# Merge the qc and metadata files
annotation_df <- merge(metadata, qc, by = 0, all.x = TRUE, sort = TRUE)
rownames(annotation_df) <- annotation_df[,"sample_name"]

# Change the rownames of the tpm
sample_names <- annotation_df[match(rownames(tpm), annotation_df[,'Row.names']), "sample_name"]
rownames(tpm) <- sample_names

# Filter out QC samples not in the TPM table
annotation_df <- annotation_df[ rownames(annotation_df) %in% rownames(tpm), ]


# Sort the rows of qc and tpm to ensure equality
annotation_df <- annotation_df[ order(row.names(annotation_df)), ]
tpm <- tpm[ order(row.names(tpm)), ]


pheno_data <- new("AnnotatedDataFrame", annotation_df)
rownames(pheno_data) <- pheno_data$sample_name

sceset <- scater::newSCESet(
  tpmData = t(tpm),
  phenoData = pheno_data
)

# Controls
control_labels <- c("Bulk control", "Empty", "rnA control", "RNA control")
sceset$is_control <- sceset$cell_type %in% control_labels

# Get the ERCC and MT values
ercc <- featureNames(sceset)[grepl("ERCC.", featureNames(sceset))]
mt_list <- mt_genes$Ensembl.Gene.ID

# Define minimum expression limit
is_exprs(sceset) <- exprs(sceset) > 2

sceset <- scater::calculateQCMetrics(
  sceset,
  feature_controls = list(ERCC = ercc, MT = mt_list)
)

hist(sceset$num_processed, breaks = 100)
hist(sceset$percent_mapped, breaks = 100)

sceset@phenoData@data$pct_tpm_feature_controls_MT

ggscatmat(
  sceset@phenoData@data,
  columns = c("num_processed", "total_features", 
              "percent_mapped", "pct_tpm_feature_controls_MT"),
  color = "cell_type"
)

sceset@phenoData@data[, c('num_processed')]
