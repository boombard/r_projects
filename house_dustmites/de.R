library(limma)
library(pheatmap)
library(dplyr)
library(DESeq2)

blister_AB_data <- read.csv('tracer/blister_AB/cell_data.csv', 
                            header = 1, sep = ',')
blood_GD_data <- read.csv('tracer/blood_GD/cell_data.csv', 
                            header = 1, sep = ',')

rm(list=ls())
load('preprocessed.RData')

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

# KS test

log_expression <- log10(exprs(sc_qc) + 1)

phenoData <- pData(sc_qc)

# AB T-Cells
cd3_blood_filter <- pData(sc_qc)$cell_type == "CD3" & pData(sc_qc)$tissue == "Blood"
cd3_blister_filter <- pData(sc_qc)$cell_type == "CD3" & pData(sc_qc)$tissue == "Blister"

# Clonal filters for blister AB
blister_AB_clonal <- select(filter(blister_AB_data, group_size > 1),
                         cell_name)
blister_AB_clonal_filter <- is.element(phenoData$Row.names,
                                       blister_AB_clonal$cell_name)
blister_AB_nonclonal_filter <- cd3_blister_filter & !blister_AB_clonal_filter
phenoData$clonal <- is.element(phenoData$Row.names,
                               blister_AB_clonal$cell_name)

# GD blood cells
gd_blood_filter <- pData(sc_qc)$cell_type == "T-Cell" & pData(sc_qc)$tissue == "Blood"
gd_blister_filter <- pData(sc_qc)$cell_type == "T-Cell" & pData(sc_qc)$tissue == "Blister"

# Clonal filters for blood GD
blood_GD_clonal <- select(filter(blood_GD_data, group_size > 1),
                          cell_name)
blood_GD_clonal_filter <- is.element(phenoData$Row.names,
                                     blood_GD_clonal$cell_name)
blood_GD_nonclonal_filter <- gd_blood_filter & ! blood_GD_clonal_filter
phenoData$clonal <- phenoData$clonal | is.element(phenoData$Row.names, blood_GD_clonal$cell_name)

# SP and DP filters
spdp_filter <- pData(sc_qc)$cell_type == "SP_and_DP"
dn_filter <- pData(sc_qc)$cell_type == "DN"

unique(pData(sc_qc)$cell_type)

blister_AB_clonal_filter
blister_AB_nonclonal_filter

source('../r_utilities/de_utilities.R')
ordered_pVals <- ksDE(log_expression, blood_GD_clonal_filter,
                      blood_GD_nonclonal_filter, fdr = F)

View(ordered_pVals)

# pVals <- rep(1, nrow(exprs(sc_qc)))
# for (i in 1:nrow(exprs(sc_qc))) {
#   res <- ks.test(
#     log_expression[i, cd3_blood_filter],
#     log_expression[i, gd_blood_filter]
#   )
#   pVals[i] <- res$p.value
# }
# 
# pVals <- p.adjust(pVals, method = "bonferroni")
# 
# pVal_df <- data.frame(pVals, rownames(log_expression))
# rownames(pVal_df) <- rownames(log_expression)
# rownames(pVal_df)[1:5]
# 
# # Order the df
# ordered_pVals <- pVal_df[order(pVal_df$pVals), ]

write.table(ordered_pVals[ordered_pVals$pVals < 0.05, ], 
            file = "data/tcell_GD_AB_blood_de.txt",
            sep = ",", quote = F, row.names = F, col.names = F)

pheatmap(log_expression[top_de, cd3_blood_filter | cd3_blister_filter],
         cutree_cols = 2,
         show_rownames = F,
         annotaion_col = pData(sc_qc)["tissue", cd3_blood_filter | cd3_blister_filter])
