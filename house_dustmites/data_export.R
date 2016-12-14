load('preprocessed.R')

cell_type <- "T-Cell"
tissue <- "Blood"
export_cols <- c("Row.names")

export_filter <- pData(sc_qc)$cell_type == cell_type & pData(sc_qc)$tissue <- tissue

write.table(
  pData(sc_qc)[export_filter, export_cols],
  "data/qc_blister_gd_files.csv", sep = ",", 
  quote = F, col.names= F, row.names = F)
