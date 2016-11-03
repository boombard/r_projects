
# Replace ensemble gene IDs with symbols

replace_ensembl_gene <- function(df, ensembl_map, axis = 0) {
  symbols <- ensembl_map[match(colnames(df), ensembl_map$Ensembl.Gene.ID),
                        'Associated.Gene.Name']
  new_names <- ifelse(is.na(symbols), colnames(df), as.character(symbols))  
  
  if (axis == 0) {
    rownames(df) <- new_names
  }
  else if (axis == 1) {
    colnames(df) <- new_names
  }
  
  return(df)
}