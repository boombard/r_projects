
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