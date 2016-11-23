
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