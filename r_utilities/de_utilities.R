ksDE <- function (expression, condition1, condition2, fdr = T) {
  pVals <- rep(1, nrow(expression))
  for (i in 1:nrow(expression)) {
    res <- ks.test(
      expression[i, condition1],
      expression[i, condition2]
    )
    pVals[i] <- res$p.value
  }
  if (fdr) {
    pVals <- p.adjust(pVals, method = "bonferroni")
  }
  pVal_df <- data.frame('pValues' = pVals,
                        row.names = rownames(expression))
  # rownames(pVal_df) <- rownames(expression)
  
  # Order the df
  ordered_pVals <- pVal_df[order(pVal_df$pValues), , drop = F]
  return(ordered_pVals);
}