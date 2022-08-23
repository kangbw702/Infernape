
# test differential APA events between two cell groups

single_gene_test_1vsothers <- function(gene.to.test, counts, anno.tbl, attr.tbl, ctype.colname, base_grp, cut.low.pct = 0.05) {

  # find associated peaks
  peak.to.test = rownames(anno.tbl)[which(anno.tbl$gene == gene.to.test)]
  # order by mode.pos (small -> large)
  peak.to.test = peak.to.test[order(anno.tbl[peak.to.test, "mode.pos"])]
  strand = utils::head(anno.tbl[peak.to.test, 'strand'], 1)
  mode.pos = anno.tbl[peak.to.test, 'mode.pos']

  # generate df.test
  df.test = data.frame(cellid = colnames(counts), cl = attr.tbl[, ctype.colname], stringsAsFactors = FALSE)
  for (i in 1:length(peak.to.test)) df.test[, paste0("raw|", peak.to.test[i])] = counts[peak.to.test[i],]
  df.test[, paste0("raw|", gene.to.test)] = rowSums(df.test[, 3:(3+length(peak.to.test)-1)])



}
