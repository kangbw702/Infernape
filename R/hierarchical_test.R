
# test differential APA events between two cell groups

single_gene_test <- function(gene.to.test, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct) {

  # find associated peaks
  peak.to.test = rownames(anno.tbl)[which(anno.tbl$gene == gene.to.test)]
  # order by mode.pos (small -> large)
  peak.to.test = peak.to.test[order(anno.tbl[peak.to.test, "mode.pos"])]
  strand = utils::head(anno.tbl[peak.to.test, 'strand'], 1)
  mode.pos = anno.tbl[peak.to.test, 'mode.pos']

  # generate df.test
  df.test = data.frame(cellid = colnames(cnt), cl = attr.tbl[, ctype.colname], stringsAsFactors = FALSE)
  for (i in 1:length(peak.to.test)) df.test[, paste0("raw|", peak.to.test[i])] = cnt[peak.to.test[i], ]
  df.test[, paste0("raw|", gene.to.test)] = rowSums(df.test[, 3:(3 + length(peak.to.test) - 1)])

  # generate aggregated count/percentage table for chi-square test
  cnt.tbl = pct.tbl = matrix(NA, ncol = 2, nrow = length(peak.to.test))
  rownames(cnt.tbl) = rownames(pct.tbl) = peak.to.test
  if ( is.null(alt_grp)) colnames(cnt.tbl) = colnames(pct.tbl) = c(base_grp, paste0('non-', base_grp))
  if (!is.null(alt_grp)) colnames(cnt.tbl) = colnames(pct.tbl) = c(base_grp, alt_grp)

  for (i in 1:length(peak.to.test)){

    df.test[, peak.to.test[i]] = df.test[, (2+i)] / df.test[, paste0("raw|", gene.to.test)]
    c1 = sum(df.test[which(df.test$cl == base_grp), 2+i])
    if ( is.null(alt_grp)) c2 = sum(df.test[which(df.test$cl != base_grp), 2+i])
    if (!is.null(alt_grp)) c2 = sum(df.test[which(df.test$cl ==  alt_grp), 2+i])
    cnt.tbl[i,] = c(c1, c2)
    p1 = mean(df.test[df.test$cl == base_grp, 2+i] > 0)
    if ( is.null(alt_grp)) p2 = mean(df.test[df.test$cl != base_grp, 2+i] > 0)
    if (!is.null(alt_grp)) p2 = mean(df.test[df.test$cl ==  alt_grp, 2+i] > 0)
    pct.tbl[i, ] = c(p1, p2)

  }

  if (all(colSums(cnt.tbl) != 0) && all(rowSums(cnt.tbl) != 0)) {

    # chisq test
    if (all(cnt.tbl == 0 )) pval.chisq = 1 else { chisq <- stats::chisq.test(cnt.tbl); pval.chisq = chisq$p.value }

    # multinomial-dirichlet LRT
    subdata.1 = as.matrix(df.test[df.test$cl == base_grp, 3:(3 + length(peak.to.test) - 1)])
    if ( is.null(alt_grp)) subdata.2 = as.matrix(df.test[df.test$cl != base_grp, 3:(3 + length(peak.to.test) - 1)])
    if (!is.null(alt_grp)) subdata.2 = as.matrix(df.test[df.test$cl ==  alt_grp, 3:(3 + length(peak.to.test) - 1)])

    alt.llk = c(f.logL(subdata.1), f.logL(subdata.2))
    names(alt.llk) <- colnames(cnt.tbl)

    if (all(!is.na(alt.llk))) {
      subdata.0 = as.matrix(df.test[, 3:(3 + length(peak.to.test) - 1)]) # under the null
      null.llk = f.logL(subdata.0)
      ll.alt <- sum(alt.llk, na.rm = T)
      ll.null = null.llk
      lr <- -2 * (ll.null - ll.alt)
      pval.DM <- stats::pchisq(lr, length(peak.to.test) * (2-1), lower.tail = FALSE)
    } else {
      pval.DM = NA
    }

    # rank-based WARM
    if (strand == '-') pos.rel = (-rank(mode.pos) + length(peak.to.test)) / (length(peak.to.test) - 1)
    if (strand == '+') pos.rel = ( rank(mode.pos) - 1) / (length(peak.to.test) - 1)
    warm.grp1  = sum(pos.rel * cnt.tbl[, base_grp] / sum(cnt.tbl[, base_grp]))
    if ( is.null(alt_grp)) warm.grp2  = sum(pos.rel * cnt.tbl[, paste0('non-', base_grp)] / sum(cnt.tbl[, paste0('non-', base_grp)]))
    if (!is.null(alt_grp)) warm.grp2  = sum(pos.rel * cnt.tbl[, alt_grp] / sum(cnt.tbl[, alt_grp]))
    delta.warm = warm.grp1 - warm.grp2

    # check simultaneous low pct in either cell-type tested
    low.pct.flag = apply(pct.tbl, 2, function(t) {all(t < cut.low.pct)})
    low.pct.flag = any(low.pct.flag)

    # MPRO
    comb = t(utils::combn(1:length(peak.to.test), 2))
    prp.tbl = t(t(cnt.tbl) / colSums(cnt.tbl))
    pd.vec = NULL

    for (i in 1:nrow(comb)) {
      this.prp = prp.tbl[comb[i, ], ]
      this.pd = this.prp[1, 1] + this.prp[2, 2] - this.prp[2, 1] - this.prp[1, 2]
      pd.vec = c(pd.vec, this.pd)
    }

    id.pd = which.max(abs(pd.vec))
    strand.multiplier = 2*(strand == '-') - 1
    max.abs.prpchg.id = ifelse(strand == '-', paste0(comb[id.pd, 1], ' vs ', comb[id.pd, 2]), paste0(comb[id.pd, 2], ' vs ', comb[id.pd, 1]))
    max.prpchg = pd.vec[id.pd] * strand.multiplier
    prp.dif.1 = ifelse(strand == '-', prp.tbl[comb[id.pd, 1], 1] - prp.tbl[comb[id.pd, 2], 1], prp.tbl[comb[id.pd, 2], 1] - prp.tbl[comb[id.pd, 1], 1])
    prp.dif.2 = ifelse(strand == '-', prp.tbl[comb[id.pd, 1], 2] - prp.tbl[comb[id.pd, 2], 2], prp.tbl[comb[id.pd, 2], 2] - prp.tbl[comb[id.pd, 1], 2])
    # prp.dif.1 - prp.dif.2 == max.prpchg

    res = list(cnt.tbl = cnt.tbl,
               pct.tbl = pct.tbl,
               pval.chisq = pval.chisq,
               pval.DM = pval.DM,
               warm.base = warm.grp1,
               warm.alt = warm.grp2,
               delta.warm = delta.warm,
               low.pct.flag = low.pct.flag,
               MPRO.id = max.abs.prpchg.id,
               MPRO = max.prpchg,
               prp.dif.base = prp.dif.1,
               prp.dif.alt = prp.dif.2,
               pks = peak.to.test)

    return(res)

  }

}


################################################################################

gene_test_wrap = function(genes.multi.pk, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct, cut.pval, cut.MPRO) {

  res.tbl = t(sapply(1:length(genes.multi.pk), function(t) {

    if (t %% 1000 == 0) print(t)

    res = single_gene_test(genes.multi.pk[t], cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct)

    if (length(res) == 0) return(rep(NA, 11)) else {
      return(c(res[['pval.chisq']],
               res[['pval.DM']],
               res[['warm.base']],
               res[['warm.alt']],
               res[['delta.warm']],
               res[['low.pct.flag']],
               res[['MPRO']],
               res[['MPRO.id']],
               res[['prp.dif.base']],
               res[['prp.dif.alt']],
               paste(res[['pks']], collapse = ', ')))
    }

  }
  )
  ) # end transpose

  message("The number of rows in the result table: ", nrow(res.tbl))
  colnames(res.tbl) = c('pval.chisq', 'pval.DM', 'warm.base', 'warm.alt', 'delta.warm', 'low.pct.flag', 'MPRO', 'MPRO.id', 'prp.dif.base', 'prp.dif.alt', 'pks')
  rownames(res.tbl) = genes.multi.pk

  # assign proper data types
  res.tbl = as.data.frame(res.tbl, stringsAsFactors = FALSE)

  res.tbl <- convert.magic(res.tbl, c(rep('num', 5), 'fac', 'num', 'char', 'num', 'num', 'char'))
  res.tbl$gene = as.character(genes.multi.pk)

  # multiple testing correction
  res.sub.tbl = res.tbl[res.tbl$low.pct.flag == FALSE, c('gene', 'pval.chisq', 'pval.DM')]
  res.sub.tbl$pval.chisq.adj = stats::p.adjust(res.sub.tbl$pval.chisq, method = 'BH')
  res.sub.tbl$pval.DM.adj = stats::p.adjust(res.sub.tbl$pval.DM, method = 'BH')
  res.tbl = res.tbl %>% dplyr::left_join(res.sub.tbl[, c('gene', 'pval.chisq.adj', 'pval.DM.adj')], by = 'gene')
  MPRO = NULL
  res.tbl = res.tbl %>% dplyr::arrange(dplyr::desc(MPRO))

  # pval.inter
  res.tbl = res.tbl[!is.na(res.tbl$MPRO) & (!is.na(res.tbl$pval.DM) | !is.na(res.tbl$pval.chisq)), ]
  res.tbl$pval.inter = ifelse(!is.na(res.tbl$pval.DM) & !is.na(res.tbl$pval.chisq), pmax(res.tbl$pval.DM, res.tbl$pval.chisq), ifelse(is.na(res.tbl$pval.DM), res.tbl$pval.chisq, res.tbl$pval.DM))
  res.sub.tbl = res.tbl[res.tbl$low.pct.flag == FALSE, c('gene', 'pval.inter')]
  res.sub.tbl$pval.adj = stats::p.adjust(res.sub.tbl$pval.inter, method = 'BH')
  res.tbl = res.tbl %>% dplyr::left_join(res.sub.tbl[, c('gene', 'pval.adj')], by = 'gene')

  # filter by effect sizes
  res.tbl$sig = res.tbl$low.pct.flag == FALSE & res.tbl$pval.adj < cut.pval & abs(res.tbl$MPRO) > cut.MPRO
  message("The number of significant genes: ", sum(res.tbl$sig))
  message(paste0("long.vs.short: ", sum(res.tbl$MPRO > 0 & res.tbl$sig == T), ",", sum(res.tbl$MPRO < 0 & res.tbl$sig == T)))

  return(res.tbl)


}


################################################################################

single_utr_test <- function(utr.to.test, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct) {

  # find associated peaks
  peak.to.test = rownames(anno.tbl)[which(anno.tbl$gene.trans.ID == utr.to.test)]
  # order by mode.pos (small -> large)
  peak.to.test = peak.to.test[order(anno.tbl[peak.to.test, "mode.pos"])]
  strand = utils::head(anno.tbl[peak.to.test, 'strand'], 1)
  mode.pos = anno.tbl[peak.to.test, 'mode.pos']

  # generate df.test
  df.test = data.frame(cellid = colnames(cnt), cl = attr.tbl[, ctype.colname], stringsAsFactors = FALSE)
  for (i in 1:length(peak.to.test)) df.test[, paste0("raw|", peak.to.test[i])] = cnt[peak.to.test[i], ]
  df.test[, paste0("raw|", utr.to.test)] = rowSums(df.test[, 3:(3 + length(peak.to.test) - 1)])

  # generate aggregated count/percentage table for chi-square test
  cnt.tbl = pct.tbl = matrix(NA, ncol = 2, nrow = length(peak.to.test))
  rownames(cnt.tbl) = rownames(pct.tbl) = peak.to.test
  if ( is.null(alt_grp)) colnames(cnt.tbl) = colnames(pct.tbl) = c(base_grp, paste0('non-', base_grp))
  if (!is.null(alt_grp)) colnames(cnt.tbl) = colnames(pct.tbl) = c(base_grp, alt_grp)

  for (i in 1:length(peak.to.test)){

    df.test[, peak.to.test[i]] = df.test[, (2+i)] / df.test[, paste0("raw|", utr.to.test)]
    c1 = sum(df.test[which(df.test$cl == base_grp), 2+i])
    if ( is.null(alt_grp)) c2 = sum(df.test[which(df.test$cl != base_grp), 2+i])
    if (!is.null(alt_grp)) c2 = sum(df.test[which(df.test$cl ==  alt_grp), 2+i])
    cnt.tbl[i,] = c(c1, c2)
    p1 = mean(df.test[df.test$cl == base_grp, 2+i] > 0)
    if ( is.null(alt_grp)) p2 = mean(df.test[df.test$cl != base_grp, 2+i] > 0)
    if (!is.null(alt_grp)) p2 = mean(df.test[df.test$cl ==  alt_grp, 2+i] > 0)
    pct.tbl[i, ] = c(p1, p2)

  }

  if (all(colSums(cnt.tbl) != 0) && all(rowSums(cnt.tbl) != 0)) {

    # chisq test
    if (all(cnt.tbl == 0 )) pval.chisq = 1 else { chisq <- stats::chisq.test(cnt.tbl); pval.chisq = chisq$p.value }

    # multinomial-dirichlet LRT
    subdata.1 = as.matrix(df.test[df.test$cl == base_grp, 3:(3 + length(peak.to.test) - 1)])
    if ( is.null(alt_grp)) subdata.2 = as.matrix(df.test[df.test$cl != base_grp, 3:(3 + length(peak.to.test) - 1)])
    if (!is.null(alt_grp)) subdata.2 = as.matrix(df.test[df.test$cl ==  alt_grp, 3:(3 + length(peak.to.test) - 1)])

    alt.llk = c(f.logL(subdata.1), f.logL(subdata.2))
    names(alt.llk) <- colnames(cnt.tbl)

    if (all(!is.na(alt.llk))) {
      subdata.0 = as.matrix(df.test[, 3:(3 + length(peak.to.test) - 1)]) # under the null
      null.llk = f.logL(subdata.0)
      ll.alt <- sum(alt.llk, na.rm = T)
      ll.null = null.llk
      lr <- -2 * (ll.null - ll.alt)
      pval.DM <- stats::pchisq(lr, length(peak.to.test) * (2-1), lower.tail = FALSE)
    } else {
      pval.DM = NA
    }

    # rank-based WARM
    if (strand == '-') pos.rel = (-rank(mode.pos) + length(peak.to.test)) / (length(peak.to.test) - 1)
    if (strand == '+') pos.rel = ( rank(mode.pos) - 1) / (length(peak.to.test) - 1)
    warm.grp1  = sum(pos.rel * cnt.tbl[, base_grp] / sum(cnt.tbl[, base_grp]))
    if ( is.null(alt_grp)) warm.grp2  = sum(pos.rel * cnt.tbl[, paste0('non-', base_grp)] / sum(cnt.tbl[, paste0('non-', base_grp)]))
    if (!is.null(alt_grp)) warm.grp2  = sum(pos.rel * cnt.tbl[, alt_grp] / sum(cnt.tbl[, alt_grp]))
    delta.warm = warm.grp1 - warm.grp2

    # check simultaneous low pct in either cell-type tested
    low.pct.flag = apply(pct.tbl, 2, function(t) {all(t < cut.low.pct)})
    low.pct.flag = any(low.pct.flag)

    # MPRO
    comb = t(utils::combn(1:length(peak.to.test), 2))
    prp.tbl = t(t(cnt.tbl) / colSums(cnt.tbl))
    pd.vec = NULL

    for (i in 1:nrow(comb)) {
      this.prp = prp.tbl[comb[i, ], ]
      this.pd = this.prp[1, 1] + this.prp[2, 2] - this.prp[2, 1] - this.prp[1, 2]
      pd.vec = c(pd.vec, this.pd)
    }

    id.pd = which.max(abs(pd.vec))
    strand.multiplier = 2*(strand == '-') - 1
    max.abs.prpchg.id = ifelse(strand == '-', paste0(comb[id.pd, 1], ' vs ', comb[id.pd, 2]), paste0(comb[id.pd, 2], ' vs ', comb[id.pd, 1]))
    max.prpchg = pd.vec[id.pd] * strand.multiplier
    prp.dif.1 = ifelse(strand == '-', prp.tbl[comb[id.pd, 1], 1] - prp.tbl[comb[id.pd, 2], 1], prp.tbl[comb[id.pd, 2], 1] - prp.tbl[comb[id.pd, 1], 1])
    prp.dif.2 = ifelse(strand == '-', prp.tbl[comb[id.pd, 1], 2] - prp.tbl[comb[id.pd, 2], 2], prp.tbl[comb[id.pd, 2], 2] - prp.tbl[comb[id.pd, 1], 2])
    # prp.dif.1 - prp.dif.2 == max.prpchg

    res = list(cnt.tbl = cnt.tbl,
               pct.tbl = pct.tbl,
               pval.chisq = pval.chisq,
               pval.DM = pval.DM,
               warm.base = warm.grp1,
               warm.alt = warm.grp2,
               delta.warm = delta.warm,
               low.pct.flag = low.pct.flag,
               MPRO.id = max.abs.prpchg.id,
               MPRO = max.prpchg,
               prp.dif.base = prp.dif.1,
               prp.dif.alt = prp.dif.2,
               pks = peak.to.test)

    return(res)

  }

}

################################################################################

utr_test_wrap = function(utr.multi.pk, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct, cut.pval, cut.MPRO) {

  res.tbl = t(sapply(1:length(utr.multi.pk), function(t) {

    if (t %% 1000 == 0) print(t)

    res = single_utr_test(utr.multi.pk[t], cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct)

    if (length(res) == 0) return(rep(NA, 11)) else {
      return(c(res[['pval.chisq']],
               res[['pval.DM']],
               res[['warm.base']],
               res[['warm.alt']],
               res[['delta.warm']],
               res[['low.pct.flag']],
               res[['MPRO']],
               res[['MPRO.id']],
               res[['prp.dif.base']],
               res[['prp.dif.alt']],
               paste(res[['pks']], collapse = ', ')))
    }

  }
  )
  ) # end transpose

  message("The number of rows in the result table: ", nrow(res.tbl))
  colnames(res.tbl) = c('pval.chisq', 'pval.DM', 'warm.base', 'warm.alt', 'delta.warm', 'low.pct.flag', 'MPRO', 'MPRO.id', 'prp.dif.base', 'prp.dif.alt', 'pks')
  rownames(res.tbl) = utr.multi.pk

  # assign proper data types
  res.tbl = as.data.frame(res.tbl, stringsAsFactors = FALSE)

  res.tbl <- convert.magic(res.tbl, c(rep('num', 5), 'fac', 'num', 'char', 'num', 'num', 'char'))
  res.tbl$Transcript.ID = as.character(utr.multi.pk)
  res.tbl$gene = sapply(as.character(utr.multi.pk), function(t) strsplit(t, ':')[[1]][1])

  # multiple testing correction
  res.sub.tbl = res.tbl[res.tbl$low.pct.flag == FALSE, c('Transcript.ID', 'pval.chisq', 'pval.DM')]
  res.sub.tbl$pval.chisq.adj = stats::p.adjust(res.sub.tbl$pval.chisq, method = 'BH')
  res.sub.tbl$pval.DM.adj = stats::p.adjust(res.sub.tbl$pval.DM, method = 'BH')
  res.tbl = res.tbl %>% dplyr::left_join(res.sub.tbl[, c('Transcript.ID', 'pval.chisq.adj', 'pval.DM.adj')], by = 'Transcript.ID')
  MPRO = NULL
  res.tbl = res.tbl %>% dplyr::arrange(dplyr::desc(MPRO))

  # pval.inter
  res.tbl = res.tbl[!is.na(res.tbl$MPRO) & (!is.na(res.tbl$pval.DM) | !is.na(res.tbl$pval.chisq)), ]
  res.tbl$pval.inter = ifelse(!is.na(res.tbl$pval.DM) & !is.na(res.tbl$pval.chisq), pmax(res.tbl$pval.DM, res.tbl$pval.chisq), ifelse(is.na(res.tbl$pval.DM), res.tbl$pval.chisq, res.tbl$pval.DM))
  res.sub.tbl = res.tbl[res.tbl$low.pct.flag == FALSE, c('Transcript.ID', 'pval.inter')]
  res.sub.tbl$pval.adj = stats::p.adjust(res.sub.tbl$pval.inter, method = 'BH')
  res.tbl = res.tbl %>% dplyr::left_join(res.sub.tbl[, c('Transcript.ID', 'pval.adj')], by = 'Transcript.ID')

  # filter by effect sizes
  res.tbl$sig = res.tbl$low.pct.flag == FALSE & res.tbl$pval.adj < cut.pval & abs(res.tbl$MPRO) > cut.MPRO
  message("The number of significant utrs: ", sum(res.tbl$sig))
  message(paste0("long.vs.short: ", sum(res.tbl$MPRO > 0 & res.tbl$sig == T), ",", sum(res.tbl$MPRO < 0 & res.tbl$sig == T)))

  return(res.tbl)

}

################################################################################


single_btw_test <- function(gene.to.test, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct) {

  # find associated peaks
  peak.to.test = anno.tbl[anno.tbl$gene == gene.to.test, c('gene.trans.ID', 'strand', 'mode.pos')]
  peak.to.test$polyA_ID = rownames(peak.to.test)
  # order by mode.pos (small -> large)
  mode.pos = NULL
  peak.to.test = peak.to.test %>% dplyr::arrange(mode.pos)
  strand = utils::head(peak.to.test$strand, 1)
  trans = as.character(unique(peak.to.test$gene.trans.ID))
  n.trans = length(trans)

  # utr by cell count sub-matrix
  peak.to.test$gene.trans.ID.copy = peak.to.test$gene.trans.ID
  per.utr.list = split(peak.to.test[, -ncol(peak.to.test)], peak.to.test[, ncol(peak.to.test)])
  sub.cnt.list = lapply(per.utr.list, function(t) colSums(matrix(cnt[t$polyA_ID, ], nrow = nrow(t))))
  sub.cnt = do.call(what = 'rbind', args = sub.cnt.list)

  # generate df.test
  df.test = data.frame(cellid = colnames(cnt), cl = attr.tbl[, ctype.colname], stringsAsFactors = FALSE)
  for (i in 1:n.trans) df.test[, paste0("raw|", trans[i])] = sub.cnt[trans[i], ]
  df.test[, paste0("raw|", gene.to.test)] = rowSums(df.test[, 3:(3 + n.trans - 1)])

  # generate aggregated count/percentage table for chi-square test
  cnt.tbl = pct.tbl = matrix(NA, ncol = 2, nrow = n.trans)
  rownames(cnt.tbl) = rownames(pct.tbl) = trans
  if ( is.null(alt_grp)) colnames(cnt.tbl) = colnames(pct.tbl) = c(base_grp, paste0('non-', base_grp))
  if (!is.null(alt_grp)) colnames(cnt.tbl) = colnames(pct.tbl) = c(base_grp, alt_grp)

  for (i in 1:n.trans){

    df.test[, trans[i]] = df.test[, (2+i)] / df.test[, paste0("raw|", gene.to.test)]
    c1 = sum(df.test[which(df.test$cl == base_grp), 2+i])
    if ( is.null(alt_grp)) c2 = sum(df.test[which(df.test$cl != base_grp), 2+i])
    if (!is.null(alt_grp)) c2 = sum(df.test[which(df.test$cl ==  alt_grp), 2+i])
    cnt.tbl[i,] = c(c1, c2)
    p1 = mean(df.test[df.test$cl == base_grp, 2+i] > 0)
    if ( is.null(alt_grp)) p2 = mean(df.test[df.test$cl != base_grp, 2+i] > 0)
    if (!is.null(alt_grp)) p2 = mean(df.test[df.test$cl ==  alt_grp, 2+i] > 0)
    pct.tbl[i, ] = c(p1, p2)

  }

  if (all(colSums(cnt.tbl) != 0) && all(rowSums(cnt.tbl) != 0)) {

    # chisq test
    if (all(cnt.tbl == 0 )) pval.chisq = 1 else { chisq <- stats::chisq.test(cnt.tbl); pval.chisq = chisq$p.value }

    # multinomial-dirichlet LRT
    subdata.1 = as.matrix(df.test[df.test$cl == base_grp, 3:(3 + n.trans - 1)])
    if ( is.null(alt_grp)) subdata.2 = as.matrix(df.test[df.test$cl != base_grp, 3:(3 + n.trans - 1)])
    if (!is.null(alt_grp)) subdata.2 = as.matrix(df.test[df.test$cl ==  alt_grp, 3:(3 + n.trans - 1)])

    alt.llk = c(f.logL(subdata.1), f.logL(subdata.2))
    names(alt.llk) <- colnames(cnt.tbl)

    if (all(!is.na(alt.llk))) {
      subdata.0 = as.matrix(df.test[, 3:(3 + n.trans - 1)]) # under the null
      null.llk = f.logL(subdata.0)
      ll.alt <- sum(alt.llk, na.rm = T)
      ll.null = null.llk
      lr <- -2 * (ll.null - ll.alt)
      pval.DM <- stats::pchisq(lr, n.trans * (2-1), lower.tail = FALSE)
    } else {
      pval.DM = NA
    }

    # rank-based WARM
    if (strand == '-') pos.rel = (-rank(mode.pos) + n.trans) / (n.trans - 1)
    if (strand == '+') pos.rel = ( rank(mode.pos) - 1) / (n.trans - 1)
    warm.grp1  = sum(pos.rel * cnt.tbl[, base_grp] / sum(cnt.tbl[, base_grp]))
    if ( is.null(alt_grp)) warm.grp2  = sum(pos.rel * cnt.tbl[, paste0('non-', base_grp)] / sum(cnt.tbl[, paste0('non-', base_grp)]))
    if (!is.null(alt_grp)) warm.grp2  = sum(pos.rel * cnt.tbl[, alt_grp] / sum(cnt.tbl[, alt_grp]))
    delta.warm = warm.grp1 - warm.grp2

    # check simultaneous low pct in either cell-type tested
    low.pct.flag = apply(pct.tbl, 2, function(t) {all(t < cut.low.pct)})
    low.pct.flag = any(low.pct.flag)

    # MPRO
    comb = t(utils::combn(1:n.trans, 2))
    prp.tbl = t(t(cnt.tbl) / colSums(cnt.tbl))
    pd.vec = NULL

    for (i in 1:nrow(comb)) {
      this.prp = prp.tbl[comb[i, ], ]
      this.pd = this.prp[1, 1] + this.prp[2, 2] - this.prp[2, 1] - this.prp[1, 2]
      pd.vec = c(pd.vec, this.pd)
    }

    id.pd = which.max(abs(pd.vec))
    strand.multiplier = 2*(strand == '-') - 1
    max.abs.prpchg.id = ifelse(strand == '-', paste0(comb[id.pd, 1], ' vs ', comb[id.pd, 2]), paste0(comb[id.pd, 2], ' vs ', comb[id.pd, 1]))
    max.prpchg = pd.vec[id.pd] * strand.multiplier
    prp.dif.1 = ifelse(strand == '-', prp.tbl[comb[id.pd, 1], 1] - prp.tbl[comb[id.pd, 2], 1], prp.tbl[comb[id.pd, 2], 1] - prp.tbl[comb[id.pd, 1], 1])
    prp.dif.2 = ifelse(strand == '-', prp.tbl[comb[id.pd, 1], 2] - prp.tbl[comb[id.pd, 2], 2], prp.tbl[comb[id.pd, 2], 2] - prp.tbl[comb[id.pd, 1], 2])
    # prp.dif.1 - prp.dif.2 == max.prpchg

    res = list(cnt.tbl = cnt.tbl,
               pct.tbl = pct.tbl,
               pval.chisq = pval.chisq,
               pval.DM = pval.DM,
               warm.base = warm.grp1,
               warm.alt = warm.grp2,
               delta.warm = delta.warm,
               low.pct.flag = low.pct.flag,
               MPRO.id = max.abs.prpchg.id,
               MPRO = max.prpchg,
               prp.dif.base = prp.dif.1,
               prp.dif.alt = prp.dif.2,
               trans = trans)

    return(res)

  }

}

################################################################################

btw_test_wrap = function(genes.multi.utr, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct, cut.pval, cut.MPRO) {

  res.tbl = t(sapply(1:length(genes.multi.utr), function(t) {

    if (t %% 1000 == 0) print(t)

    res = single_btw_test(genes.multi.utr[t], cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct)

    if (length(res) == 0) return(rep(NA, 11)) else {
      return(c(res[['pval.chisq']],
               res[['pval.DM']],
               res[['warm.base']],
               res[['warm.alt']],
               res[['delta.warm']],
               res[['low.pct.flag']],
               res[['MPRO']],
               res[['MPRO.id']],
               res[['prp.dif.base']],
               res[['prp.dif.alt']],
               paste(res[['trans']], collapse = ', ')))
    }

  }
  )
  ) # end transpose

  message("The number of rows in the result table: ", nrow(res.tbl))
  colnames(res.tbl) = c('pval.chisq', 'pval.DM', 'warm.base', 'warm.alt', 'delta.warm', 'low.pct.flag', 'MPRO', 'MPRO.id', 'prp.dif.base', 'prp.dif.alt', 'trans')
  rownames(res.tbl) = genes.multi.utr

  # assign proper data types
  res.tbl = as.data.frame(res.tbl, stringsAsFactors = FALSE)

  res.tbl <- convert.magic(res.tbl, c(rep('num', 5), 'fac', 'num', 'char', 'num', 'num', 'char'))
  res.tbl$gene = as.character(genes.multi.utr)

  # multiple testing correction
  res.sub.tbl = res.tbl[res.tbl$low.pct.flag == FALSE, c('gene', 'pval.chisq', 'pval.DM')]
  res.sub.tbl$pval.chisq.adj = stats::p.adjust(res.sub.tbl$pval.chisq, method = 'BH')
  res.sub.tbl$pval.DM.adj = stats::p.adjust(res.sub.tbl$pval.DM, method = 'BH')
  res.tbl = res.tbl %>% dplyr::left_join(res.sub.tbl[, c('gene', 'pval.chisq.adj', 'pval.DM.adj')], by = 'gene')
  MPRO = NULL
  res.tbl = res.tbl %>% dplyr::arrange(dplyr::desc(MPRO))

  # pval.inter
  res.tbl = res.tbl[!is.na(res.tbl$MPRO) & (!is.na(res.tbl$pval.DM) | !is.na(res.tbl$pval.chisq)), ]
  res.tbl$pval.inter = ifelse(!is.na(res.tbl$pval.DM) & !is.na(res.tbl$pval.chisq), pmax(res.tbl$pval.DM, res.tbl$pval.chisq), ifelse(is.na(res.tbl$pval.DM), res.tbl$pval.chisq, res.tbl$pval.DM))
  res.sub.tbl = res.tbl[res.tbl$low.pct.flag == FALSE, c('gene', 'pval.inter')]
  res.sub.tbl$pval.adj = stats::p.adjust(res.sub.tbl$pval.inter, method = 'BH')
  res.tbl = res.tbl %>% dplyr::left_join(res.sub.tbl[, c('gene', 'pval.adj')], by = 'gene')

  # filter by effect sizes
  res.tbl$sig = res.tbl$low.pct.flag == FALSE & res.tbl$pval.adj < cut.pval & abs(res.tbl$MPRO) > cut.MPRO
  message("The number of significant genes: ", sum(res.tbl$sig))
  message(paste0("long.vs.short: ", sum(res.tbl$MPRO > 0 & res.tbl$sig == T), ",", sum(res.tbl$MPRO < 0 & res.tbl$sig == T)))

  return(res.tbl)

}


################################################################################

# function to calculate logL under the alternative in multinomial-dirichlet LRT test

f.logL = function(dat) {

  tryCatch({
    tmp.1 <- MGLM::MGLMfit(dat, dist = "DM")
    if(is.nan(sum(tmp.1@SE)) == FALSE && sum(tmp.1@SE) != Inf && sum(tmp.1@estimate > 1e-6) == length(tmp.1@estimate)) return(tmp.1@logL)
    else return(NA)
  }, error = function(e) return(NA)
  )

}

################################################################################

# Column data type conversion

convert.magic <- function(obj, types) {
  for (i in 1:ncol(obj)){
    FUN <- switch(types[i], char = as.character, num = as.numeric, fac = as.factor)
    obj[, i] <- FUN(obj[, i])
  }
  obj
}




