
# test for APA differential events between two cell groups (all legitimate genes/3'UTRs)

#' Differential APA testing
#'
#' @param counts.dir Directory for count matrix and corresponding peak site names and cell barcodes.
#' @param attr.file Directory for the cell attributes table.
#' @param anno.file Directory for the peak annotation table.
#' @param utr3.file Directory for the 3'UTR annotation table (mapping each 3'UTR region to a unique transcript ID).
#' @param ctype.colname Column name that defines cell groups in the cell attributes table.
#' @param base_grp Base cell group name.
#' @param alt_grp Alternative cell group name. If 'NULL', take all cell groups except base as alternative.
#' @param cut.low.pct Threshold for low expressing percentage. Default to be 0.05.
#' @param cut.pval Threshold for p-value. Default to be 0.05.
#' @param cut.MPRO Threshold for effect size (MPRO) filtering. Default to be 0.2.
#' @param test.type Type of testing. Value can be 'gene', 'utr', 'btw'.
#' @param out.dir Directory to save test results.
#'
#' @return Differential APA testing result table.
#' @export
#'
#' @examples
#' ls()
#'
#'
#'
#'
apa_test <- function(counts.dir,
                     attr.file,
                     anno.file,
                     utr3.file,
                     ctype.colname,
                     base_grp,
                     alt_grp = NULL,
                     cut.low.pct = 0.05,
                     cut.pval = 0.05,
                     cut.MPRO = 0.2,
                     test.type,
                     out.dir
                     )
  {

  # prepare the subset of counts, attr, anno tables that cope with the test

  cnt = Matrix::readMM(paste0(counts.dir, '/matrix.mtx'))
  colnames(cnt) = unlist(utils::read.table(paste0(counts.dir, '/barcodes.tsv'), sep = '\t', header = F))
  rownames(cnt) = unlist(utils::read.table(paste0(counts.dir, '/sitenames.tsv'), sep = '\t', header = F))

  # manipulate anno.tbl
  anno.tbl = utils::read.csv(anno.file, stringsAsFactors = FALSE)

  # assign each peak a Transcript.ID
  anno.tbl$Transcript.ID = sapply(1:nrow(anno.tbl), assign.trans.ID, utr3.file = utr3.file, anno.tbl = anno.tbl)

  # add summary of number of peaks/transcripts
  gene = gene.trans.ID = NULL
  anno.tbl$gene.trans.ID = paste0(anno.tbl$gene, ":", anno.tbl$Transcript.ID)
  anno.tbl.agg1 = anno.tbl %>% dplyr::group_by(gene) %>% dplyr::summarise(nn.pk = dplyr::n())
  anno.tbl.agg2 = anno.tbl %>% dplyr::group_by(gene.trans.ID) %>% dplyr::summarise(nn.pk.per.trans = dplyr::n(), gene = utils::head(gene, 1))
  anno.tbl.agg3 = anno.tbl.agg2 %>% dplyr::group_by(gene) %>% dplyr::summarise(nn.trans.per.gene = dplyr::n())
  anno.tbl = anno.tbl %>% dplyr::left_join(anno.tbl.agg1, by = 'gene') %>% dplyr::left_join(anno.tbl.agg2[, 1:2], by = 'gene.trans.ID')  %>% dplyr::left_join(anno.tbl.agg3, by = 'gene')

  rownames(anno.tbl) = anno.tbl$polyA_ID

  message('Peak site names in anno table are consistent with count matrix? ', all(rownames(cnt) == rownames(anno.tbl)))

  # manipulate attr.tbl
  attr.tbl = utils::read.csv(attr.file, stringsAsFactors = FALSE, row.names = 1)
  attr.tbl = attr.tbl[!duplicated(attr.tbl$cellbc), ]
  rownames(attr.tbl) = attr.tbl$cellbc
  shared.cb = intersect(colnames(cnt), rownames(attr.tbl))
  attr.tbl = attr.tbl[shared.cb, ]
  cnt = cnt[, shared.cb]

  message('Cell barcodes in attr table are consistent with count matrix? ', all(colnames(cnt) == rownames(attr.tbl)))

  out = TRUE

  if (!test.type %in% c('gene', 'utr', 'btw')) { warning("Wrong test type is given!"); out =F; stop }

  if (test.type == 'gene') {

    genes.multi.pk = unique(anno.tbl[anno.tbl$nn.pk > 1, 'gene'])
    message("The number of multi-peak genes to test: ", length(genes.multi.pk))
    if (length(genes.multi.pk) >= 1) {
      res = gene_test_wrap(genes.multi.pk, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct, cut.pval, cut.MPRO)
    } else {
      message('No multi-peak genes found.')
      out = F
    }

  }

  if (test.type == 'utr') {

    utr.multi.pk = unique(anno.tbl[anno.tbl$nn.pk.per.trans > 1, 'gene.trans.ID'])
    message("The number of multi-peak utrs to test: ", length(utr.multi.pk))
    if (length(utr.multi.pk) >= 1) {
      res = utr_test_wrap(utr.multi.pk, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct, cut.pval, cut.MPRO)
    } else {
      message('No multi-peak utrs found.')
      out = F
    }

  }

  if (test.type == 'btw') {

    genes.multi.utr = unique(anno.tbl[anno.tbl$nn.trans.per.gene > 1, 'gene'])
    message("The number of multi-utr genes to test: ", length(genes.multi.utr))
    if (length(genes.multi.utr) >= 1) {
      res = btw_test_wrap(genes.multi.utr, cnt, anno.tbl, attr.tbl, ctype.colname, base_grp, alt_grp, cut.low.pct, cut.pval, cut.MPRO)
    } else {
      message('No multi-utr genes found.')
      out = F
    }

  }

  if (out) {
    if ( is.null(alt_grp)) utils::write.csv(res, paste0(out.dir, '/', test.type, '.', base_grp, '.vs.others.csv'))
    if (!is.null(alt_grp)) utils::write.csv(res, paste0(out.dir, '/', test.type, '.', base_grp, '.vs.', alt_grp, '.csv'))
    message('Write out done!')
  }

}



assign.trans.ID <- function(pid, utr3.file, anno.tbl) {

  utr3.anno = utils::read.csv(utr3.file, stringsAsFactors = FALSE, row.names = 1)
  ids = which(utr3.anno$gene_symbol == anno.tbl[pid, 'gene'] & utr3.anno$start <= anno.tbl[pid, 'mode.pos'] & utr3.anno$end >= anno.tbl[pid, 'mode.pos'])

  if (length(ids) >  0) return(as.character(utr3.anno$Transcript.ID[ids]))

  if (length(ids) == 0) { # assign ID of the closest downstream 3'UTR

    tmp = utr3.anno[utr3.anno$gene_symbol == anno.tbl[pid, 'gene'], ]
    tmp = tmp[order(tmp$start), ]

    if (anno.tbl[id, 'strand'] == '-') {
      id = max(which(anno.tbl[pid, 'mode.pos'] - tmp$end > 0))
    } else {
      id = min(which(tmp$start - anno.tbl[pid, 'mode.pos'] > 0))
    }

    return(as.character(tmp$Transcript.ID[pid]))
  } # na ids

}


