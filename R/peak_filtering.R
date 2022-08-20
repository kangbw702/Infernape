
#' Filter peaks
#'
#' @param anno.raw.file Peak annotation table.
#' @param q Vector of length 2 which defines the PAS searching window.
#'
#' @return A peak annotation table after filtering.
#' @export
#'
#' @examples
#' ls()
#'
peak_filtering <- function(anno.raw.file, q) {

  # PAS
  anno = utils::read.csv(anno.raw.file, row.names = 1, stringsAsFactors = FALSE)

  anno$PAS = sapply(anno$PAS_full, function(t) {
    PAS.vec = as.numeric(strsplit(t, ',')[[1]]); return (any(PAS.vec >= q[1] & PAS.vec <= q[2]))
  })

  message(paste0("n.PAS: ", sum(anno$PAS)))

  # polyA motif
  motif.col = names(anno)[grepl('^pA_motif_', names(anno))]
  motif.short = gsub('pA_motif_', '', motif.col)

  for (i in 1:length(motif.col)) {
    anno[, motif.short[i]] = sapply(as.character(anno[, motif.col[i]]), function(t) {
      motif.vec = as.numeric(strsplit(t, ',')[[1]]); return (any(motif.vec >= (q[1] - 20) & motif.vec <= (q[2] - 20)))
    })
    message(paste0("n.", motif.col[i], ": ", sum(anno[, motif.short[i]], na.rm = T)))
  }

  anno$motif = apply(anno[, motif.short], 1, any)
  message(paste0("n.motif: ", sum(anno$motif, na.rm = T)))

  # internal priming
  anno$As = sapply(as.character(anno$pA_stretch), function(t) {
    As.vec = as.numeric(strsplit(t, ',')[[1]]); return (any(As.vec >= q[1] & As.vec <= q[2]))
  })

  message(paste0("n.nonA_stretch: ", nrow(anno) - sum(anno$As, na.rm = T)))

  # novel PAS
  anno$novel = anno$PAS == F & anno$motif == T & (anno$As==F |is.na(anno$As))
  message(paste0("n.novel: ",sum(anno$novel, na.rm = T)))

  # final
  anno$select = anno$PAS == T | anno$motif == T & (anno$As==F |is.na(anno$As))
  message(paste0("n.keep: ", sum(anno$select, na.rm = T)))

  # n.PAS per gene
  tmp  = anno[anno$select == T & !is.na(anno$select), ]
  gene = NULL
  tmp = tmp %>% dplyr::group_by(gene) %>% dplyr::summarise(n.pas = dplyr::n())
  message(paste0("n.genes: ", length(unique(tmp$gene))))

  # pdf('n.pas.per.gene.pdf'); plot(table(tmp$n.pas), col = 'red3'); dev.off()
  anno.filter = anno[anno$select == TRUE & !is.na(anno$select), ]
  return(anno.filter)

}


################################################################################



#' Adaptively determine filtering parameters using single PAS genes as negative controls.
#'
#' @param anno.raw.file Peak annotation table.
#' @param pas.reference.file Known PAS table.
#' @param q.low Probability corresponding to quantile low bound
#' @param q.high Probability corresponding to quantile upper bound
#' @param plot.output.dir Directory to output plot.
#'
#' @return Quantile lower/upper bounds and a plot showing distribution of peak mode-PAS distance for negative control genes.
#' @export
#'
#' @examples
#' ls()
#'
#'
rule_PAS_mode_dist <- function(anno.raw.file, pas.reference.file, q.low, q.high, plot.output.dir) {

  # load raw peak annotation table
  pk = utils::read.table(anno.raw.file, header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
  pk$gene.id = paste0(pk$gene, ':', pk$seq, ':', pk$strand)

  # load reference table with single-PAS genes
  ref = utils::read.csv(pas.reference.file, stringsAsFactors = FALSE)
  ref = ref[ref$uniGene, c('mgi_symbol', 'Position', 'Chrom', 'Strand')]
  ref$gene.id = paste0(ref$mgi_symbol, ':', ref$Chrom, ':', ref$Strand)

  # find true negative controls
  # only keep single-peak gene in pk anno table
  ref$n.gene.nc = sapply(ref$gene.id, function(t) nrow(pk[pk$gene.id == t,]))
  ref = ref[ref$n.gene.nc == 1, ]

  # only keep genes with n.PAS=1 in pk
  tmp = pk[pk$gene.id %in% ref$gene.id, ]
  tmp = tmp[which(tmp$n_PAS == 1), ]
  ref = ref %>% dplyr::left_join(tmp[, c('gene.id', 'PAS_full', 'sigma')], by = 'gene.id')
  ref = ref[!is.na(ref$PAS_full), ]
  ref$d = as.numeric(ref$PAS_full)

  message("Number of negative control genes: ", nrow(ref))

  # plot
  d = ..count.. = NULL
  ggplot2::ggplot(ref, ggplot2::aes(x = d)) +
  ggplot2::geom_density(ggplot2::aes(y = ..count..), col = 'white', fill = 'orange', alpha = 0.5) +
  ggplot2::geom_vline(xintercept = stats::quantile(ref$d, c(q.low, q.high))) +
  ggplot2::theme_bw()
  ggplot2::ggsave(paste0(plot.output.dir, '/PAS.mode.distance.pdf'))

  return(list(n.nc = nrow(ref), q = round(stats::quantile(ref$d, c(q.low, q.high)), 0)))

}
