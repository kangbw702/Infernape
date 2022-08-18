#' Detect peaks from sorted bam files
#'
#' @param genome.ref genome reference gtf file.
#'
#' @return message showing the total number of genes in genome reference.
#' @export
#'
#' @examples
#'
#'
peak_calling = function(genome.ref) {

  # read in gtf files and extract information: gene symbol, start and end positions, chromosome, strand
  tictoc::tic("Import genome reference ...")
  ref.df = utils::read.csv(genome.ref)
  feature = chromosome = start = end = strand = gene_symbol = NULL
  gene.ref.df = ref.df %>% dplyr::filter(feature == 'gene') %>% dplyr::select(chromosome, start, end, strand, gene_symbol)
  n.genes = nrow(gene.ref.df)
  message(paste0("Total number of genes: ", n.genes))
  tictoc::toc()

}
