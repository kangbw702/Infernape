#' Detect peaks from sorted bam files
#'
#' @param genome.ref Genome reference gtf file.
#' @param bam Sorted bam file.
#' @param batch.start ID of the first gene in this batch.
#' @param batch.end ID of the last gene in this batch.
#' @param ncores The number of cores.
#' @param d A parameter for smoothing.
#' @param h A parameter for smoothing.
#' @param d.cut A parameter for merging close peaks.
#' @param hr The estimated half range of peak mode - PAS interval. Used for annotating peak modes.
#' @param min.mode.prop Minimal relative threshold (wrt highest peak in 3'UTR region) to filter out small peaks.
#' @param min.mode.cutoff Minimal absolute coverage threshold to filter out small peaks.
#' @param output.path Output directory.
#' @param suffix Suffix as unique label for output table.
#' @param pas.reference.file Known PAS table.
#' @param pas.reference.file Known PAS table
#' @param genome Full genome sequence. Required to be a S4 object.
#' @param pas.search.cut.1 Distance (nt) from the peak mode to upstream boundary of the searching interval for PAS
#' @param pas.search.cut.2 Distance (nt) from the peak mode to downstream boundary of the searching interval for PAS
#' @param polystretch_length Length of consecutive A sequence
#' @param max_mismatch Maximal tolerance of mismatch
#' @param motif.search.cut Window width for searching specified motifs
#' @param invert_strand default FALSE
#' @param q Vector of length 2 which defines the PAS searching window.
#' @param whitelist.file Whitelist of cell barcodes.
#' @param start.cid The first peak cluster ID to analyze.
#' @param end.cid The last peak cluster ID to analyze.
#'
#' @return A peak table saved as \[output.path\]/peaks_\[suffix\].txt.
#' A peak annotation table saved as \[output.path\]/anno_\[suffix\].txt.
#' A peak annotation table after fltering by PAS and motifs saved as \[output.path\]/anno_filtered_\[suffix\].txt.
#' A peak by cell UMI count sparse matrix.
#'
#' @export
#'
#' @examples
#' ls()
#'
#'
Infernape <- function(genome.ref,
                      bam,
                      batch.start = 1,
                      batch.end = 10,
                      ncores = 1,
                      d = 31,
                      h = 5,
                      d.cut = 50,
                      hr = 160,
                      min.mode.prop = 0.05,
                      min.mode.cutoff = 10,
                      output.path,
                      suffix,
                      pas.reference.file,
                      genome,
                      pas.search.cut.1 = 0,
                      pas.search.cut.2 = 300,
                      polystretch_length = 13,
                      max_mismatch = 1,
                      motif.search.cut = 300,
                      invert_strand = FALSE,
                      q = c(110, 200),
                      whitelist.file,
                      start.cid = NULL,
                      end.cid = NULL
) {

  # Peak calling
  peak.sites = peak_calling(genome.ref, bam, batch.start, batch.end, ncores, d, h, d.cut, hr, min.mode.prop, min.mode.cutoff)

  message("\n\nNumber of raw peaks identified: ", nrow(peak.sites))
  utils::write.csv(peak.sites, paste0(output.path, '/peaks_', suffix, '.txt'))
  message("\n\nPeak table is output? ", file.exists(paste0(output.path, '/peaks_', suffix, '.txt')))
  print(utils::head(peak.sites))

  # Peak annotation
  peak.sites.file = paste0(output.path, '/peaks_', suffix, '.txt')
  anno = peak_annotation(peak.sites.file, pas.reference.file, genome, pas.search.cut.1, pas.search.cut.2, polystretch_length, max_mismatch, motif.search.cut, invert_strand)

  message("\n\nPeak annotation completes.")
  utils::write.csv(anno, paste0(output.path, '/anno_', suffix, '.txt'))
  print(utils::head(anno))

  # Peak filtering
  anno.raw.file = paste0(output.path, '/anno_', suffix, '.txt')
  anno.filter = peak_filtering(anno.raw.file, q)
  utils::write.csv(anno.filter, paste0(output.path, '/anno_filtered_', suffix, '.csv'))
  print(utils::head(anno.filter))

  # peak counting
  peak.sites.file = paste0(output.path, '/anno_filtered_', suffix, '.csv')
  mat = peak_counting(bamfile = bam, whitelist.file, peak.sites.file, ncores = ncores, start.cid = NULL, end.cid = NULL)

  cnt.out.path = paste0(output.path, '/cnt_mat')
  if (!dir.exists(cnt.out.path)) dir.create(cnt.out.path)
  Matrix::writeMM(mat, file = paste0(cnt.out.path, "/matrix.mtx"))
  utils::write.table(colnames(mat), file = paste0(cnt.out.path, "/barcodes.tsv"),  quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(rownames(mat), file = paste0(cnt.out.path, "/sitenames.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)



  return ('DONE!')

}

