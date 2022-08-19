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
#'
#' @return A peak table saved as \[output.path\]/peaks.\[suffix\].txt.
#' @export
#'
#' @examples
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
                      suffix
) {

  # Peak calling
  peak.sites = peak_calling(genome.ref, bam, batch.start, batch.end, ncores, d, h, d.cut, hr, min.mode.prop, min.mode.cutoff)

  message("\n\nNumber of raw peaks identified: ", nrow(peak.sites))
  utils::write.csv(peak.sites, paste0(output.path, '/peaks_', suffix, '.txt'))
  message("\n\nPeak table is output? ", file.exists(paste0(output.path, '/peaks_', suffix, '.txt')))
  print(utils::head(peak.sites))

  # Peak annotation


  return (1)

}

