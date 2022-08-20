
#' Count reads for each peak site
#'
#' @param bamfile BAM file.
#' @param whitelist.file Whitelist of cell barcodes.
#' @param peak.sites.file Peak site table after filtering.
#' @param ncores Number of cores.
#' @param CBtag Cell barcode tag in BAM file.
#' @param UMItag UMI barcode tag in BAM file.
#' @param start.cid The first peak cluster ID to analyze.
#' @param end.cid The last peak cluster ID to analyze.
#'
#' @return A peak site by cell UMI count matrix.
#' @export
#'
#' @examples
#' ls()
#'
#'
peak_counting <- function(bamfile=bamfile[i],
                          whitelist.file = whitelist.file[i],
                          peak.sites.file = peak.sites.file,
                          ncores = ncores,
                          CBtag = 'CB',
                          UMItag = 'UB',
                          start.cid = 1,
                          end.cid = 26316
                          ) {



}
