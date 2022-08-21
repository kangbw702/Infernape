
#' Count reads for each peak site
#'
#' @param bamfile BAM file.
#' @param whitelist.file Whitelist of cell barcodes.
#' @param peak.sites.file Peak site table after filtering.
#' @param ncores Number of cores.
#' @param start.cid The first peak cluster ID to analyze.
#' @param end.cid The last peak cluster ID to analyze.
#'
#' @return A peak site by cell UMI count matrix (sparse).
#' @export
#'
#' @examples
#' ls()
#'
#'
peak_counting <- function(bamfile,
                          whitelist.file ,
                          peak.sites.file,
                          ncores = ncores,
                          start.cid = NULL,
                          end.cid = NULL
                          ) {

  lock <- tempfile()
  whitelist.bc = utils::read.csv(whitelist.file, header = T, stringsAsFactors = FALSE, row.names = 1)
  whitelist.bc = whitelist.bc$x
  n.bcs = length(whitelist.bc)
  message("There are ", n.bcs, " whitelist barcodes.")

  # read in filtered peak annotation table
  peak.sites = utils::read.csv(peak.sites.file, header = TRUE, stringsAsFactors = FALSE)
  peak.sites$strand = ifelse(peak.sites$strand == "+", 1, -1)
  cids = unique(peak.sites$cluster_ID)
  message(paste0('Total - n.genes: ', length(unique(peak.sites$gene)), '; n.peak-clusters: ', length(cids), '; n.peaks: ', nrow(peak.sites)))
  # n.genes: 14243; n.peak-clusters: 21142; n.peaks: 35233

  # optional: select a subset of cids
  if (!is.null(start.cid) & !is.null(end.cid)) cids = cids[start.cid:min(end.cid, length(cids))]
  message(paste0('This run - n.peak-clusters: ', length(cids), '; n.peaks: ', nrow(peak.sites[peak.sites$cluster_ID %in% cids, ])))

  # Set up multiple workers
  system.name = Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == 'Windows') {
    new_cl = TRUE
    cluster = parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores = ncores)
  }

  `%dopar%` <- foreach::`%dopar%` # %dopar% doesn't need to be attached
  cluster.id = NULL
  mat.to.write <- foreach::foreach(cluster.id = cids, .combine = 'rbind', .packages = c("magrittr")) %dopar% {

    peak.sites.sub = peak.sites[peak.sites$cluster_ID == cluster.id, ]
    cluster_ID = gene = seq = from = to = mu = sigma = polyA_ID = NULL
    peak.sites.sub = peak.sites.sub %>% dplyr::select(gene, seq, from, to, strand, mu, sigma, pi, polyA_ID, cluster_ID)
    strand = peak.sites.sub$strand[1]; seq.name = peak.sites.sub$seq[1]
    n.pks = nrow(peak.sites.sub)
    message(paste0(which(unique(peak.sites$cluster_ID) == cluster.id), 'th cluster: ', peak.sites.sub$cluster_ID[1], '; seq:', seq.name, '; strand: ', strand, '; n.pks:', n.pks))

    isMinusStrand = if(strand == 1) FALSE else TRUE
    which = GenomicRanges::GRanges(seqnames = seq.name, ranges = IRanges::IRanges(min(peak.sites.sub$from), max(peak.sites.sub$to)))
    param = Rsamtools::ScanBamParam(tag = c('CB', 'UB'), which = which, flag = Rsamtools::scanBamFlag(isMinusStrand = isMinusStrand))
    aln = GenomicAlignments::readGAlignments(bamfile, param = param)

    if (length(aln) > 0) {

      GenomicRanges::mcols(aln)$CB_UB = paste0(unlist(GenomicRanges::mcols(aln)$CB), "_", unlist(GenomicRanges::mcols(aln)$UB))

      # remove null CB/UB
      na.cb = which(unlist(GenomicRanges::mcols(aln)$CB == 'NA' | is.na(GenomicRanges::mcols(aln)$CB)))
      na.ub = which(unlist(GenomicRanges::mcols(aln)$UB == 'NA' | is.na(GenomicRanges::mcols(aln)$UB)))
      na.cb.ub = dplyr::union(na.cb, na.ub)
      if (length(na.cb.ub)>0) aln = aln[-na.cb.ub]

      if (length(aln) > 0) {

        # check CB nchar consistency
        nchar.CBtag = nchar(unlist(GenomicRanges::mcols(aln[1])['CB']))
        nchar.white = nchar(whitelist.bc[1])
        if (nchar.CBtag != nchar.white) warning('CB lengths differs in BAM and whitelist!'); stop

        # only keep CB in whitelist
        whitelist.pos = which(unlist(GenomicRanges::mcols(aln)['CB']) %in% whitelist.bc)

        if (length(whitelist.pos) > 0) {

          aln = aln[whitelist.pos]

          # remove intronic reads that do not overlap with [min(from), max(to)]
          range.M = GenomicAlignments::cigarRangesAlongReferenceSpace(GenomicAlignments::cigar(aln),  pos = GenomicAlignments::start(aln), ops = "M", with.ops = TRUE, reduce.ranges = F)
          gene.subject = IRanges::IRanges(start = min(peak.sites.sub$from), end = max(peak.sites.sub$to))
          overlap.flag = sapply(range.M, function(t) any(IRanges::overlapsAny(t, gene.subject)))

          if (sum(overlap.flag) > 0) {

            aln = aln[overlap.flag]

            # de-duplication
            aln.df = as.data.frame(aln)
            aln.df$id = 1:nrow(aln.df)
            aln.df.split = split(aln.df, f = aln.df$CB_UB)
            ids.overlap = sapply(aln.df.split, function(t) {
            t = t[order(t$start), ]
            m = ifelse(nrow(t) %% 2 == 1, (nrow(t) + 1) / 2, nrow(t) / 2)
            return(t[m, 'id'])
            })

            aln = aln[ids.overlap]

            if (n.pks == 1) {

              aln.per.cell = GenomicRanges::split(aln, unlist(GenomicRanges::mcols(aln)['CB']))
              polyA.GR = GenomicRanges::GRanges(seqnames = peak.sites.sub$seq, IRanges::IRanges(start = as.integer(peak.sites.sub$from), end = as.integer(peak.sites.sub$to)))
              barcodes.gene = names(aln.per.cell)
              res = sapply(barcodes.gene, function(x) GenomicRanges::countOverlaps(polyA.GR, aln.per.cell[[x]])) # integer of length = n.cells

              # Reorder the columns of the res matrix to match the whitelist barcodes
              res.mat = matrix(0L, nrow = n.pks, ncol = n.bcs)
              res.mat[, match(barcodes.gene, whitelist.bc)] = res
              rownames(res.mat) = as.character(peak.sites.sub$polyA_ID)
              res.mat = Matrix::Matrix(res.mat, sparse = TRUE)

            } # end n.pks == 1

            if (n.pks > 1) {

              # re-weight
              peak.sites.sub$pi = peak.sites.sub$pi / sum(peak.sites.sub$pi)
              mid = (GenomicRanges::start(aln) + GenomicRanges::end(aln))/2

              assign.id = sapply(mid, function(t) {
              if (mid < min(peak.sites.sub$from)) return (1)
              if (mid > max(peak.sites.sub$to)) return (nrow(peak.sites.sub))
              if (mid <= max(peak.sites.sub$to) & mid >= min(peak.sites.sub$from)) {
              prob = stats::dnorm(t, peak.sites.sub$mu, peak.sites.sub$sigma) * peak.sites.sub$pi
              return (which.max(prob))
              }
              })

              GenomicRanges::mcols(aln)$assign.id = assign.id
              aln.per.cell = GenomicRanges::split(aln, unlist(GenomicRanges::mcols(aln)['CB']))
              barcodes.gene = names(aln.per.cell)

              # res: matrix of dim = n.pks x n.cells
              res = sapply(barcodes.gene, function(x) {
              ids = unlist(GenomicRanges::mcols(aln.per.cell[[x]])['assign.id'])
              return(sapply(1:nrow(peak.sites.sub), function(t) sum(ids == t)))
              })

              # Reorder the columns of the res matrix to match the whitelist barcodes
              res.mat = matrix(0L, nrow = n.pks, ncol = n.bcs)
              res.mat[, match(barcodes.gene, whitelist.bc)] = res
              rownames(res.mat) = as.character(peak.sites.sub$polyA_ID)
              res.mat = Matrix::Matrix(res.mat, sparse = TRUE)

            } # end n.pks > 1

          } else { # all overlap.flag == FALSE
            res.mat = matrix(0L, nrow = n.pks, ncol = n.bcs)
            rownames(res.mat) = as.character(peak.sites.sub$polyA_ID)
            res.mat = Matrix::Matrix(res.mat, sparse = TRUE)
          }

        } else { # length(whitelist.pos) == 0
          res.mat = matrix(0L, nrow = n.pks, ncol = n.bcs)
          rownames(res.mat) = as.character(peak.sites.sub$polyA_ID)
          res.mat = Matrix::Matrix(res.mat, sparse = TRUE)
        }

      } else { # length(aln) == 0 after removing NA
        res.mat = matrix(0L, nrow = n.pks, ncol = n.bcs)
        rownames(res.mat) = as.character(peak.sites.sub$polyA_ID)
        res.mat = Matrix::Matrix(res.mat, sparse = TRUE)
      }

    } else { # length(aln) == 0
      res.mat = matrix(0L, nrow = n.pks, ncol = n.bcs)
      rownames(res.mat) = as.character(peak.sites.sub$polyA_ID)
      res.mat = Matrix::Matrix(res.mat, sparse = TRUE)
    }

    # Return sparse matrix for peaks of each peak-cluster
    return(res.mat)

  } # end dopar

  # Shut down cluster if on Windows
  if (new_cl) parallel::stopCluster(cluster)

  #if (!dir.exists(output.dir)) dir.create(output.dir)

  colnames(mat.to.write) = whitelist.bc
  return(mat.to.write)

}
