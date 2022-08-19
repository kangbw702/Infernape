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
#' @return A peak table.
#' @export
#'
#' @examples
#'
#'
peak_calling = function(genome.ref,
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

  # read in gtf files and extract information: gene symbol, start and end positions, chromosome, strand
  tictoc::tic("Import genome reference ...")
  ref.df = utils::read.csv(genome.ref)
  feature = chromosome = start = end = strand = gene_symbol = NULL
  gene.ref.df = ref.df %>% dplyr::filter(feature == 'gene') %>% dplyr::select(chromosome, start, end, strand, gene_symbol)
  n.genes = nrow(gene.ref.df)
  message(paste0("Total number of genes: ", n.genes))
  batch.end = min(n.genes, batch.end)
  message(paste("Genes to process in this batch: from", batch.start, "to", batch.end))
  tictoc::toc()

  # build gr (used later in annotation)
  gr_ref = GenomicRanges::GRanges(seqnames = ref.df$chromosome, ranges = IRanges::IRanges(ref.df$start, ref.df$end))
  gr_ref$feature = ref.df$feature
  gr_ref$gene_symbol = ref.df$gene_symbol

  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {

    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)

  } else {

    doParallel::registerDoParallel(cores = ncores)

  }

  tictoc::tic('Detecting peaks ... ')
  res = ii = NULL
  `%dopar%` <- foreach::`%dopar%` # %dopar% doesn't need to be attached
  res = foreach::foreach(ii = batch.start:batch.end, .packages = c("GenomicRanges"), .combine = "rbind") %dopar% {

    # which gene
    gene.name <- gene.ref.df[ii, "gene_symbol"]
    seq.name <- gene.ref.df[ii, "chromosome"]
    gene.start <- gene.ref.df[ii, "start"]
    gene.end <- gene.ref.df[ii, "end"]
    strand <- gene.ref.df[ii, "strand"]

    message(paste0("\n\n", ii, 'th gene ', gene.name, ' ', seq.name,' ', gene.start,' ', gene.end,' ', strand))

    # coverage table
    isMinusStrand <- if(strand == 1) FALSE else TRUE
    which <- GenomicRanges::GRanges(seqnames = seq.name, ranges = IRanges::IRanges(gene.start, gene.end))
    param <- Rsamtools::ScanBamParam(tag = c('CB', 'UB'), which = which, flag = Rsamtools::scanBamFlag(isMinusStrand = isMinusStrand))
    aln <- GenomicAlignments::readGAlignments(bam, param = param)

    if (length(aln) > 0) {

      GenomicRanges::mcols(aln)$CB_UB = paste0(unlist(GenomicRanges::mcols(aln)$CB), "_", unlist(GenomicRanges::mcols(aln)$UB))

      # remove null CB/UB
      na.cb = which(unlist(GenomicRanges::mcols(aln)$CB == 'NA' | is.na(GenomicRanges::mcols(aln)$CB)))
      na.ub = which(unlist(GenomicRanges::mcols(aln)$UB == 'NA' | is.na(GenomicRanges::mcols(aln)$UB)))
      na.cb.ub = dplyr::union(na.cb, na.ub)
      if (length(na.cb.ub)>0) aln = aln[-na.cb.ub]

      if (length(aln) > 0) {

        # de-duplicate
        aln.df = as.data.frame(aln)
        aln.df$id = 1:nrow(aln.df)
        aln.df.split = split(aln.df, f = aln.df$CB_UB)

        ids = sapply (
          aln.df.split, function(t) {
            t = t[order(t$start), ]
            m = ifelse(nrow(t) %% 2 == 1, (nrow(t) + 1) / 2, nrow(t) / 2)
            return(t[m, 'id'])
          }
        )

        aln_dd = aln[ids]
        aln_cov = GenomicRanges::coverage(aln_dd)[seq.name][[1]]

        if (length(aln_cov@values) > 1) {

          gene.end = min(gene.end, max(end(aln_cov)))
          data <- data.frame(pos = seq(gene.start, gene.end), coverage = S4Vectors::runValue(aln_cov)[S4Vectors::findRun(gene.start:gene.end, aln_cov)])
          message(paste0("Number of rows in raw coverage: ", nrow(data)))

          # smoothing and detect modes
          dat.mode = find.modes(dat = data, d = d, h = h, min.mode.cutoff = min.mode.cutoff, min.mode.prop = min.mode.prop, ref = ref.df, gene.name = gene.name)

          mode.select = dat.mode[which(dat.mode$select_mode == TRUE), c('pos','coverage')]
          names(mode.select) = c('mode.pos', 'mode.val')

          if (nrow(mode.select) > 0) {

            mode.select = mode.select %>% dplyr::mutate(gene = gene.name, seq = seq.name, strand = strand, gene.start = gene.start, gene.end = gene.end)
            mode.select$strand = plyr::mapvalues(x = mode.select$strand, from = c("1", "-1"), to = c("+", "-"))

            # merge modes that are too close to each other
            mode.select$distance = mode.select$mode.pos - dplyr::lag(mode.select$mode.pos, 1)
            mode.select <- merge_close_pk(mode.select, d.cut = d.cut)

            # annotate: UTR3_mode and CDS_mode

            # single position gr list for all modes
            gr_mode = GenomicRanges::GRanges(paste0(mode.select$seq,":", mode.select$mode.pos, ":", mode.select$strand))

            mode.select$UTR3_mode = ""
            my_UTR3 <- gr_ref[gr_ref$feature == "three_prime_UTR" & gr_ref$gene_symbol == gene.name]
            all_UTR3_mode_hits <- GenomicAlignments::findOverlaps(gr_mode, my_UTR3, type = 'any')
            idx_to_annotate_UTR3_mode <- S4Vectors::queryHits(all_UTR3_mode_hits)
            mode.select$UTR3_mode[unique(idx_to_annotate_UTR3_mode)] <- "YES"

            # Annotating CDS by mode
            mode.select$CDS_mode <- ""
            #my_CDS <- gr_ref[gr_ref$feature == "CDS" ]
            my_CDS <- gr_ref[gr_ref$feature == "CDS" & gr_ref$gene_symbol == gene.name]
            all_CDS_mode_hits <-  GenomicAlignments::findOverlaps(gr_mode, my_CDS, type = 'any')
            idx_to_annotate_CDS_mode <- S4Vectors::queryHits(all_CDS_mode_hits)
            mode.select$CDS_mode[unique(idx_to_annotate_CDS_mode)] <- "YES"

            # Annotating 3'UTR by mode + extended range
            mode.select$UTR3_hr <- ""
            gr_hr = GenomicRanges::GRanges(paste0(mode.select$seq,":", mode.select$mode.pos - ifelse(strand == 1, 0, hr), "-", mode.select$mode.pos + ifelse(strand == 1, hr, 0), ":", mode.select$strand))
            all_UTR3_hr_hits = GenomicAlignments::findOverlaps(gr_hr, my_UTR3, type = 'any')
            idx_to_annotate_UTR3_hr = S4Vectors::queryHits(all_UTR3_hr_hits)
            mode.select$UTR3_hr[unique(idx_to_annotate_UTR3_hr)] <- "YES"

            # flag isolated peak and peak clusters
            # initialize peak region
            mode.pos = NULL
            mode.select = mode.select %>% dplyr::mutate(start = mode.pos - 300, end = mode.pos + 300)
            mode.select$start = pmax(mode.select$start, gene.start)
            mode.select$end = pmin(mode.select$end, gene.end)
            mode.select$width = mode.select$end - mode.select$start

            # clustered vs isolated peaks
            mode.select$pre_end = dplyr::lag(mode.select$end, 1)
            mode.select$is.isolate = mode.select$start > mode.select$pre_end

            mode.select$cluster = 1
            if (nrow(mode.select) > 1) {
              for (i in 2:nrow(mode.select)) {
                if (!mode.select$is.isolate[i]) mode.select$cluster[i] = mode.select$cluster[i-1]
                else mode.select$cluster[i] = mode.select$cluster[i-1] + 1
              }
            }

            mode.select$cluster = as.factor(mode.select$cluster)

            # aggregate peak regions by cluster
            mode.select.agg = mode.select %>% dplyr::group_by(cluster) %>% dplyr::summarise(K = dplyr::n(), start.c = min(start), end.c = max(end))
            mode.select.agg$width.c = mode.select.agg$end.c - mode.select.agg$start.c
            mode.select = mode.select %>% dplyr::left_join(mode.select.agg, by = 'cluster')

            # add polyA_ID, cluster_ID
            polyA.ids <- paste0(mode.select$gene,":",mode.select$seq,":",mode.select$start,"-",mode.select$end,":",mode.select$strand)
            cluster.ids <- paste0(mode.select$gene, ':', mode.select$cluster)
            mode.select$polyA_ID = polyA.ids
            mode.select$cluster_ID = cluster.ids

            # label CDS helper peaks and filter out CDS/non-utr3 peaks
            keep.clusterid = unique(mode.select$cluster_ID[mode.select$UTR3_hr == 'YES'])
            mode.select = mode.select[mode.select$cluster_ID %in% keep.clusterid, ]

            if (length(keep.clusterid) > 0) {

              mode.select$helper = FALSE
              mode.select$helper[which(mode.select$UTR3_hr != 'YES')] = TRUE

              # fit peak regions
              mode.select$cluster_ID.copy = mode.select$cluster_ID
              message(paste0("Number of peak modes: ", nrow(mode.select)))
              if (nrow(mode.select) > 20) message(paste0("Large number of peak modes: ", max(table(mode.select$cluster_ID))))

              per_clusterid = split(x = mode.select[, -ncol(mode.select)], f = mode.select[, ncol(mode.select)])
              per_clusterid.fit = lapply(per_clusterid, Fit_regions, coverage.data = data, sd.init = 100)
              tbl.fit <- do.call(what = rbind, args = per_clusterid.fit)

              # discard peaks that fail to fit and abnormal from:to
              keep.ids = !is.na(tbl.fit$sigma) & (tbl.fit$from < tbl.fit$to)
              tbl.fit = tbl.fit[keep.ids,]

              if (sum(keep.ids) > 0) {

                tbl.fit = tbl.fit %>% dplyr::arrange(mode.pos)

                # aggregate peak regions by cluster
                tbl.fit[, c('K','start.c','end.c', 'width.c')] = list(NULL)
                from = to = NULL
                tbl.fit.agg = tbl.fit %>% dplyr::group_by(cluster) %>% dplyr::summarise(K = dplyr::n(), start.c = min(from), end.c = max(to))
                tbl.fit.agg$width.c = tbl.fit.agg$end.c - tbl.fit.agg $start.c
                tbl.fit = tbl.fit %>% dplyr::left_join(tbl.fit.agg, by = 'cluster')

                # update polyA_ID, cluster_ID
                polyA.ids <- paste0(tbl.fit$gene, ":", tbl.fit$seq, ":", tbl.fit$from, "-", tbl.fit$to, ":", tbl.fit$strand)
                tbl.fit$polyA_ID = polyA.ids

                # label CDS helper peaks and filter out CDS/non-utr3 peaks
                keep.clusterid.after.fit = unique(tbl.fit$cluster_ID[tbl.fit$UTR3_hr == 'YES'])
                tbl.fit = tbl.fit[tbl.fit $cluster_ID %in% keep.clusterid.after.fit, ]

                print(as.data.frame(tbl.fit))

                if (length(keep.clusterid.after.fit) > 0) {

                  tbl.fit$helper = FALSE
                  tbl.fit$helper[which(tbl.fit$UTR3_hr != 'YES')] = TRUE

                  # delete intermediate columns
                  tbl.fit[, c('sum.pi','start','end', 'width', 'pre_end', 'is.isolate', 'cluster')] = list(NULL)

                  message(paste0("Final number of peaks: ", nrow(tbl.fit)))

                  return (tbl.fit)

                } # if length(keep.clusterid.after.fit) > 0
              } # if sum(keep.ids) > 0
            } # if length(keep.clusterid) > 0
          } # if mode.select is non-null
        } # if length(aln_cov@values) > 1
      } # if aln is non-null after removing NA
    } # if raw aln is non-null
  } # foreach

  tictoc::toc()

  tictoc::tic('Final filtering and writing out files')

  # final filtering
  peak.sites <- res

  # Remove any duplicates
  polyA_ID = NULL
  peak.sites %>% dplyr::distinct(polyA_ID, .keep_all = TRUE) -> peak.sites
  n.updated.sites = nrow(peak.sites)
  message("There are ", n.updated.sites, " sites after filtering.")

  return(peak.sites)
  tictoc::toc()

}



## HELPER FUNCTION

##############################################################################

find.modes = function(dat, d, h, min.mode.cutoff = 200, min.mode.prop = 0.01, ref, gene.name) {
  n = nrow(dat)
  # smoothing
  dat$sm = sapply(1:n, gaussian_kernel_smoother, dat=dat, d=d, h=h)
  # find modes
  # apply cutoff. Too small peak modes are removed
  ref.utr = as.data.frame(ref[ref$gene_symbol == gene.name & ref$feature == 'three_prime_UTR', ])
  if (nrow(ref.utr) > 0) {
    loc = NULL
    for (i in 1:nrow(ref.utr)) {
      loc = c(loc, ref.utr$start[i]: ref.utr$end[i])
    } # loop
    dat.utr = dat[dat$pos %in% unique(loc), ]
    max.mode.val <- max(dat.utr$coverage)
  } # if

  if (nrow(ref.utr) == 0) max.mode.val <- max(dat$coverage)

  max.mode.cutoff <- max(min.mode.cutoff, min.mode.prop*max.mode.val)
  dat$sign = (dat$sm - dplyr::lag(dat$sm, 1)) >= 0
  dat$sign_change = dat$sign != dplyr::lag(dat$sign, 1)
  ids = which(dat$sign_change == T & dat$sign == F & dat$sm > max.mode.cutoff) - 1
  dat$select_mode = (1:n) %in% ids
  return (dat)
}


##############################################################################

gaussian_kernel_smoother = function(id, dat, d, h) {
  n = nrow(dat)
  d.half = (d-1)/2
  if (id <= d.half) d.half = id-1
  if (id >= n-d.half+1) d.half = n-id
  weights = exp(-((-d.half) : d.half)^2/2/h^2)
  weights = weights/sum(weights)
  return (sum(dat$coverage[(id-d.half) : (id+d.half)] * weights))
}

##############################################################################

merge_close_pk = function(t, d.cut = 50){
  n = nrow(t)
  too.close.grp = mode.pos = gene = strand = seq = gene.start = gene.end = NULL
  t$too.close.grp = 1
  if (n == 1) return (t[, c('mode.pos','gene','seq','strand', 'gene.start', 'gene.end')]) else {
    for (i in 2:n) {
      if (t$distance[i] <= d.cut) t$too.close.grp[i] = t$too.close.grp[i-1]
      else t$too.close.grp[i] = t$too.close.grp[i-1] + 1
    }
    t.agg = t %>% dplyr::group_by(too.close.grp) %>% dplyr::summarise(
      mode.pos = mean(mode.pos), gene = utils::head(gene,1), seq = utils::head(seq,1), strand = utils::head(strand,1),
      gene.start = utils::head(gene.start,1), gene.end = utils::head(gene.end,1))
    t.agg$mode.pos = as.integer(t.agg$mode.pos)
    return (t.agg[,-1])
  }
}

##############################################################################

Fit_regions = function(tbl.clus, coverage.data, sd.init=100) {
  K = nrow(tbl.clus)
  mode = tbl.clus$mode.pos
  start = tbl.clus$start.c[1]
  end = tbl.clus$end.c[1]
  gene.start = tbl.clus$gene.start[1]
  gene.end = tbl.clus$gene.end[1]
  fit.data <- data.frame(x = start:end, y = coverage.data[coverage.data$pos %in% start:end,"coverage"])

  res = NULL
  if (K == 1) { # isolated
    nls.res <- NULL
    tryCatch({nls.res <- stats::nls(y~k*exp(-1/2*(x-mu)^2/sigma^2), start=c(mu=mode, sigma=sd.init, k=max(fit.data$y)),
                                    data=fit.data, control = stats::nls.control(warnOnly = T))}, error=function(err){})

    if (!is.null(nls.res) & tryCatch({tmp = summary(nls.res)}, error = function(err){'err'}) != 'err') {
      v <- summary(nls.res)$parameters[, "Estimate"]
      v[2] = abs(v[2])
      from <- max(gene.start, floor(v[1] - 3*v[2]))
      to <- min(gene.end, floor(v[1] + 3*v[2]))
      res = cbind(tbl.clus, data.frame(mu=v[1], sigma=v[2], pi=1, from=from, to=to, width.fit=to-from))
    } else {
      res = cbind(tbl.clus, data.frame(mu=NA, sigma=NA, pi=NA, from=NA, to=NA, width.fit=NA))
    }
  } # K=1

  if (K > 1) { # clustered

    pkid.utr = which(tbl.clus$helper == FALSE)
    pkid.lb = max(1, min(pkid.utr) - 2)
    pkid.ub = min(K, max(pkid.utr) + 2)

    start = tbl.clus$start[pkid.lb]
    end = tbl.clus$end[pkid.ub]
    tbl.clus = tbl.clus[pkid.lb:pkid.ub, ]
    K = nrow(tbl.clus)
    # update cluster start / end position, width, number of pks and coverage data
    tbl.clus$start.c = start; tbl.clus$end.c = end; tbl.clus$width.c = end - start; tbl.clus$K = K
    mode = tbl.clus$mode.pos
    fit.data = data.frame(x = start:end, y = coverage.data[coverage.data$pos %in% start:end, "coverage"])

    sd.init.vec = rep(sd.init, K)
    w.init.vec = rep(1/K, K)

    if (any(tbl.clus$helper == T)) { # helper exist
      ee <- NULL
      tryCatch({ee <- mixture.EM.mu.fix(fit.data$x, fit.data$y, K, mu.true=mode, sd.init=sd.init.vec, w.init=w.init.vec, max.iter=3000, tol=1e-5)}, error=function(err){})
      if(!is.null(ee)) {
        pi_hat = ee[[1]]
        sd_hat = ee[[2]]
        from = pmax(gene.start, floor(mode - 3*sd_hat))
        to = pmin(gene.end, floor(mode + 3*sd_hat))
        res = cbind(tbl.clus, data.frame(mu=mode, sigma=sd_hat, pi=pi_hat, from=from, to=to, width.fit=to-from))
      } else {
        res = cbind(tbl.clus, data.frame(mu=NA, sigma=NA, pi=NA, from=NA, to=NA, width.fit=NA))
      }
    } else { # no helper
      ee = NULL
      tryCatch({ee <- mixture.EM(fit.data$x, fit.data$y, K, mu.init=mode, sd.init=sd.init.vec, w.init=w.init.vec, max.iter=3000, tol=1e-5)}, error=function(err){})
      if(!is.null(ee)) {
        pi_hat = ee[[1]]
        mu_hat = ee[[2]]
        sd_hat = ee[[3]]
        from = pmax(gene.start, floor(mu_hat - 3*sd_hat))
        to = pmin(gene.end, floor(mu_hat + 3*sd_hat))
        res = cbind(tbl.clus, data.frame(mu=mu_hat, sigma=sd_hat, pi=pi_hat, from=from, to=to, width.fit=to-from))
      } else {
        res = cbind(tbl.clus, data.frame(mu=NA, sigma=NA, pi=NA, from=NA, to=NA, width.fit=NA))
      }
    } # no helper
  } # K>1

  return (res)
}

################################################################################

# Gaussian mixture fit
# likelihood matrix, nrow(X) by K
compute.L = function(X, K, mu, sd) {
  L = matrix(NA, nrow=length(X), ncol= K)
  for (i in 1:K) L[, i] = stats::dnorm(X, mean=mu[i], sd = sd[i])
  return (L)
}

##############################################################################

compute.log.lik <- function(L, y, K, w) {
  for (i in 1:K) L[,i] = L[,i]*w[i]
  return(sum(log(rowSums(L))*y))
}

##############################################################################

# K need to match the dimension of parameters
mixture.EM <- function(X, y, K, mu.init, sd.init, w.init, max.iter=3000, tol=1e-5) {

  # initialize parameters and likelihood matrix
  w.curr <- w.init; mu.curr <- mu.init; sd.curr <- sd.init
  L_curr = compute.L(X, K, mu.curr, sd.curr)

  # store log-likelihoods for each iteration
  log_liks <- c()
  ll <- compute.log.lik(L_curr, y, K, w.curr)
  log_liks <- c(log_liks, ll)
  delta.ll <- 1
  iter = 1

  while(delta.ll > tol & iter <= max.iter) {
    para.curr <- EM.iter(w.curr, mu.curr, X, y, K, L_curr)
    w.curr = para.curr[[1]]
    mu.curr = para.curr[[2]]
    sd.curr = para.curr[[3]]
    L_curr = compute.L(X, K, mu.curr, sd.curr)
    ll <- compute.log.lik(L_curr, y, K, w.curr)
    log_liks <- c(log_liks, ll)
    delta.ll <- log_liks[length(log_liks)]  - log_liks[length(log_liks)-1]
    iter = iter + 1
  }

  return(list(w.curr, mu.curr, sd.curr, log_liks, iter))
}

##############################################################################

EM.iter <- function(w.curr, mu.curr, X, y, K, L) {
  # E-step: compute E_{Z|X,w0}[I(Z_i = k)] or r_Zi(k)
  z_ik <- L
  for(i in 1:K) z_ik[,i] <- w.curr[i]*z_ik[,i]
  z_ik <- z_ik / rowSums(z_ik)

  # M-step, update w, sigma
  w.next <- colSums(z_ik*y)/sum(z_ik*y)
  mu.next = colSums(t(t(X*z_ik*y)/colSums(z_ik*y)))
  xsq_mat = (t(t(matrix(rep(X,K), ncol=K)) - mu.curr))^2
  sd.next = sqrt(colSums(xsq_mat*t(t(z_ik*y)/colSums(z_ik*y))))

  return(list(w.next, mu.next, sd.next))
}

##############################################################################

mixture.EM.mu.fix <- function(X, y, K, mu.true, sd.init, w.init, max.iter=3000, tol=1e-5) {

  # initialize parameters and likelihood matrix
  w.curr <- w.init; mu.curr <- mu.true; sd.curr <- sd.init
  L_curr = compute.L(X, K, mu.curr, sd.curr)

  # store log-likelihoods for each iteration
  log_liks <- c()
  ll <- compute.log.lik(L_curr, y, K, w.curr)
  log_liks <- c(log_liks, ll)
  delta.ll <- 1
  iter = 1

  while(delta.ll > tol & iter <= max.iter) {
    para.curr <- EM.iter.mu.fix(w.curr, mu.curr, X, y, K, L_curr)
    w.curr = para.curr[[1]]
    sd.curr = para.curr[[2]]
    L_curr = compute.L(X, K, mu.curr, sd.curr)
    ll <- compute.log.lik(L_curr, y, K, w.curr)
    log_liks <- c(log_liks, ll)
    delta.ll <- log_liks[length(log_liks)]  - log_liks[length(log_liks)-1]
    iter = iter + 1
  }

  return(list(w.curr, sd.curr, log_liks, iter))
}

##############################################################################

EM.iter.mu.fix <- function(w.curr, mu.curr, X, y, K, L) {
  # E-step: compute E_{Z|X,w0}[I(Z_i = k)] or r_Zi(k)
  z_ik <- L
  for(i in 1:K) z_ik[,i] <- w.curr[i]*z_ik[,i]
  z_ik <- z_ik / rowSums(z_ik)

  # M-step, update w, sigma
  w.next <- colSums(z_ik*y)/sum(z_ik*y)
  xsq_mat = (t(t(matrix(rep(X,K), ncol=K)) - mu.curr))^2
  sd.next = sqrt(colSums(xsq_mat*t(t(z_ik*y)/colSums(z_ik*y))))

  return(list(w.next, sd.next))
}

## END ##

