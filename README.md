
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Infernape

<!-- badges: start -->
<!-- badges: end -->

The goal of Infernape is to identify and quantify APA events from
scRNA-seq data.

## Installation

You can install the development version of Infernape from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kangbw702/Infernape")
```

## Example

This is a basic example which shows you how to solve a common problem.
The data used for this example can be downloaded from
<https://www.dropbox.com/sh/pp9hoe128lfci7u/AABCtyOjxB8Ejb_ObcBw7k9ya?dl=0>.

Function `Infernape_cnt` outputs raw peaks, peak annotation (before and
after filtering), and peak by cell UMI count matrix.

``` r
library(Infernape)
genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
Infernape_cnt(genome.ref = '../data/ref.csv',
              bam = '../data/TEGLU16.bam',
              batch.start = 2901,
              batch.end = 2920,
              ncores = 1,
              d = 31,
              h = 5,
              d.cut = 50,
              hr = 160,
              min.mode.prop = 0.05,
              min.mode.cutoff = 5,
              output.path = '../result',
              pas.reference.file = '../data/PAS_withinfo.csv',
              genome = genome,
              pas.search.cut.1 = 0,
              pas.search.cut.2 = 300,
              polystretch_length = 13,
              max_mismatch = 1,
              motif.search.cut = 300,
              invert_strand = FALSE,
              q = c(110, 200),
              whitelist.file = "../data/whitelist.TEGLU16.csv",
              start.cid = NULL,
              end.cid = NULL
)
#> Total number of genes: 21544
#> Genes to process in this batch: from 2901 to 2920
#> Import genome reference ...: 2.308 sec elapsed
#> Warning: package 'GenomicRanges' was built under R version 4.1.2
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> Warning: package 'S4Vectors' was built under R version 4.1.3
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: GenomeInfoDb
#> Warning: package 'GenomeInfoDb' was built under R version 4.1.2
#> 
#> 
#> 2901th gene Cdc23 chr18 34630947 34651735 -1
#> Number of rows in raw coverage: 20789
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>   mode.pos  gene   seq strand gene.start gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 34632942 Cdc23 chr18      -   34630947 34651735       YES              YES
#>      start      end width pre_end is.isolate cluster
#> 1 34632642 34633242   600      NA         NA       1
#>                          polyA_ID cluster_ID helper       mu    sigma pi
#> 1 Cdc23:chr18:34632770-34633099:-    Cdc23:1  FALSE 34632935 54.77055  1
#>       from       to width.fit K  start.c    end.c width.c
#> 1 34632770 34633099       329 1 34632770 34633099     329
#> Final number of peaks: 1
#> 
#> 
#> 2902th gene Cdc25a chr9 109875579 109893895 1
#> Number of rows in raw coverage: 18317
#> 
#> 
#> 2903th gene Cdc25b chr2 131186958 131199183 1
#> Number of rows in raw coverage: 12226
#> 
#> 
#> 2904th gene Cdc25c chr18 34732993 34751533 -1
#> Number of rows in raw coverage: 18541
#> 
#> 
#> 2905th gene Cdc26 chr4 62392346 62408642 -1
#> Number of rows in raw coverage: 16297
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 2
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>   mode.pos  gene  seq strand gene.start gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 62393966 Cdc26 chr4      -   62392346 62408642       YES              YES
#> 2 62394747 Cdc26 chr4      -   62392346 62408642       YES              YES
#>      start      end width  pre_end is.isolate cluster
#> 1 62393666 62394266   600       NA         NA       1
#> 2 62394447 62395047   600 62394266       TRUE       2
#>                         polyA_ID cluster_ID helper       mu    sigma pi
#> 1 Cdc26:chr4:62393814-62394126:-    Cdc26:1  FALSE 62393971 52.10829  1
#> 2 Cdc26:chr4:62394584-62394902:-    Cdc26:2  FALSE 62394743 52.95152  1
#>       from       to width.fit K  start.c    end.c width.c
#> 1 62393814 62394126       312 1 62393814 62394126     312
#> 2 62394584 62394902       318 1 62394584 62394902     318
#> Final number of peaks: 2
#> 
#> 
#> 2906th gene Cdc27 chr11 104502513 104550483 -1
#> Number of rows in raw coverage: 47971
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>    mode.pos  gene   seq strand gene.start  gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 104502919 Cdc27 chr11      -  104502513 104550483       YES              YES
#>       start       end width pre_end is.isolate cluster
#> 1 104502619 104503219   600      NA         NA       1
#>                            polyA_ID cluster_ID helper        mu   sigma pi
#> 1 Cdc27:chr11:104502513-104503349:-    Cdc27:1  FALSE 104502892 152.596  1
#>        from        to width.fit K   start.c     end.c width.c
#> 1 104502513 104503349       836 1 104502513 104503349     836
#> Final number of peaks: 1
#> 
#> 
#> 2907th gene Cdc34 chr10 79682195 79688939 1
#> Number of rows in raw coverage: 6745
#> The following `from` values were not present in `x`: -1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>   mode.pos  gene   seq strand gene.start gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 79688222 Cdc34 chr10      +   79682195 79688939       YES              YES
#>      start      end width pre_end is.isolate cluster
#> 1 79687922 79688522   600      NA         NA       1
#>                          polyA_ID cluster_ID helper       mu    sigma pi
#> 1 Cdc34:chr10:79688071-79688416:+    Cdc34:1  FALSE 79688244 57.41739  1
#>       from       to width.fit K  start.c    end.c width.c
#> 1 79688071 79688416       345 1 79688071 79688416     345
#> Final number of peaks: 1
#> 
#> 
#> 2908th gene Cdc34b chr11 94741782 94743033 1
#> Number of rows in raw coverage: 1252
#> 
#> 
#> 2909th gene Cdc37 chr9 21133204 21149982 -1
#> Number of rows in raw coverage: 16779
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>   mode.pos  gene  seq strand gene.start gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 21139566 Cdc37 chr9      -   21133204 21149982       YES              YES
#>      start      end width pre_end is.isolate cluster
#> 1 21139266 21139866   600      NA         NA       1
#>                         polyA_ID cluster_ID helper       mu   sigma pi     from
#> 1 Cdc37:chr9:21139434-21139714:-    Cdc37:1  FALSE 21139575 46.6527  1 21139434
#>         to width.fit K  start.c    end.c width.c
#> 1 21139714       280 1 21139434 21139714     280
#> Final number of peaks: 1
#> 
#> 
#> 2910th gene Cdc37l1 chr19 28990352 29026681 1
#> Number of rows in raw coverage: 36330
#> The following `from` values were not present in `x`: -1
#> Number of peak modes: 2
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used

#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>   mode.pos    gene   seq strand gene.start gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 29014961 Cdc37l1 chr19      +   28990352 29026681       YES              YES
#> 2 29026497 Cdc37l1 chr19      +   28990352 29026681       YES              YES
#>      start      end width  pre_end is.isolate cluster
#> 1 29014661 29015261   600 28995471       TRUE       2
#> 2 29026197 29026681   484 29015261       TRUE       3
#>                            polyA_ID cluster_ID helper       mu    sigma pi
#> 1 Cdc37l1:chr19:29014815-29015107:+  Cdc37l1:2  FALSE 29014961 48.63785  1
#> 2 Cdc37l1:chr19:29026299-29026681:+  Cdc37l1:3  FALSE 29026496 65.70711  1
#>       from       to width.fit K  start.c    end.c width.c
#> 1 29014815 29015107       292 1 29014815 29015107     292
#> 2 29026299 29026681       382 1 29026299 29026681     382
#> Final number of peaks: 2
#> 
#> 
#> 2911th gene Cdc40 chr10 40831621 40883311 -1
#> Number of rows in raw coverage: 51691
#> The following `from` values were not present in `x`: 1
#> 
#> 
#> 2912th gene Cdc42 chr4 137318643 137357720 -1
#> Number of rows in raw coverage: 39078
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 3
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used

#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#> Warning in stats::nls(y ~ k * exp(-1/2 * (x - mu)^2/sigma^2), start = c(mu =
#> mode, : step factor 0.000488281 reduced below 'minFactor' of 0.000976562
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>    mode.pos  gene  seq strand gene.start  gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 137319836 Cdc42 chr4      -  137318643 137357720       YES              YES
#> 2 137321931 Cdc42 chr4      -  137318643 137357720       YES              YES
#> 3 137322579 Cdc42 chr4      -  137318643 137357720                YES     YES
#>       start       end width   pre_end is.isolate cluster
#> 1 137319536 137320136   600        NA         NA       1
#> 2 137321631 137322231   600 137320136       TRUE       2
#> 3 137322279 137322879   600 137322231       TRUE       3
#>                           polyA_ID cluster_ID helper        mu    sigma pi
#> 1 Cdc42:chr4:137319669-137320030:-    Cdc42:1  FALSE 137319850 60.17544  1
#> 2 Cdc42:chr4:137321723-137322107:-    Cdc42:2  FALSE 137321916 63.95951  1
#> 3 Cdc42:chr4:137322429-137322701:-    Cdc42:3  FALSE 137322566 45.35001  1
#>        from        to width.fit K   start.c     end.c width.c
#> 1 137319669 137320030       361 1 137319669 137320030     361
#> 2 137321723 137322107       384 1 137321723 137322107     384
#> 3 137322429 137322701       272 1 137322429 137322701     272
#> Final number of peaks: 3
#> 
#> 
#> 2913th gene Cdc42bpa chr1 179960585 180165603 1
#> Number of rows in raw coverage: 205019
#> The following `from` values were not present in `x`: -1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>    mode.pos     gene  seq strand gene.start  gene.end UTR3_mode CDS_mode
#> 1 180165439 Cdc42bpa chr1      +  179960585 180165603       YES         
#>   UTR3_hr     start       end width   pre_end is.isolate cluster
#> 1     YES 180165139 180165603   464 180128758       TRUE       3
#>                              polyA_ID cluster_ID helper        mu    sigma pi
#> 1 Cdc42bpa:chr1:180165218-180165603:+ Cdc42bpa:3  FALSE 180165440 73.97914  1
#>        from        to width.fit K   start.c     end.c width.c
#> 1 180165218 180165603       385 1 180165218 180165603     385
#> Final number of peaks: 1
#> 
#> 
#> 2914th gene Cdc42bpb chr12 111292972 111377718 -1
#> Number of rows in raw coverage: 84747
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>    mode.pos     gene   seq strand gene.start  gene.end UTR3_mode CDS_mode
#> 1 111293153 Cdc42bpb chr12      -  111292972 111377718       YES         
#>   UTR3_hr     start       end width pre_end is.isolate cluster
#> 1     YES 111292972 111293453   481      NA         NA       1
#>                               polyA_ID cluster_ID helper        mu    sigma pi
#> 1 Cdc42bpb:chr12:111293004-111293294:- Cdc42bpb:1  FALSE 111293150 48.35535  1
#>        from        to width.fit K   start.c     end.c width.c
#> 1 111293004 111293294       290 1 111293004 111293294     290
#> Final number of peaks: 1
#> 
#> 
#> 2915th gene Cdc42bpg chr19 6306456 6325999 1
#> Number of rows in raw coverage: 19544
#> 
#> 
#> 2916th gene Cdc42ep1 chr15 78842624 78850902 1
#> 
#> 
#> 2917th gene Cdc42ep2 chr19 5914956 5924816 -1
#> Number of rows in raw coverage: 9861
#> 
#> 
#> 2918th gene Cdc42ep3 chr17 79331296 79355091 -1
#> Number of rows in raw coverage: 23796
#> The following `from` values were not present in `x`: 1
#> Number of peak modes: 1
#> Warning in if (!is.null(nls.res) & tryCatch({: the condition has length > 1 and
#> only the first element will be used
#>   mode.pos     gene   seq strand gene.start gene.end UTR3_mode CDS_mode UTR3_hr
#> 1 79334606 Cdc42ep3 chr17      -   79331296 79355091       YES              YES
#>      start      end width pre_end is.isolate cluster
#> 1 79334306 79334906   600      NA         NA       1
#>                             polyA_ID cluster_ID helper       mu    sigma pi
#> 1 Cdc42ep3:chr17:79334472-79334775:- Cdc42ep3:1  FALSE 79334624 50.58771  1
#>       from       to width.fit K  start.c    end.c width.c
#> 1 79334472 79334775       303 1 79334472 79334775     303
#> Final number of peaks: 1
#> 
#> 
#> 2919th gene Cdc42ep4 chr11 113726850 113751318 -1
#> Number of rows in raw coverage: 24469
#> 
#> 
#> 2920th gene Cdc42ep5 chr7 4151260 4164860 -1
#> Number of rows in raw coverage: 13601
#> Detecting peaks ... : 14.41 sec elapsed
#> There are 14 sites after filtering.
#> 
#> 
#> Number of raw peaks identified: 14
#> 
#> 
#> Peak table is output? TRUE
#>    mode.pos  gene   seq strand gene.start  gene.end UTR3_mode CDS_mode UTR3_hr
#> 1  34632942 Cdc23 chr18      -   34630947  34651735       YES              YES
#> 2  62393966 Cdc26  chr4      -   62392346  62408642       YES              YES
#> 3  62394747 Cdc26  chr4      -   62392346  62408642       YES              YES
#> 4 104502919 Cdc27 chr11      -  104502513 104550483       YES              YES
#> 5  79688222 Cdc34 chr10      +   79682195  79688939       YES              YES
#> 6  21139566 Cdc37  chr9      -   21133204  21149982       YES              YES
#>                            polyA_ID cluster_ID helper        mu     sigma pi
#> 1   Cdc23:chr18:34632770-34633099:-    Cdc23:1  FALSE  34632935  54.77055  1
#> 2    Cdc26:chr4:62393814-62394126:-    Cdc26:1  FALSE  62393971  52.10829  1
#> 3    Cdc26:chr4:62394584-62394902:-    Cdc26:2  FALSE  62394743  52.95152  1
#> 4 Cdc27:chr11:104502513-104503349:-    Cdc27:1  FALSE 104502892 152.59596  1
#> 5   Cdc34:chr10:79688071-79688416:+    Cdc34:1  FALSE  79688244  57.41739  1
#> 6    Cdc37:chr9:21139434-21139714:-    Cdc37:1  FALSE  21139575  46.65270  1
#>        from        to width.fit K   start.c     end.c width.c
#> 1  34632770  34633099       329 1  34632770  34633099     329
#> 2  62393814  62394126       312 1  62393814  62394126     312
#> 3  62394584  62394902       318 1  62394584  62394902     318
#> 4 104502513 104503349       836 1 104502513 104503349     836
#> 5  79688071  79688416       345 1  79688071  79688416     345
#> 6  21139434  21139714       280 1  21139434  21139714     280
#> Number of peaks after pre-filtering: 14
#> Annotating known PAS around peaks.
#> Annotating motifs around peaks.
#>   |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   7%  |                                                                              |==========                                                            |  14%  |                                                                              |===============                                                       |  21%  |                                                                              |====================                                                  |  29%  |                                                                              |=========================                                             |  36%  |                                                                              |==============================                                        |  43%  |                                                                              |===================================                                   |  50%  |                                                                              |========================================                              |  57%  |                                                                              |=============================================                         |  64%  |                                                                              |==================================================                    |  71%  |                                                                              |=======================================================               |  79%  |                                                                              |============================================================          |  86%  |                                                                              |=================================================================     |  93%  |                                                                              |======================================================================| 100%
#> 
#> 
#> Peak annotation completes.
#>    mode.pos  gene   seq strand gene.start  gene.end        mu     sigma pi
#> 1  34632942 Cdc23 chr18      -   34630947  34651735  34632935  54.77055  1
#> 2  62393966 Cdc26  chr4      -   62392346  62408642  62393971  52.10829  1
#> 3  62394747 Cdc26  chr4      -   62392346  62408642  62394743  52.95152  1
#> 4 104502919 Cdc27 chr11      -  104502513 104550483 104502892 152.59596  1
#> 5  79688222 Cdc34 chr10      +   79682195  79688939  79688244  57.41739  1
#> 6  21139566 Cdc37  chr9      -   21133204  21149982  21139575  46.65270  1
#>        from        to width.fit                          polyA_ID cluster_ID
#> 1  34632770  34633099       329   Cdc23:chr18:34632770-34633099:-    Cdc23:1
#> 2  62393814  62394126       312    Cdc26:chr4:62393814-62394126:-    Cdc26:1
#> 3  62394584  62394902       318    Cdc26:chr4:62394584-62394902:-    Cdc26:2
#> 4 104502513 104503349       836 Cdc27:chr11:104502513-104503349:-    Cdc27:1
#> 5  79688071  79688416       345   Cdc34:chr10:79688071-79688416:+    Cdc34:1
#> 6  21139434  21139714       280    Cdc37:chr9:21139434-21139714:-    Cdc37:1
#>   lb.center ub.center           PAS_full n_PAS crude.start crude.end
#> 1  34632642  34632942 233,167,153,128,27     5    34632642  34633242
#> 2  62393666  62393966                160     1    62393666  62394266
#> 3  62394447  62394747        178,164,113     3    62394447  62395047
#> 4 104502619 104502919     206,173,159,51     4   104502619 104503219
#> 5  79688222  79688522   34,53,95,126,174     5    79687922  79688522
#> 6  21139266  21139566  234,134,102,77,17     5    21139266  21139866
#>   pA_motif_TTTAAA pA_motif_TATAAA pA_motif_GATAAA pA_motif_CATAAA
#> 1            <NA>            <NA>            <NA>             144
#> 2            <NA>            <NA>            <NA>            <NA>
#> 3             137            <NA>            <NA>            <NA>
#> 4            <NA>            <NA>            <NA>            <NA>
#> 5            <NA>            <NA>            <NA>            <NA>
#> 6            <NA>            <NA>            <NA>            <NA>
#>   pA_motif_ATTAAA pA_motif_AGTAAA pA_motif_ACTAAA pA_motif_AATGAA
#> 1            <NA>            <NA>              NA              NA
#> 2            <NA>            <NA>              NA              NA
#> 3            <NA>            <NA>              NA              NA
#> 4             133            <NA>              NA              NA
#> 5            <NA>            <NA>              NA              NA
#> 6            <NA>            <NA>              NA              NA
#>   pA_motif_AATATA pA_motif_AATAAA pA_motif_AAGAAA pA_motif_AACAAA pA_stretch
#> 1            <NA>             136            <NA>            <NA>        183
#> 2            <NA>             136            <NA>             215       <NA>
#> 3            <NA>             143            <NA>            <NA>       <NA>
#> 4             117             155            <NA>            <NA>       <NA>
#> 5            <NA>         152,157            <NA>            <NA>       <NA>
#> 6            <NA>             114            <NA>            <NA>       <NA>
#> n.PAS: 13
#> n.pA_motif_TTTAAA: 4
#> n.pA_motif_TATAAA: 1
#> n.pA_motif_GATAAA: 0
#> n.pA_motif_CATAAA: 1
#> n.pA_motif_ATTAAA: 5
#> n.pA_motif_AGTAAA: 2
#> n.pA_motif_ACTAAA: 0
#> n.pA_motif_AATGAA: 0
#> n.pA_motif_AATATA: 2
#> n.pA_motif_AATAAA: 9
#> n.pA_motif_AAGAAA: 2
#> n.pA_motif_AACAAA: 4
#> n.motif: 14
#> n.nonA_stretch: 11
#> n.novel: 0
#> n.keep: 13
#> n.genes: 9
#>    mode.pos  gene   seq strand gene.start  gene.end        mu     sigma pi
#> 1  34632942 Cdc23 chr18      -   34630947  34651735  34632935  54.77055  1
#> 2  62393966 Cdc26  chr4      -   62392346  62408642  62393971  52.10829  1
#> 3  62394747 Cdc26  chr4      -   62392346  62408642  62394743  52.95152  1
#> 4 104502919 Cdc27 chr11      -  104502513 104550483 104502892 152.59596  1
#> 5  79688222 Cdc34 chr10      +   79682195  79688939  79688244  57.41739  1
#> 6  21139566 Cdc37  chr9      -   21133204  21149982  21139575  46.65270  1
#>        from        to width.fit                          polyA_ID cluster_ID
#> 1  34632770  34633099       329   Cdc23:chr18:34632770-34633099:-    Cdc23:1
#> 2  62393814  62394126       312    Cdc26:chr4:62393814-62394126:-    Cdc26:1
#> 3  62394584  62394902       318    Cdc26:chr4:62394584-62394902:-    Cdc26:2
#> 4 104502513 104503349       836 Cdc27:chr11:104502513-104503349:-    Cdc27:1
#> 5  79688071  79688416       345   Cdc34:chr10:79688071-79688416:+    Cdc34:1
#> 6  21139434  21139714       280    Cdc37:chr9:21139434-21139714:-    Cdc37:1
#>   lb.center ub.center           PAS_full n_PAS crude.start crude.end
#> 1  34632642  34632942 233,167,153,128,27     5    34632642  34633242
#> 2  62393666  62393966                160     1    62393666  62394266
#> 3  62394447  62394747        178,164,113     3    62394447  62395047
#> 4 104502619 104502919     206,173,159,51     4   104502619 104503219
#> 5  79688222  79688522   34,53,95,126,174     5    79687922  79688522
#> 6  21139266  21139566  234,134,102,77,17     5    21139266  21139866
#>   pA_motif_TTTAAA pA_motif_TATAAA pA_motif_GATAAA pA_motif_CATAAA
#> 1              NA              NA              NA             144
#> 2              NA              NA              NA              NA
#> 3             137              NA              NA              NA
#> 4              NA              NA              NA              NA
#> 5              NA              NA              NA              NA
#> 6              NA              NA              NA              NA
#>   pA_motif_ATTAAA pA_motif_AGTAAA pA_motif_ACTAAA pA_motif_AATGAA
#> 1              NA              NA              NA              NA
#> 2              NA              NA              NA              NA
#> 3              NA              NA              NA              NA
#> 4             133              NA              NA              NA
#> 5              NA              NA              NA              NA
#> 6              NA              NA              NA              NA
#>   pA_motif_AATATA pA_motif_AATAAA pA_motif_AAGAAA pA_motif_AACAAA pA_stretch
#> 1              NA             136            <NA>            <NA>        183
#> 2              NA             136            <NA>             215       <NA>
#> 3              NA             143            <NA>            <NA>       <NA>
#> 4             117             155            <NA>            <NA>       <NA>
#> 5              NA         152,157            <NA>            <NA>       <NA>
#> 6              NA             114            <NA>            <NA>       <NA>
#>    PAS TTTAAA TATAAA GATAAA CATAAA ATTAAA AGTAAA ACTAAA AATGAA AATATA AATAAA
#> 1 TRUE     NA     NA     NA   TRUE     NA     NA     NA     NA     NA   TRUE
#> 2 TRUE     NA     NA     NA     NA     NA     NA     NA     NA     NA   TRUE
#> 3 TRUE   TRUE     NA     NA     NA     NA     NA     NA     NA     NA   TRUE
#> 4 TRUE     NA     NA     NA     NA   TRUE     NA     NA     NA   TRUE   TRUE
#> 5 TRUE     NA     NA     NA     NA     NA     NA     NA     NA     NA   TRUE
#> 6 TRUE     NA     NA     NA     NA     NA     NA     NA     NA     NA   TRUE
#>   AAGAAA AACAAA motif   As novel select
#> 1     NA     NA  TRUE TRUE FALSE   TRUE
#> 2     NA  FALSE  TRUE   NA FALSE   TRUE
#> 3     NA     NA  TRUE   NA FALSE   TRUE
#> 4     NA     NA  TRUE   NA FALSE   TRUE
#> 5     NA     NA  TRUE   NA FALSE   TRUE
#> 6     NA     NA  TRUE   NA FALSE   TRUE
#> There are 61 whitelist barcodes.
#> Total - n.genes: 9; n.peak-clusters: 13; n.peaks: 13
#> This run - n.peak-clusters: 13; n.peaks: 13
#> Warning: package 'magrittr' was built under R version 4.1.2
#> 1th cluster: Cdc23:1; seq:chr18; strand: -1; n.pks:1
#> 2th cluster: Cdc26:1; seq:chr4; strand: -1; n.pks:1
#> 3th cluster: Cdc26:2; seq:chr4; strand: -1; n.pks:1
#> 4th cluster: Cdc27:1; seq:chr11; strand: -1; n.pks:1
#> 5th cluster: Cdc34:1; seq:chr10; strand: 1; n.pks:1
#> 6th cluster: Cdc37:1; seq:chr9; strand: -1; n.pks:1
#> 7th cluster: Cdc37l1:2; seq:chr19; strand: 1; n.pks:1
#> 8th cluster: Cdc37l1:3; seq:chr19; strand: 1; n.pks:1
#> 9th cluster: Cdc42:1; seq:chr4; strand: -1; n.pks:1
#> 10th cluster: Cdc42:2; seq:chr4; strand: -1; n.pks:1
#> 11th cluster: Cdc42:3; seq:chr4; strand: -1; n.pks:1
#> 12th cluster: Cdc42bpa:3; seq:chr1; strand: 1; n.pks:1
#> 13th cluster: Cdc42bpb:1; seq:chr12; strand: -1; n.pks:1
#> [1] "DONE!"
```

Function `Infernape_apa` performs hierarchical differential APA testing.

``` r
Infernape_apa(counts.dir = '../result/cnt_mat',
              attr.file = '../data/attr.tbl.example.csv',
              anno.file = '../result/anno_filtered.csv',
              utr3.file = '../data/ref.utr3.anno.csv',
              ctype.colname = 'ctype',
              base_grp = 'c1',
              alt_grp = 'c2',
              cut.low.pct = 0.05,
              cut.pval = 0.05,
              cut.MPRO = 0.2,
              test.type = 'gene',
              out.dir = '../result'
)
#> Peak site names in anno table are consistent with count matrix? TRUE
#> Cell barcodes in attr table are consistent with count matrix? TRUE
#> The number of multi-peak genes to test: 3
#> Warning in stats::chisq.test(cnt.tbl): Chi-squared approximation may be
#> incorrect
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 8 rows are removed because the row
#> sums are 0.
#> Warning in DMD.DM.fit(data = data, init = init, weight = weight, epsilon =
#> epsilon, : The algorithm doesn't converge within 3 iterations.
#> Warning in sqrt(diag(invI)): NaNs produced
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 21 rows are removed because the row
#> sums are 0.
#> Warning in stats::chisq.test(cnt.tbl): Chi-squared approximation may be
#> incorrect
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 14 rows are removed because the row
#> sums are 0.
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 28 rows are removed because the row
#> sums are 0.
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 3 rows are removed because the row
#> sums are 0.
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 8 rows are removed because the row
#> sums are 0.
#> Warning in MGLM::MGLMfit(dat, dist = "DM"): 11 rows are removed because the row
#> sums are 0.
#> The number of rows in the result table: 3
#> The number of significant genes: 0
#> long.vs.short: 0,0
#> Write out done!
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
