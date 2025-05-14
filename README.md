
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Infernape

<!-- badges: start -->
<!-- badges: end -->

The goal of Infernape is to identify and quantify APA events from
scRNA-seq data.

Please cite the [paper](https://genome.cshlp.org/content/33/10/1774.short) "Kang, Bowei, Yalan Yang, Kaining Hu, Xiangbin Ruan, Yi-Lin Liu, Pinky Lee, Jasper Lee, Jingshu Wang, and Xiaochang Zhang. "Infernape uncovers cell type–specific and spatially resolved alternative polyadenylation in the brain." Genome Research 33, no. 10 (2023): 1774-1787."

## Installation

You can install the development version of Infernape from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kangbw702/Infernape")
```
Please install Infernape under R/4.0.0. Recommended versions of dependencies are as follows:

``` r
BiocGenerics_0.34.0
Biostrings_2.56.0
BSgenome_1.56.0
BSgenome.Mmusculus.UCSC.mm10_1.4.0
doParallel_1.0.15
foreach_1.5.0
dplyr_0.8.5
tidyr_1.1.0 
GenomicAlignments_1.24.0
GenomicRanges_1.40.0
GenomeInfoDb_1.24.0
ggplot2_3.3.1
IRanges_2.22.2
magrittr_1.5
MGLM_0.2.1
Matrix_1.4-1
plyr_1.8.6
Rsamtools_2.4.0
S4Vectors_0.26.1
tictoc_1.2
```

## Example

This is a basic example which shows you how to solve a common problem.
The data used for this example can be downloaded from [here](https://www.dropbox.com/sh/pp9hoe128lfci7u/AABCtyOjxB8Ejb_ObcBw7k9ya?dl=0).

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
```

## References

3'UTR and PAS references for Human, Rat, and Caenorhabditis_elegans can be found [here](https://www.dropbox.com/sh/3e5kwslflzdevfu/AAD5nyd-Bcf9GOSFw8I7YwhGa?dl=0). Credit to Yalan Yang.



<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
