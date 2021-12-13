Compute per-barcode escape
================
Tyler Starr
11/20/2020

-   [Setup](#setup)
-   [Calculating per-barcode escape
    score](#calculating-per-barcode-escape-score)
-   [Data Output](#data-output)

This notebook computes and summarizes per-barcode escape fractions for
the homolog library variants.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$escape_scores_dir)){
  dir.create(file.path(config$escape_scores_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [5] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0       tibble_3.0.2     
    ##  [9] ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8 yaml_2.2.0       
    ## [13] knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.5     compiler_3.6.2  
    ##  [5] dbplyr_1.4.2     tools_3.6.2      digest_0.6.23    lubridate_1.7.4 
    ##  [9] jsonlite_1.6     evaluate_0.14    lifecycle_0.2.0  gtable_0.3.0    
    ## [13] pkgconfig_2.0.3  rlang_0.4.7      reprex_0.3.0     cli_2.0.0       
    ## [17] rstudioapi_0.10  DBI_1.1.0        haven_2.2.0      xfun_0.11       
    ## [21] withr_2.1.2      xml2_1.2.2       httr_1.4.1       fs_1.3.1        
    ## [25] hms_0.5.2        generics_0.0.2   vctrs_0.3.1      grid_3.6.2      
    ## [29] tidyselect_1.1.0 glue_1.3.1       R6_2.4.1         fansi_0.4.0     
    ## [33] readxl_1.3.1     rmarkdown_2.0    modelr_0.1.5     magrittr_1.5    
    ## [37] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   htmltools_0.4.0 
    ## [41] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3   
    ## [45] munsell_0.5.0    broom_0.7.0      crayon_1.3.4

## Setup

Read in table of variant genotypes and barcode counts in each sequencing
sample, barcode-variant association table, and previously determined
barcode-level expression measurements. Add variant sequence and
expression information to the sequencing counts data table.

``` r
dt <- data.table(read.csv(file=config$variant_counts,stringsAsFactors = F))

bc_lookup <- data.table(read.csv(file=config$barcode_variant_table))

bc_expr <- data.table(read.csv(file=config$barcode_expression))

#add target and expression to dt table
dt <- merge(dt[,.(barcode,library,sample,count)],bc_lookup[,.(barcode,library,target)],by=c("library","barcode"))

dt <- merge(dt[,.(barcode,library,sample,count,target)],bc_expr[,.(barcode,library,ML_meanF)],by=c("library","barcode"))

#read dataframe with list of barcode runs and associated metainfo for calculating escape fracs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))
#rename frac_escape to Frac, to avoid confusion with escape_frac later on
barcode_runs$Frac <- barcode_runs$frac_escape

#add this info to the dt table
dt <- merge(dt[,.(library,target,barcode,ML_meanF,sample,count)],barcode_runs[,c("date","antibody","selection","sample","Frac","cells_sorted")],by=c("sample"))
```

We have some barcodes in this library that do not express – we expect
these to impact our escape fraction calculation, which relies on the
fraction of RBD+ cells that escape antibody – however, these
non-expressing cells will be in the pre-count (because our library was
not pre-RBD+ sorted), even though they are not in the RBD+ subset that
are sampled for escape. Therefore, they can bias up escape fraction
scores (in a systematic way, but still, makes it ugly that escape
fractions are biased to be larger than the theoretical limit 1). Because
they are not truly part of the “escape” experiment (filtered out by RBD+
gating), we remove them.

``` r
dt <- dt[ML_meanF > 7.5 | is.na(ML_meanF),]
```

## Calculating per-barcode escape score

Next, for each barcode at each of the ACE2 concentrations, calculate the
“escape fraction”, corresponding to the estimated fraction of cells with
that barcode that were sorted into the antibody-escape bin. Escape
fraciton for a variant (E\_v) is equal to its post-sort frequency
divided by its pre-sort frequency, multiplied by the big-F fraction,
corresponding to the fraction of cells that were sorted into the
antibody-escape bin.

``` r
#add big-N N_total_count, which corresponds to the sum of counts across all barcodes in a given sample (either pre or post selection)
dt[,N_total_count:=sum(count),by=c("library","sample")]

#post_freq for escape samples
dt[selection=="escape",post_freq:=count/N_total_count,by=c("library","sample","barcode")]

#pre_freq for reference samples
dt[selection=="reference",pre_freq:=count/N_total_count,by=c("library","sample","barcode")]

#for each escape sample, get the barcode pre_freq and pre_count from the corresponding date's reference sample
#define a function for vectorizing to pull out the associated pre-sort frequency
get_presort_freq_count <- function(day, lib, bc){
  return(list(dt[selection=="reference" & date==day & library==lib & barcode==bc,pre_freq],dt[selection=="reference" & date==day & library==lib & barcode==bc,count]))
}

#assign
dt[selection=="escape",c("pre_freq","pre_count"):=get_presort_freq_count(day=date,lib=library,bc=barcode),by=c("library","sample","barcode")]

#calculate escape frac
dt[selection=="escape",escape_frac:=Frac*post_freq/pre_freq]
```

These calculations took a while to compute, so let’s output the summary
table, and do post-calc analysis in a new notebook (so this doesn’t have
to be run every time we want to replot or change filtering)

## Data Output

``` r
dt <- dt[selection=="escape",]
dt[,post_count:=count]
dt[,expression:=ML_meanF]
dt[,.(library, target, barcode, expression, antibody, pre_count, post_count, pre_freq, post_freq, escape_frac)] %>%
  mutate_if(is.numeric, round, digits=5) %>%
  write.csv(file=config$escape_fracs_barcodes, row.names=F)
```
