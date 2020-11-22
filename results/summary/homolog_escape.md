Compute per-homolog escape
================
Tyler Starr
11/20/2020

-   [Data input](#data-input)
-   [Filtering](#filtering)
-   [Escape fraction per homolog](#escape-fraction-per-homolog)
-   [Heatmaps](#heatmaps)
-   [Output](#output)

This notebook filters the per-barcode escape fraction estimates, and
computes summary escape fractions for each homolog.

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

Session info for reproducing environment:

    sessionInfo()

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
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

Data input
----------

Read in table of escape scores for barcodes, also annotated with variant
counts and expression that we’ll use in filtering. Remove barcodes with
pre\_count of 0 or NA expression effects – these are likely barcodes
that either do not express well, are generally low-frequency, or were
actually “mutant” variants (but no mutation in the SSM) which were not
part of the sub-pool of the homologs experiments that we used for these
antibody selections.

    dt <- data.table(read.csv(file=config$escape_fracs_barcodes,stringsAsFactors = F))
    dt <- dt[pre_count>0 & !is.na(expression),]

Filtering
---------

We will use the per-barcode pre\_count and expression scores to filter
out escape scores used in computing per-homolog escape.

First, let’s look at the distribution of pre-counts across barcodes. The
median pre-count is 203. Vertical lines on the two plots below indicate
a threshold for pre\_count of 1/2 that of the median pre-count, which is
what we’ll apply below.

    hist(log10(dt$pre_count),col="gray50",main="",xlab="pre-selection sequencing count (log10 scale)");abline(v=log10(0.5*median(dt$pre_count)),col="red",lty=2)

<img src="homolog_escape_files/figure-gfm/hist_pre_count-1.png" style="display: block; margin: auto;" />
Let’s also see how escape fraction correlates with pre-count.
Theoretically, escape\_fraction should not exceed 1, though we see it
does, particularly when pre\_count is lower.

    plot(dt$pre_count,dt$escape_frac,pch=16,col="#00000060",log="x",xlab="pre-count (log scale)",ylab="escape fraction (shouldn't exceed 1)")
    abline(h=1,lty=2,col="red")
    abline(v=0.5*median(dt$pre_count),lty=2,col="red")

<img src="homolog_escape_files/figure-gfm/scatter_pre_count-1.png" style="display: block; margin: auto;" />

Remove barcode measurements for those where pre-count is less than half
the median pre-count. This corresponds to removing variants with less
than 101.5 pre-sort counts.

    dt <- dt[pre_count > 0.5*median(dt$pre_count),]

Censor NA measurements for AncSARS1a\_alt and HKU3-8 which were
generally non-expressing from our prior measurements (and so any
barcodes that eke through are probably poorly expressed/artefactual)

    dt[target%in%c("AncSARS1a_alt","HKU3-8"),escape_frac:=NA]
    dt[target%in%c("AncSARS1a_alt","HKU3-8"),expression:=NA]

Escape fraction per homolog
---------------------------

Next, let’s visualize the escape fraction scores per barcode grouped by
variant, with violin plots.

    #set factor order for homologs to display
    dt$target <- factor(dt$target, levels=config$targets_ordered)

    ggplot(dt,aes(x=target,y=escape_frac))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("escape frac by homolog")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      facet_wrap(~antibody,ncol=1)

    ## Warning: Removed 27 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 27 rows containing non-finite values (stat_summary).

<img src="homolog_escape_files/figure-gfm/escape_frac_vioplot-1.png" style="display: block; margin: auto;" />

Collapse each homolog escape fraction to its median across barcodes.

    dt[,median_escape_frac:=median(escape_frac,na.rm=T),by=c("library","target","antibody")] #median escape across all barcodes for a homolog
    dt[,n_barcodes:=sum(!is.na(escape_frac)),by=c("library","target","antibody")] #the number of barcodes on which a homolog median escape was calculated
    dt[,median_expression:=median(expression,na.rm=T),by=c("library","target")] #the average expression of that homolog
    dt_collapse <- unique(dt[,.(library,target,median_expression,antibody,median_escape_frac,n_barcodes)]) #collapse down to homolog-level data instead of barcode-level
    dt_collapse[,expression:=median_expression];dt_collapse[,escape_frac:=median_escape_frac] #rename median quantities to just now be the homolog-level quantity
    dt_collapse <- dt_collapse[,.(library,target,expression,antibody,escape_frac,n_barcodes)]
    dt_collapse[is.na(escape_frac),n_barcodes:=NA]

Histogram below shows that a few medians end up above the theoretical
max escape fraction of 1. We set these to the maximum 1.

    hist(dt_collapse$escape_frac,main="",col="gray50",xlab="homolog escape fraction",breaks=20)

<img src="homolog_escape_files/figure-gfm/hist_median_unaltered-1.png" style="display: block; margin: auto;" />

    dt_collapse[escape_frac>1,escape_frac:=1]

Also make histograms showing the typical number of barcodes on which a
homolog escape fraction was averaged across. The median number of
barcodes across all homolog escape fracs is 257.

    hist(dt_collapse$n_barcodes,main="",col="gray50",xlab="number of barcodes",breaks=20)

<img src="homolog_escape_files/figure-gfm/hist_n_barcode-1.png" style="display: block; margin: auto;" />

How does escape frac for a homolog tend to correlate with its
expression? Not that this was a particularly well-devised check, but for
most homologs (especially expression around 10), doesn’t seem to be an
overall association between expression and escape. Perhaps notable that
the only measurement that really hits the “intermediate” range of escape
frac (instead of being more binary 0/1) is the more poorly expressing
variant AncClade2\_alt1\_subs\_only (that is, a variant that we do
expect to be perhaps less stable because of epistatic incompatibilites –
this variant contains all of the 48 substitutions that occur between
AncAsia and AncClade2, but does not contain the two deletions that also
occur along this branch).

If this variant were being disproportionately selected *out* by the RBD+
gating during the sort (because it is lower expressing than the average
homolog), this would have the effect of dragging down its theoretical
maximum escape fraction, because its theoretical maximum antibody-escape
bin counts cannot reach the relative frequency of its pre-sort counts,
because some fraction of the cells containing those barcodes aren’t even
making it through the RBD+ gate and then into the antibody-escape bin.

Consistent with this premise, it is even more notable that the four
antibody escape samples for this genotype that *do* have these
intermediate escape fractions (\~0.5) are the four antibodies that were
run together on the Aria 2-2 on the 11/9/2020 sort day – the five
antibodies run on 11/6, and the one antibody (S2X227) that were run on
the Aria 2-3 on 11/9 (so even prepped in parallel, but different sorter
with differently drawn sort gates). The FACS log from the 11/9 2-2 does
indeed show that the RBD+ gating was drawn more stringently than in the
other two batches. So, taken together, I do believe quite strongly that
these intermediate escape fractions are not distinguishable from escape
fractions of 1. And therefore, I wonder across all of the homologs
whether we should even display binding as a quantitative 0 to 1 scale,
versus just binary yes/no? Even though from the FACS log, there *were*
instances where we could see clouds of cells representing different
homologs with differing levels of escape (which we were hoping we could
then maybe see in the phenotypes some small and meaningful differences
from 0 or 1), I do not think we can effectively distinguish those
genotypes when conflated by variance in expression and inherent variance
in stringency of RBD+ gating across sort batches. Even though in this
scatterplot, we can for the maximally expressing homologs see some
variation in the escape scores away from ==0 or ==1 that might be
interpretable in the end, I don’t see a super principled way of how to
de-conflate varying expression effects here.

    plot(dt_collapse$expression,dt_collapse$escape_frac,pch=16,col="#00000066",xlab="homolog expression",ylab="homolog escape fraction")

<img src="homolog_escape_files/figure-gfm/scatter_escape_v_expression-1.png" style="display: block; margin: auto;" />

Heatmaps
--------

Last, make heatmaps illustrating the fraction escape of each homolog
versus each antibody. These will be what we’ll eventually align to the
phylogeny as our data display?

First, for all homologs in the library, including ancestors. This is the
default ordering of homologs set by the factor variable.

    #set antibody order as factor from config
    dt_collapse$antibody <- factor(dt_collapse$antibody,levels=config$antibodies_ordered)

    ggplot(dt_collapse,aes(target,antibody))+geom_tile(aes(fill=escape_frac),color="black",lwd=0.1)+
      scale_fill_gradient(low="white",high="#353D41",limits=c(0,1),na.value="cyan")+
      labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold"))

<img src="homolog_escape_files/figure-gfm/heatmap_homolog-escape_all-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$escape_scores_dir,"/heatmap_homologs_all.pdf",sep="")))

Extant homologs.

    extant <- c(config$EurAf_extant,config$SARS2_extant,config$SARS1_extant,config$Clade2_extant)

    temp <- dt_collapse[target %in% extant,];temp$target <- factor(temp$target,levels=extant)

    ggplot(temp,aes(target,antibody))+geom_tile(aes(fill=escape_frac),color="black",lwd=0.1)+
      scale_fill_gradient(low="white",high="#353D41",limits=c(0,1),na.value="cyan")+
      labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold"))

<img src="homolog_escape_files/figure-gfm/heatmap_homolog-escape_extant-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$escape_scores_dir,"/heatmap_homologs_extant.pdf",sep="")))

MAP ancestors.

    ancestors <- c(config$ancestors_MAP)

    temp <- dt_collapse[target %in% ancestors,];temp$target <- factor(temp$target,levels=ancestors)

    ggplot(temp,aes(target,antibody))+geom_tile(aes(fill=escape_frac),color="black",lwd=0.1)+
      scale_fill_gradient(low="white",high="#353D41",limits=c(0,1),na.value="cyan")+
      labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold"))

<img src="homolog_escape_files/figure-gfm/heatmap_homolog-escape_MAP-anc-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$escape_scores_dir,"/heatmap_homologs_MAP-ancestors.pdf",sep="")))

MAP plus alternative ancestors.

    ancestors <- c(config$ancestors_MAP_v_alt)

    temp <- dt_collapse[target %in% ancestors,];temp$target <- factor(temp$target,levels=ancestors)

    ggplot(temp,aes(target,antibody))+geom_tile(aes(fill=escape_frac),color="black",lwd=0.1)+
      scale_fill_gradient(low="white",high="#353D41",limits=c(0,1),na.value="cyan")+
      labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold"))

<img src="homolog_escape_files/figure-gfm/heatmap_homolog-escape_all-anc-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$escape_scores_dir,"/heatmap_homologs_all-ancestors.pdf",sep="")))

Output
------

Save our final per-homolog escape fraction estimates for each antibody.

    dt_collapse %>%
      mutate_if(is.numeric, round, digits=5) %>%
      write.csv(file=config$escape_fracs_homologs, row.names=F)
