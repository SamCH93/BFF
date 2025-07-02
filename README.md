## A Bayes Factor Framework for Unified Parameter Estimation and Hypothesis Testing

This repository contains code and data related to the preprint

>  Pawel, S. (2024). A Bayes Factor Framework for Unified Parameter Estimation and Hypothesis Testing. <https://doi.org/10.48550/arxiv.2403.09350>
  
A BibTeX entry is given by

```BibTeX
@article{Pawel2024,
  year = {2024},
  author = {Samuel Pawel},
  title = {A {Bayes} Factor Framework for Unified Parameter Estimation and Hypothesis Testing},
  url = {https://doi.org/10.48550/arxiv.2403.09350},
  note = {Preprint}
}
```

### Reproducibility

Our results can be reproduced by installing the necessary packages

``` r
## CRAN packages
pkgs <- c("brms", "ggplot2", "remotes", "knitr", "haven", "metadat")
install.packages(pkgs)

## self-developed package for meta-analysis with BFs
install.packages(pkgs = "metabf_0.1.tar.gz", repos = NULL)

## INLA (see also https://www.r-inla.org/download-install)
install.packages("INLA",
                 repos = c(getOption("repos"), 
                           INLA = "https://inla.r-inla-download.org/R/stable"), 
                 dep = TRUE) 
```

downloading the file `data.RData` from <https://osf.io/nb56d>, saving it under
`data/data.RData`, and then rerunning the code in `paper/bff.R`. To recompile
the manuscript make sure to have LaTeX installed (tested only with TeX Live
2023/Debian) and then run

``` sh
cd paper
make
```

which should produce `paper/bff.pdf`. The R and R package versions that were
used when the paper was successfully compiled before submission are visible in
the following output

``` r
sessionInfo()

#> R version 4.4.1 (2024-06-14)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=de_CH.UTF-8       
#>  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
#> [10] LC_TELEPHONE=C             LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Zurich
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] haven_2.5.5   ggplot2_3.5.2 brms_2.22.0   Rcpp_1.0.14   metabf_0.1    metadat_1.4-0
#> [7] knitr_1.50    INLA_25.06.07 Matrix_1.7-2 
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6         tensorA_0.36.2.1     xfun_0.52            #> lattice_0.22-6      
#>  [5] mathjaxr_1.6-0       vctrs_0.6.5          tools_4.4.1          #> generics_0.1.3      
#>  [9] parallel_4.4.1       sandwich_3.1-0       tibble_3.2.1         #> proxy_0.4-27        
#> [13] pkgconfig_2.0.3      KernSmooth_2.23-24   checkmate_2.3.1      #> distributional_0.4.0
#> [17] RcppParallel_5.1.10  cubature_2.1.0       lifecycle_1.0.4      #> compiler_4.4.1      
#> [21] stringr_1.5.1        fmesher_0.1.5        Brobdingnag_1.2-9    #> munsell_0.5.1       
#> [25] codetools_0.2-19     class_7.3-22         bayesplot_1.11.1     #> pillar_1.10.1       
#> [29] MASS_7.3-64          classInt_0.4-10      bridgesampling_1.1-2 #> abind_1.4-5         
#> [33] multcomp_1.4-25      nlme_3.1-167         posterior_1.6.1      #> tidyselect_1.2.1    
#> [37] mvtnorm_1.2-5        stringi_1.8.4        sf_1.0-19            #> dplyr_1.1.4         
#> [41] forcats_1.0.0        splines_4.4.1        grid_4.4.1           #> colorspace_2.1-1    
#> [45] cli_3.6.4            magrittr_2.0.3       loo_2.8.0            #> survival_3.7-0      
#> [49] TH.data_1.1-2        e1071_1.7-14         withr_3.0.2          #> scales_1.3.0        
#> [53] backports_1.5.0      sp_2.1-3             estimability_1.5     #> matrixStats_1.3.0   
#> [57] emmeans_1.10.1       hms_1.1.3            zoo_1.8-12           #> coda_0.19-4.1       
#> [61] evaluate_0.24.0      rstantools_2.4.0     rlang_1.1.5          #> xtable_1.8-4        
#> [65] glue_1.8.0           DBI_1.2.3            R6_2.6.1             units_0.8-5   
#> 
cat(paste(Sys.time(), Sys.timezone(), "\n"))

#> 2025-07-02 11:41:20.327366 Europe/Zurich 
```
