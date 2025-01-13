## A Bayes Factor Framework for Unified Parameter Estimation and Hypothesis Testing

This repository contains code and data related to the preprint

  Pawel, S. (2024). A Bayes Factor Framework for Unified Parameter Estimation
  and Hypothesis Testing. <https://doi.org/10.48550/arxiv.2403.09350>
  
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
pkgs <- c("brms", "ggplot2", "remotes", "knitr", "haven")
install.packages(pkgs)

## other packages
install.packages(pkgs = "metabf_0.1.tar.gz", repos = NULL)
remotes::install_github("wviechtb/metadat") # Bartos data only in dev version
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
#> Running under: Ubuntu 24.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Zurich
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] haven_2.5.4   ggplot2_3.5.1 INLA_24.06.27 sp_2.1-3      Matrix_1.7-1 
#>  [6] brms_2.21.0   Rcpp_1.0.13-1 metabf_0.1    metadat_1.3-0 knitr_1.48   
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6         tensorA_0.36.2.1     xfun_0.49           
#>  [4] QuickJSR_1.1.3       inline_0.3.19        lattice_0.22-5      
#>  [7] mathjaxr_1.6-0       vctrs_0.6.5          tools_4.4.1         
#> [10] generics_0.1.3       stats4_4.4.1         parallel_4.4.1      
#> [13] sandwich_3.1-0       proxy_0.4-27         tibble_3.2.1        
#> [16] pkgconfig_2.0.3      KernSmooth_2.23-24   checkmate_2.3.1     
#> [19] distributional_0.4.0 RcppParallel_5.1.9   cubature_2.1.0      
#> [22] lifecycle_1.0.4      compiler_4.4.1       stringr_1.5.1       
#> [25] fmesher_0.1.5        Brobdingnag_1.2-9    munsell_0.5.1       
#> [28] codetools_0.2-19     class_7.3-22         bayesplot_1.11.1    
#> [31] pillar_1.10.1        MASS_7.3-60.2        classInt_0.4-10     
#> [34] StanHeaders_2.32.7   bridgesampling_1.1-2 abind_1.4-5         
#> [37] multcomp_1.4-25      nlme_3.1-165         posterior_1.5.0     
#> [40] rstan_2.32.6         tidyselect_1.2.1     mvtnorm_1.2-5       
#> [43] stringi_1.8.4        sf_1.0-19            dplyr_1.1.4         
#> [46] forcats_1.0.0        splines_4.4.1        grid_4.4.1          
#> [49] colorspace_2.1-1     cli_3.6.3            magrittr_2.0.3      
#> [52] loo_2.7.0            survival_3.7-0       pkgbuild_1.4.4      
#> [55] e1071_1.7-14         TH.data_1.1-2        withr_3.0.2         
#> [58] scales_1.3.0         backports_1.5.0      estimability_1.5    
#> [61] matrixStats_1.3.0    emmeans_1.10.1       gridExtra_2.3       
#> [64] hms_1.1.3            zoo_1.8-12           coda_0.19-4.1       
#> [67] rstantools_2.4.0     rlang_1.1.4          xtable_1.8-4        
#> [70] glue_1.8.0           DBI_1.2.3            jsonlite_1.8.8      
#> [73] R6_2.5.1             units_0.8-5  

cat(paste(Sys.time(), Sys.timezone(), "\n"))

#> 2025-01-13 14:36:47.005593 Europe/Zurich
```
