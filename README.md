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
pkgs <- c("ReplicationSuccess", "brms", "ggplot2", "remotes", "knitr")
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

and then rerunning the code in `paper/bff.R`. To recompile the manuscript make
sure to have LaTeX installed (tested only with TeX Live 2022/dev/Debian) and
then run

``` sh
cd paper
make
```

which should produce `paper/bff.pdf`. The R and R package versions that were
used when the paper was successfully compiled before submission are visible in
the following output

``` r
sessionInfo()

#> R version 4.3.3 (2024-02-29)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 22.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
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
#>  [1] ggplot2_3.5.0            INLA_24.02.09            sp_2.1-3                
#>  [4] Matrix_1.6-3             brms_2.20.4              Rcpp_1.0.12             
#>  [7] metabf_0.1               ReplicationSuccess_1.3.2 metadat_1.3-0           
#> [10] knitr_1.45              
#> 
#> loaded via a namespace (and not attached):
#>   [1] DBI_1.1.3            gridExtra_2.3        inline_0.3.19       
#>   [4] sandwich_3.1-0       rlang_1.1.1          magrittr_2.0.3      
#>   [7] multcomp_1.4-25      e1071_1.7-13         matrixStats_1.0.0   
#>  [10] compiler_4.3.3       loo_2.6.0            callr_3.7.3         
#>  [13] vctrs_0.6.5          reshape2_1.4.4       stringr_1.5.0       
#>  [16] pkgconfig_2.0.3      crayon_1.5.2         fastmap_1.1.1       
#>  [19] backports_1.4.1      ellipsis_0.3.2       utf8_1.2.3          
#>  [22] threejs_0.3.3        promises_1.2.1       markdown_1.12       
#>  [25] ps_1.7.5             MatrixModels_0.5-3   xfun_0.40           
#>  [28] jsonlite_1.8.7       later_1.3.1          parallel_4.3.3      
#>  [31] prettyunits_1.1.1    R6_2.5.1             dygraphs_1.1.1.6    
#>  [34] stringi_1.7.12       StanHeaders_2.32.5   estimability_1.4.1  
#>  [37] rstan_2.32.5         zoo_1.8-12           base64enc_0.1-3     
#>  [40] bayesplot_1.11.0     httpuv_1.6.11        splines_4.3.3       
#>  [43] igraph_2.0.1.1       tidyselect_1.2.0     abind_1.4-5         
#>  [46] codetools_0.2-19     miniUI_0.1.1.1       processx_3.8.2      
#>  [49] pkgbuild_1.4.2       lattice_0.22-5       tibble_3.2.1        
#>  [52] plyr_1.8.9           withr_2.5.0          shiny_1.7.5         
#>  [55] bridgesampling_1.1-2 posterior_1.5.0      coda_0.19-4         
#>  [58] survival_3.5-7       sf_1.0-14            units_0.8-4         
#>  [61] proxy_0.4-27         RcppParallel_5.1.7   xts_0.13.2          
#>  [64] pillar_1.9.0         tensorA_0.36.2.1     KernSmooth_2.23-22  
#>  [67] checkmate_2.2.0      DT_0.31              stats4_4.3.3        
#>  [70] shinyjs_2.1.0        distributional_0.4.0 generics_0.1.3      
#>  [73] mathjaxr_1.6-0       rstantools_2.4.0     munsell_0.5.0       
#>  [76] scales_1.3.0         gtools_3.9.5         xtable_1.8-4        
#>  [79] class_7.3-22         glue_1.6.2           cubature_2.1.0      
#>  [82] emmeans_1.8.9        tools_4.3.3          shinystan_2.6.0     
#>  [85] colourpicker_1.3.0   mvtnorm_1.2-3        grid_4.3.3          
#>  [88] QuickJSR_1.1.3       crosstalk_1.2.1      colorspace_2.1-0    
#>  [91] nlme_3.1-163         cli_3.6.1            fansi_1.0.4         
#>  [94] fmesher_0.1.5        Brobdingnag_1.2-9    dplyr_1.1.4         
#>  [97] gtable_0.3.4         digest_0.6.33        classInt_0.4-10     
#> [100] TH.data_1.1-2        farver_2.1.1         htmlwidgets_1.6.2   
#> [103] htmltools_0.5.6      lifecycle_1.0.3      mime_0.12           
#> [106] shinythemes_1.2.0    MASS_7.3-60.0.1 
 
cat(paste(Sys.time(), Sys.timezone(), "\n"))

#> 2024-03-12 15:11:38.508182 Europe/Zurich 
```
