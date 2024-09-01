# longHM

<!-- badges: start -->
<!-- badges: end -->

High-dimensional mediation analysis in longitudinal studies.


## Workflow Overview

![](man/figures/workflow.png)

## Installation

You can install the development version of longHM like so:

``` r
library(devtools)
install_github("Cui-yd/longHM")
```

## Example

This is a basic example for two different methods: linear mixed-effects model based method and generalized estimating equation based method:

``` r
library(longHM)

data(example_lmm)
results_lmm = longMediation(Y = example_lmm$Y, X = example_lmm$X, M = example_lmm$M, 
                            COV = example_lmm$COV, id = example_lmm$id, 
                            wave = example_lmm$wave, topN = c(2,2), 
                            method = "lmm", verbose = TRUE)


data(example_gee)
results_gee = longMediation(Y = example_gee$Y, X = example_gee$X, M = example_gee$M,
                            COV = example_gee$COV, id = example_gee$id,
                            wave = example_gee$wave, topN = c(2,2), 
                            method = "gee", verbose = TRUE)

```

Please do not hesitate to contact me (amber_cui@sjtu.edu.cn) if you meet any problem. Suggestions or comments are also welcome :D
