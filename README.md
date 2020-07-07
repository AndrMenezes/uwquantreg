
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `uwquantreg`: unit-Weibull quantile regression

<!-- badges: start -->

<!-- badges: end -->

The goal of `uwquantreg` is to provide tools for fitting unit-Weibull
quantile regression proposed by [Mazucheli et
al. (2020)](https://www.tandfonline.com/doi/abs/10.1080/02664763.2019.1657813?journalCode=cjas20).
The Zero-or-One inflation model are also available.

## Installation

You can install the development version of `uwquantreg` from
[GitHub](https://github.com/AndrMenezes/uwquantreg) with:

``` r
# install.packages("devtools")
devtools::install_github("AndrMenezes/uwquantreg")
```

## Example

The packages follows the structure of `glm` objects. The main functions
are `uwquantreg` and `zouwquantreg`.

``` r
library(uwquantreg)
#> 
#> Attaching package: 'uwquantreg'
#> The following object is masked from 'package:datasets':
#> 
#>     trees

fit <- uwquantreg(phpws ~ mhdi + incpc + region + log(pop), 
                  data = water, 
                  tau = 0.5, 
                  link = "logit")
fit
#> 
#>  unit-Weibull quantile regression model
#> Call:  uwquantreg(formula = phpws ~ mhdi + incpc + region + log(pop), 
#>     tau = 0.5, data = water, link = "logit")
#> 
#> Mu coefficients (quantile model with logit link and tau = 0.5): 
#> (Intercept)         mhdi        incpc       region     log(pop)  
#>     -7.8176      13.6422      -0.0004      -0.3526       0.1555  
#> 
#> Model with constant shape parameter:
#>  phi  
#> 1.18
summary(fit)
#> 
#>  Wald-tests for unit-Weibull quantile regression model
#> 
#> Call:  uwquantreg(formula = phpws ~ mhdi + incpc + region + log(pop), 
#>     tau = 0.5, data = water, link = "logit")
#> 
#> Mu coefficients: (quantile model with logit link and tau = 0.5): 
#>               Estimate Std. Error Z value Pr(>|z|)    
#> (Intercept) -7.8175550  0.3136689 -24.923  < 2e-16 ***
#> mhdi        13.6421970  0.5133700  26.574  < 2e-16 ***
#> incpc       -0.0004465  0.0001037  -4.306 1.67e-05 ***
#> region      -0.3525859  0.0489174  -7.208 5.69e-13 ***
#> log(pop)     0.1555402  0.0154771  10.050  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Model with constant shape:
#>     Estimate Std. Error Z value Pr(>|z|)    
#> phi  1.18031    0.01643   71.85   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
confint(fit)
#>                     2.5 %        97.5 %
#> (Intercept) -8.4323346771 -7.2027752861
#> mhdi        12.6360102356 14.6483838346
#> incpc       -0.0006497558 -0.0002432488
#> region      -0.4484621923 -0.2567096379
#> log(pop)     0.1252055979  0.1858747188
#> phi          1.1481128532  1.2125099849
```

The currently methods implemented are

``` r
methods(class = "uwquantreg")
#>  [1] coef      confint   fitted    hnp       logLik    predict   print    
#>  [8] residuals summary   vcov     
#> see '?methods' for accessing help and source code
```
