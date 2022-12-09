# ModernBayes

The goal of ModernBayes is to implement a few state-of-the-art Bayesian inference algorithms including Stein variational gradient descent (SVGD) and projected Stein variational Newton (pSVN).

## Installation 
Run the following code to install the package:
```{r}
install.packages('devtools')
devtools::install_github('AChenAC/ModernBayes', build_vignettes = T)
library("ModernBayes")
```
## Examples
Run the following code to get some examples
```{r}
browseVignettes("ModernBayes")
```
