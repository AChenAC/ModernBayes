# ModernBayes

The goal of ModernBayes is to implement the state-of-the-art Bayesian inference algorithm - Stein variational gradient descent (SVGD) and compare it with the gold standard Bayesian inference - Markov chain Monte Carlo (MCMC). 

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
