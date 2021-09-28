# GHSurv R package
The `GHSurv` R package can be used to calculate the maximum likelihood estimates and confidence intervals for the parametric General Hazards (GH) regression model studied in:

> Rubio, F. J., Remontet, L., Jewell, N. P., & Belot, A. (2019). On a general structure for hazard-based regression models: 
an application to population-based cancer research. Statistical Methods in Medical Research, 28(8), 2404-2417.
https://doi.org/10.1177%2F0962280218782293

This implementation allows for using the GH model in the Overall (hazard regression models) and Relative (excess hazard regression models) 
survival frameworks, using the commands `GHMLE` and `GEHMLE`. 
The `GHSurv` R package implements the GH model with several parametric baseline hazards, which are specified through the option `hstr` (Overall and Relative survival): lognormal (LNGH, LNGEH), log-logistic (LLGH,LLGEH), Gamma (GGH, GGEH), [Power Generalised Weibull](https://rpubs.com/FJRubio/PGW) (PGWGH, PGWGEH) , [Exponentiated Weibull](https://rpubs.com/FJRubio/EWD) (EWGH, EWGEH), and [Generalised Gamma](https://rpubs.com/FJRubio/GG) (GGGH, GGGEH) baseline hazards.

These models are fitted using the R commands `nlminb` and `optim`. Thus, the user needs to specify the initial points and to check the convergence of the
optimisation step, as usual.

To install the `GHSurv` R package use:

```
library(devtools)
install_github("FJRubio67/GHSurv")

library(GHSurv)
?GHMLE
?GEHMLE
```


Two illustrative examples (for the Overall and Relative survival frameworks) on the use of the `GHSurv` R package using the [Simulacrum](https://rpubs.com/FJRubio/GHSimulacrum) data set (available in this repository in the file `data_lung.Rda`) can be found at:

1. [Simulacrum data: lung cancer and Overall Survival GH model](https://rpubs.com/FJRubio/GHSimulacrum)

2. [Simulacrum data: lung cancer and Net Survival GH model](https://rpubs.com/FJRubio/GEHSimulacrum)
