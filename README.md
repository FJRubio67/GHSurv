# GHSurv R package
The `GHSurv` R package implements the parametric General Hazards (GH) regression model studied in:

> Rubio, F. J., Remontet, L., Jewell, N. P., & Belot, A. (2019). On a general structure for hazard-based regression models: 
an application to population-based cancer research. Statistical Methods in Medical Research, 28(8), 2404-2417.

The GH model contains the Proportional hazards, accelerated hazards, and accelerated failure time models as particular cases. 

This implementation allows for using the GH model in the Overall (for hazard regression models) and Relative (for excess hazard regression models) 
survival frameworks, using the commands `GHMLE` and `GEHMLE`. 
These models are fitted using the R commands `nlminb` and `optim`. Thus, the user needs to specify the initial points and check the convergence of the
optimisation step, as usual.
For more information on these frameworks, we refer the reader to the paper
[On a general structure for hazard-based regression models: an application to population-based cancer research](https://doi.org/10.1177%2F0962280218782293).

To install the `GHSurv` R package use:

```
library(devtools)
install_github("FJRubio67/GHSurv")
```


Two illustrative examples (for the Overall and Relative survival frameworks) on the use of the `GHSurv` R package using the [Simulacrum](https://rpubs.com/FJRubio/GHSimulacrum) data set can be found at:

1. [Simulacrum data: lung cancer and Overall Survival GH model](https://rpubs.com/FJRubio/GHSimulacrum)

2. [Simulacrum data: lung cancer and Net Survival GH model](https://rpubs.com/FJRubio/GEHSimulacrum)
