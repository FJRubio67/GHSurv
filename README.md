# GHSurv
The `GHSurv` R package implements the parametric General Hazards (GH) regression studied in 

> Rubio, F. J., Remontet, L., Jewell, N. P., & Belot, A. (2019). On a general structure for hazard-based regression models: 
an application to population-based cancer research. Statistical methods in medical research, 28(8), 2404-2417.

This implementation allows for both using the GH model in the Overall (for hazard regression models) and Relative (for excess hazard regression models) 
survival frameworks, using the commands `GHMLE` and . For more information on these frameworks, we refer the reader to the paper
[On a general structure for hazard-based regression models: an application to population-based cancer research](https://doi.org/10.1177%2F0962280218782293).

The GH model is formulated in terms of the hazard function:
\[
h(t; {\bf x}, \alpha, \beta, \theta) = h_0(t \exp\{ {\tilde{\bf x}}^{\top}\alpha\}; \theta) \exp\{{\bf x}^{\top}\beta\},
\]
where $h_0(;\theta)$ is a parametric baseline hazard, ${\bf x}\in {\mathbb R^p}$ are the hazard-level covariates (effects), 
${\tilde{\bf x}}\in{\mathbb R}^{\tilde{p}}$ are the time-scale effects, and $\theta\in\theta$ are the parameters of the baseline hazard $h_0$. 
This model contains the Proportional hazards ($\alpha=0$), accelerated hazards ($\beta=0$), 
and accelerated failure time ($\alpha = \beta$, $\tilde{\bf x} = {\bf x}$) models as particular cases.

Two illustrative examples on the use of the `GHSurv` R package using the [Simulacrum](https://rpubs.com/FJRubio/GHSimulacrum) data set can be found at:

1. [Simulacrum data: lung cancer and Overall Survival GH model](https://rpubs.com/FJRubio/GHSimulacrum)

2. [Simulacrum data: lung cancer and Net Survival GH model](https://rpubs.com/FJRubio/GEHSimulacrum)
