# MCAR

Efficient methods based on Integrated Nested Laplace Approximations (INLA) for spatial generalized linear mixed models (GMMM).

## Models
### Single assay
$$
p\big(y_{ij}|\lambda_{ij}, \phi_j\big) = \begin{cases}
\phi_j + (1-\phi_j)\exp\{-\exp\{s_{ij}\lambda_{ij}\}\}, & y_{ij} = 0 \\
(1-\phi_j)\exp\left\lbrace\big(y_{ij}\log(s_{ij}\lambda_{ij}) - \exp\{s_{ij}\lambda_{ij}\}\big)\right\rbrace, & y_{ij} > 1, 
\end{cases},\ i = 1,\ldots,n_j, \ j = 1,\ldots, J,
$$

where the Poisson rates are modelled as

$$
\log \lambda_j = \mathbf{X}_{j}^T\beta + \psi_j,\ j = 1, \ldots, J,
$$

where $\boldsymbol{\psi}_j$ are Gaussian Markov Random Fields (GMRF) defined as in [1]

$$
\psi_j \sim \mathcal{N}_{n_j}\Big(\mathbf{0}, \tau_j^{-1}\big(\pi_j\mathbf{I} + (1-\pi_j)\mathbf{Q}_j\big)^{-1}\Big)
$$

where $\pi_j \in [0,1)$ weights between independent and spatially structured noise. 

### Multiple assays

## References
[1] Leroux, B. G., Lei, X., and Breslow, N. Estimation of disease rates in small areas: a new mixed model
for spatial dependence. In _Statistical models in epidemiology, the environment, and clinical trials_, pages
179–191. Springer, 2000.
