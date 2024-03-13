# INLAomics

Efficient methods based on Integrated Nested Laplace Approximations (INLA) for spatial generalized linear mixed models (GMMM). 

# Models
All models are implemented in [`INLA`](https://www.r-inla.org/) using the `inla.rgeneric.define()` method. The relevant scripts are under `./INLA/`

## Single assay
Gaussian Markov Random Fields (GMRF) as defined in [1] are considered

$$
\psi_j \sim \mathcal{N}_{n_j}\Big(\mathbf{0}, \tau_j^{-1}\big(\pi_j\mathbf{I} + (1-\pi_j)\mathbf{Q}_j\big)^{-1}\Big)
$$

where $\pi_j \in [0,1)$ weights between independent and spatially structured noise. 

## Multiple assays
For the non-conditional multivariate CAR (MCAR) we can follow [`INLAMSM`](https://github.com/becarioprecario/INLAMSM/tree/master) in implementation with a modified precision matrix following [1].

### RNAs


## References
[1] Leroux, B. G., Lei, X., and Breslow, N. Estimation of disease rates in small areas: a new mixed model
for spatial dependence. In _Statistical models in epidemiology, the environment, and clinical trials_, pages
179–191. Springer, 2000.
