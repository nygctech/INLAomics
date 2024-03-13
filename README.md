# INLAomics

Spatial generalized linear mixed models (GMMM) methods for multiomic analysis using Integrated Nested Laplace Approximations (INLA). 

# Models
All models are implemented using the R-package [`INLA`](https://www.r-inla.org/) using the `inla.rgeneric.define()` method. The relevant scripts are under `./INLA/`

## Single assay
Gaussian Markov Random Fields (GMRF) as defined in [1] are considered
$$
\psi \sim \mathcal{N}_{n}\Big(\mathbf{0}, \tau^{-1}\big(\pi\mathbf{I} + (1-\pi)(\mathbf{D}-\mathbf{W})\big)^{-1}\Big)
$$
where $\pi \in [0,1)$ weights between independent and spatially structured noise. 

## Multiple assays
For the non-conditional multivariate CAR (MCAR) we can utilize parts of [`INLAMSM`](https://github.com/becarioprecario/INLAMSM/tree/master) in implementation with a modified precision matrix following [1]. Thus for the joint modeling of RNAs we refer to Sections 2.3 and 2.4 of [2]. The relevant scripts `./INLA/MCAR.R` and `./INLA/indepMCAR.R`.

The conditional GMRF, i.e., protein | RNA, is based on the approach of [3]. For the case with a single RNA in the conditioning set the relevant scripts are implemented in `./INLA/CCAR.R` and `./INLA/spotCCAR.R` where the latter restricts the cross assay effect to be from spot to spot. For the case with $G$ genes in the conditing set, the suggested extension of [3] is
$$
p(\psi^{(1)}, \psi^{(2)}_1, \ldots, \psi^{(2)}_G) 
$$

# References
[1] Leroux, B. G., Lei, X., and Breslow, N. "Estimation of disease rates in small areas: a new mixed model
for spatial dependence". In _Statistical models in epidemiology, the environment, and clinical trials_, pages
179–191. Springer, 2000.

[2] Francisco, F., Gómez-Rubio, V., and Martinez-Beneito,  M. A. "Bayesian multivariate spatial models for lattice data with INLA." arXiv preprint arXiv:1909.10804 (2019).

[3] Xiaoping, J., Carlin, B. P., and Banerjee, S. "Generalized hierarchical multivariate CAR models for areal data." Biometrics 61.4 (2005): 950-961.