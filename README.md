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

The conditional GMRF, i.e., protein | RNA, is based on the approach of [3]. For the case with a single RNA in the conditioning set the relevant scripts are implemented in `./INLA/CCAR.R` and `./INLA/spotCCAR.R` where the latter restricts the cross assay effect to be from spot to spot. Note that `./STAN/CCAR.stan` is the corresponding implementatin of `./INLA/CCAR.R` in stan using a Poisson likelihood. For the case with $G$ genes in the conditing set, the suggested extension of [3] is

$$
p(\psi^{(1)}, \psi^{(2)}_1, \ldots, \psi^{(2)}_G)= p(\psi^{(1)} | \psi^{(2)}_1, \ldots, \psi^{(2)}_G) p(\psi^{(2)}_1, \ldots, \psi^{(2)}_G),
$$

where the protein GMRF is 

$$
p(\psi^{(1)} | \psi^{(2)}_1, \ldots, \psi^{(2)}_G) = \mathcal{N}_n \bigg( \sum _{i=1}^G \big(\eta _{0,i}\mathbf{I} + \eta _{1,i}\mathbf{W}\big)\psi_i^{(2)}, \tau_1^{-1}\big(\pi_1\mathbf{I} + (1-\pi_1)(\mathbf{D}-\mathbf{W})\big)^{-1} \bigg).
$$

The relecant scripts are implemented in `./INLA/MCCAR.R` and `./INLA/spotMCCAR.R`.

# Analysing the SPOTS data
The data generated in [4] is considered, where we have added cell annotations to two replicates of spleen tissue sections. The files needed to recreate the analysis are

```
.
├── GSE198353_spleen_rep_1.csv
├── GSE198353_spleen_rep_1_filtered_feature_bc_matrix.h5
├── GSE198353_spleen_rep_2.csv
├── GSE198353_spleen_rep_2_filtered_feature_bc_matrix.h5
├── GSE198353_spleen_replicate_1_spatial.tar.gz
├── GSE198353_spleen_replicate_2_spatial.tar.gz
├── spatial
│   ├── qc_aligned_fiducials_image.jpg
│   ├── qc_detected_tissue_image.jpg
│   ├── scalefactors_json.json
│   ├── tissue_hires_image.png
│   ├── tissue_lowres_image.png
│   └── tissue_positions_list.csv
├── spatial2
│   ├── qc_aligned_fiducials_image.jpg
│   ├── qc_detected_tissue_image.jpg
│   ├── scalefactors_json.json
│   ├── tissue_hires_image.png
│   ├── tissue_lowres_image.png
│   └── tissue_positions_list.csv
```
`...1_spatial.tar.gz` and `...2_spatial.tar.gz` are our own annotations, the remaining files can be found using [GSE198353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353)

The Figure below outlines estimation of $\eta_0$ (left) and $\eta_1$ (right) when restricting the candidate models to a set of proteins and their paired genes. Subsequentily, from the model with the relevant pair we add on genes in the conditing set based on a variable selection procedure. The set of RNAs that are conditioned on is then expanded until there is a drop in the Deviance Information Criterion (DIC). Solid dots represent significant effects, in the sense that their $95$% credible sets does not cover $0$. Code to recreate CD3 rows are found in `./scripts/SPOTS/ProtVsGenes.R`

![github-small](https://github.com/nygctech/INLAomics/blob/main/ProtVsGenes.png)



# References
[1] Leroux, B. G., Lei, X., and Breslow, N. "Estimation of disease rates in small areas: a new mixed model
for spatial dependence". In _Statistical models in epidemiology, the environment, and clinical trials_, pages
179–191. Springer, 2000.

[2] Francisco, F., Gómez-Rubio, V., and Martinez-Beneito,  M. A. "Bayesian multivariate spatial models for lattice data with INLA." arXiv preprint arXiv:1909.10804 (2019).

[3] Xiaoping, J., Carlin, B. P., and Banerjee, S. "Generalized hierarchical multivariate CAR models for areal data." Biometrics 61.4 (2005): 950-961.

[4] Ben-Chetrit, N., Niu, X., Swett, A. D., Sotelo, J., Jiao, M. S., Stewart, C. M., ... & Landau, D. A. (2023). Integration of whole transcriptome spatial profiling with protein markers. Nature biotechnology, 41(6), 788-793.
