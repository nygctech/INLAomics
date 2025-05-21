# Models

All models are implemented using the R-package [`INLA`](https://www.r-inla.org/) via the `inla.rgeneric.define()` method which allows the implementation of latent effects. This is possible because `R` can be embedded into other front-end programs in a way that allows `R` code to be run from `C`. More details of this feature can be accessed by running `inla.doc('rgeneric')` in `R`.

The definition of a new latent effect requires the definition of the latent effect as a Gaussian Markov Random Field (GMRF). This means that the mean and precision matrix of the GRMF need to be defined, as well as the hyperparameters $\theta$ involved and their prior distributions. Additionally, a binary representation of the precision matrix—used to exploit conditional independence properties (i.e., a graph)—needs to be specified.

An excellent resource on modeling with `INLA` is [Bayesian inference with INLA](https://becarioprecario.bitbucket.io/inla-gitbook/index.html), where chapter 11 provides an in-depth explanation of `rgeneric`.


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

The relevant scripts are implemented in `./INLA/MCCAR.R` and `./INLA/spotMCCAR.R`($\eta_{1,i} = 0\  \forall i$).