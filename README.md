# INLAomics
INLAomics is a hierarchical Bayesian model for analysing multiomic Spatial data using Integrated Nested Laplace Approximations (INLA). Biorxiv preprint [`INLAomics for Scalable and Interpretable Spatial Multiomic Data Integration`](https://www.biorxiv.org/content/10.1101/2025.05.02.651831v1.abstract).

![github-small](https://github.com/nygctech/INLAomics/blob/main/main_1.png)


## Versions & Installation
All analysis in manuscript is carried out in `R` V. 4.3.1 with packages `R-INLA` V. 23.12.17 and `R-stan` V. 2.26.23 (`Stan` V. 2.26.1). For instructions on installation we refer to [`R-INLA`](https://www.r-inla.org/download-install) and [`mc-stan`](https://mc-stan.org/install/).

All models are implemented through the R-package [`INLA`](https://www.r-inla.org/) using the `inla.rgeneric.define()` method. The relevant scripts are under `./INLA/`

## Usage
In `./tutorial` there is a script which gives some details of the implemented INLA methods and required format of the data.

# Analysing the SPOTS data
The data generated in [1] is considered where files can be accessed at [GSE198353](www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353). We have added cell annotations to two replicates of spleen tissue sections found in `./data`.

## Spleen
The necessary files are 
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
`...1_spatial.tar.gz` and `...2_spatial.tar.gz` are our own annotations. 

Figure d above illustrates how `INLAomoics` can be utilized to construct genes-to-protein networks based on specific parameter estimates that encodes assay-to-assay effects.
Code to recreate CD3 rows are found in `./scripts/SPOTS/ProtVsGenes.R` with runtime on Apple M2 approximately 6h. Instructions for recreating any of the other rows of Figure d is outlined in the script file.

To recreate the results of figure e & f, see `./scripts/SPOTS/SpleenPred.R` 

## Breast cancer
The necessary files are 
```
.
├── GSE198353_mmtv_pymt.csv
├── GSE198353_mmtv_pymt_ADT.csv.gz
├── GSE198353_mmtv_pymt_GEX_filtered_feature_bc_matrix.h5
├── GSE198353_mmtv_pymt_spatial.tar.gz
└── spatial
    ├── aligned_fiducials.jpg
    ├── detected_tissue_image.jpg
    ├── scalefactors_json.json
    ├── tissue_hires_image.png
    ├── tissue_lowres_image.png
    └── tissue_positions_list.csv
```
Example code can be found in `./scripts/SPOTS/BreastPrediction.R`

# Other datasets
## Visium10x, tonsil
[Visium10x datasets](https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard)
```
.
├── raw_feature_bc_matrix
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── spatial
    ├── aligned_fiducials.jpg
    ├── aligned_tissue_image.jpg
    ├── cytassist_image.tiff
    ├── detected_tissue_image.jpg
    ├── scalefactors_json.json
    ├── spatial_enrichment.csv
    ├── tissue_hires_image.png
    ├── tissue_lowres_image.png
    └── tissue_positions.csv
```
Example code can be found in `./scripts/visium/tonsil.R`

## Highplex data
[GSE213264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213264)
```
.
├── GSM6578059_mousecolon_RNA.tsv.gz
├── GSM6578061_mousekidney_RNA.tsv.gz
├── GSM6578062_humantonsil_RNA.tsv.gz
├── GSM6578064_humanthymus_RNA.tsv.gz
├── GSM6578065_humanskin_RNA.tsv.gz
├── GSM6578068_mousecolon_protein.tsv.gz
├── GSM6578070_mousekidney_protein.tsv.gz
├── GSM6578071_humantonsil_protein.tsv.gz
├── GSM6578073_humanthymus_protein.tsv.gz
└── GSM6578074_humanskin_protein.tsv.gz
```
Example code can be found in `./scripts/Highplex/highplex.R`

# Simulation study
The scripts for carrying out the simulation studies are found in `./scripts/simulation/`. The script for comparison between INLA and MCMC is outlined `InlaVsMcmc.R`, where one Monte-Carlo round takes approximately 30 minutes. The script for comparison between INLAomics and univariate correlation is found in `InlaVsCorr.R`, where one Monte-Carlo round takes approximately 25 seconds.

# References
[1] Ben-Chetrit, N., Niu, X., Swett, A. D., Sotelo, J., Jiao, M. S., Stewart, C. M., ... & Landau, D. A. (2023). Integration of whole transcriptome spatial profiling with protein markers. Nature biotechnology, 41(6), 788-793.
