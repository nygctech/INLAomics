# Benchmarking

This folder contains example scripts for prediction of protein expression using [`totalVI`](https://docs.scvi-tools.org/en/1.3.2/user_guide/models/totalvi.html) and [`scArches`](https://docs.scarches.org/en/latest/#). Our implementation of these methods follow the implementations of [`MultiomeBenchmarking`](https://github.com/QuKunLab/MultiomeBenchmarking)[1].

We only include some specific examples, but the scripts are easily generalizable to any data/protein in the paper.

## Requirements
Pip freeze files can be found under `./totalvi/requirements.txt` and `./scarches/requirements.txt`. We note some difficulties in setting up the conda environments of [1] on macOS, however we ran `totalVI` on macOS with `Python v. 3.13.5` and we ran `scArches` on Linux (LinuxMint 22.1) `Python v. 3.13.5`. 

## References
[1]: Hu, Y., Wan, S., Luo, Y., Li, Y., Wu, T., Deng, W., ... & Qu, K. (2024). Benchmarking algorithms for single-cell multi-omics prediction and integration. Nature Methods, 21(11), 2182-2194.