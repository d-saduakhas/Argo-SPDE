# **Incorporating Correlated Nugget Effects in Multivariate SPDE Models: An Application to Argo Ocean Data**

This repository contains the reproducible code for ["**Incorporating Correlated Nugget Effects in Multivariate SPDE Models: An Application to Argo Ocean Data**"](arxiv.org)

## **Requirements**

1.  R 4.3.0 and above. The latest version of [ngme2 package](https://davidbolin.github.io/ngme2/)

2.  We assume you submit the jobs to HPC with SLURM scheduler and provide the `script_generator` bash files for automatic writing and submission of jobs scripts. Though running the individual models does not require HPC.

## Pipeline

**Model fitting**

1.  To run Gaussian models use `Code/fit_argo_gauss.R` and for non-Gaussian models use `Code/fit_argo_nig.R`
    - Assuming your run on HPC with SLURM scheduler, execute `script_generator_fit_argo.sh` to generate and automatically submit the jobs for all models for one GridID.
2.  To collect the results use `Code/collect_results.R`

**Uncertainty Quantification**

1.  To generate the QQ plots use `Code/qq_plots.R`
2.  The collected results are available in the corresponding folders for each pressure level.

**Results**

1.  The results are available in the `Results` folder. They include the `main_results.csv` and the `cv_results.csv` files. The individual model fits are not included in the repository due to the total size of the files, they are available upon request.

**Misc**

Simulation study

- Use `simulation_study_gauss.R` and `simulation_study_nig.R` to replicate the simulation study shown in paper. If executed in parallel for different settings on the HPC, use `script_generator_simulation_gauss.sh` to automatically write and submit multiple jobs at once for all the simulation settings.

List of plot functions (run in MATLAB)

- `supplement_plots_generate.m` - the automatic script to generate alll the required plots for the supplement material. Has the following individual functions:

  - `plot_correlation_common_cb.m` - plots the Pearson correlation

  - `plot_pair_nig_final.m` - plots the $\eta, \mu$ values for the different models

  - `plot_parameter.m` - plots the $\rho,\sigma,\kappa$ values for the all models

  - `plot_rho_epsilon.m` - plots the $\rho_\varepsilon$

**Data Snapshot used**

Title : Profile directory file of the Argo Global Data Assembly Center

Description : The directory file describes all individual profile files of the argo GDAC ftp site.

Project : ARGO

Format version : 2.0

Date of update : 20230827072418

FTP root number 1 : <ftp://ftp.ifremer.fr/ifremer/argo/dac>

FTP root number 2 : <ftp://usgodae.org/pub/outgoing/argo/dac>
GDAC node : CORIOLIS
