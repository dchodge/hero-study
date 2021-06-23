# Analysis of the HERO Study

## Summary of project

This repository contains an R package which performs the data cleaning, model fitting and plotting for the analysis of the SARS-CoV-2 HERO study performed at Sheffield Teaching Hospital. The mcmc sampling is done in cmdstanr (not rstan).

## Installation guide

*Quick-install instructions*

Clone the repository and use `devtools::load_all()` to load everything into R. Then go through the vignettes in the `/vignettes` folder.

*Step-by-step installation instructions*

* Before cloning, make sure you have installed [github lfs](https://git-lfs.github.com/) on your computer.
* In R, make sure you have installed and loaded [devtools](https://devtools.r-lib.org/).
* Clone this repository onto your computer (e.g. git clone https://github.com/dchodge/hero-study.git)
* Load R and make sure your working directory is `/YOUR-PATH/hero-study`
* Use `devtools::load_all()` to install the package. This package has a number of dependency packages (tidybayes, cmdstanr, posterior, rstan, rstanarm etc.) which must be install in order for the installation to work. `devtools::load_all()` will throw an error if there are any package dependencies missing on your machine. Just use `install.packages(“MISSING_PACKAGE_NAME”)` to install them. 
* Now consult the vignettes in the `/vignettes` folder to run the code.

## Overview of files

There are two Markdown files in the `/vignettes` folder (study_run_figs.Rmd and study_other):
* `/vignettes/study_run_figs.Rmd` — runs the three models in the study and plots the figures
* `/vignettes/study_other.Rmd `— produces everything else needed for the study (tables, point estimations etc)

The guts of the code is in the `/R` folder. The scripts are as follows:
* `/R/data_clean.R` — cleans the raw serological data into a R-friendly version
* `/R/data_model.R` — creates datasets for each of the models
* `/R/m1_prev.R` — All the code needed to run, process the posterior distributions, and plot the figures for the prevalence model
* `/R/m2_abkin.R` — All the code needed to run, process the posterior distributions, and plot the figures for the antibody kinetics model
* `/R/m3_hetsen.R` — All the code needed to run, process the posterior distributions, and plot the figures for the hetergeneous sensitivity and specificity model
* R/supp.R` — code to make the supplementary figs

The mcmc samples are in the `/data` folder:
*`/data/datafit.RData` Is all the data required to fit the mcmc models. This is saved with the `save_all_datafit()` function, and loaded into the environment with the `load_all_datafit()` function.
*`/data/fit_mcmc_prev.RData` are the raw mcmc samples from the cmdstanr model for the prevalence model shown in the Figure 1.
 *`/data/fit_mcmc_asymp.RData` are the raw mcmc samples from the cmdstanr model for the asymptomatic model shown in the Figure SX. 
*`/data/fit_mcmc_start.RData` are the raw mcmc samples from the cmdstanr model for the starting titre value model shown in the Figure 2.
 *`/data/fit_mcmc_change.RData` are the raw mcmc samples from the cmdstanr model for the change in titre value model shown in the Figure 2.
*`/data/fit_mcmc_sens.RData` are the raw mcmc samples from the cmdstanr model for the heterogenous sensitivity model shown in Figure 3.
 *`/data/fit_mcmc_spec.RData` are the raw mcmc samples from the cmdstanr model for the heterogenous specificity model shown in Figure SX.

The stan code and the raw mcmc outputs is contained in the `outputs/` folder.

The figures used in the final article are in the `outputs/figs/` folder.

## Linked publications

This research is still work-in-progress, so the manuscript is not available yet.

## Contact details

 Please email david.hodgson@lshtm.ac.uk with any queries relating to this code.