# Analysis of the HERO Study

This repository contains an R package which performs the data cleaning, model fitting and plotting for the analysis of the SARS-CoV-2 HERO study performed at Sheffield Teaching Hospital. The mcmc sampling is done in rstan and cmdstanr. 

There are two Markdown files in the vignettes folder (study_run_figs.Rmd and study_other).
* study_run_figs.Rmd - runs the three models in the study and plots the figures
* study_other.Rmd - produces everything else needed for the study (tables, point estimations etc)

The guts of the code is in the R folder. The scripts are as follows:
* data_clean.R — cleans the raw serological data into a R-friendly version
* data_model.R — creates datasets for each of the models
* m1_prev.R — All the code needed to run, process the posterior distributions, and plot the figures for the prevalence model
* m2_abkin.R — All the code needed to run, process the posterior distributions, and plot the figures for the antibody kinetics model
* m3_hetsen.R — All the code needed to run, process the posterior distributions, and plot the figures for the hetergeneous sensitivity and specificity model
* supp.R — code to make the supplementary figs