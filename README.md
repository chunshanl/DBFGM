# Dynamic Bayesian Functional Graphical Models
_____________________________

Author: Chunshan Liu

Contact: chunshanl@hotmail.com

This is the R code to reproduce results in the following paper:
>Dynamic Bayesian Functional Graphical Models, Chunshan Liu, Daniel R. Kowal, and Marina Vannucci
______________________________

## Introduction

The repository contains the R code for the proposed Dynamic Bayesian Functional Graphical Model. The proposed model is a Bayesian graphical model for multivariate functional data where each curve is a function over time. Each node in the resulting graph is a function inside a time interval.

An example of applying the proposed model and doing posterior analysis can be found in Application_study. The repository also contains R code to reproduce all results in the simulation study and application study.

## Files:

- Simulation_study.R
  - This script can reproduce MCMC results in the simulation study. It includes a script to randomly generate precision matrices and data sets.

- Application_study.R
  - This script can reproduce results in the application study on the sea surface temperature data. 
  - This script contains an example to run the proposed dynamic Bayesian functional graphical model.

- MCMC_algorithms
  - This folder contains the MCMC algorithms of proposed dynamic Bayesian Graphical models.

- Helper functions
  - This folder contains helper functions for the MCMC algorithms and for regenerating the results presented in the paper. They are called in the script and don't need to be run independently.
  - call_DBFGM.R: This is the functions to call the MCMC algorithm of the proposed dynamic Bayesian functional graphical model in this paper. It generate initial values, and pass data and hyperparameters into the MCMC algorithms.
  - MCMC_changepoint_DBFGM.R: MCMC algorithm of the proposed dynamic Bayesian functional graphical model. It is called by call_DBFGM.R.

- sst.csv
  - The sea surface temperature data set. It is available on the [ERA5 website] (https://cds.climate.copernicus.eu/cdsapp#!/search?type=dataset).

## Acknowledgements

The code provided here includes code associated with the following publications:

> Wang H. (2015). Scaling it up: Stochastic search structure learning in graphical models. Bayesian Analysis, 10(2), 351-377.