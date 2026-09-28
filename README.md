# OED_BIODICAPT_SIMULATION

## Overview
This repository explores **Optimal Experimental Designs (OEDs)** for the BIODICAPT project. It focuses on the optimisation of adaptive sampling designs to improve the predictions of SDM models and to study the effects of actionable variables.

> [!TIP]
> For the report containing full context, method and analysis of results, see [**this page**](https://ll-mnhn.github.io/OED_biodicapt_simulation/report/report.html).

## Context
**BIODICAPT** is a French research project that aims at monitoring biodiversity of agricultural lands on a large (national) scale through the use of various recording devices.

The initiative collects data (field surveys, drone images and acoustic recordings) in agricultural plots. It is deployed in ~125 research plots, which are the initial sampled units. Through a partnership with 500 ENI, we have access to ~500 more agricultural plots. However due to logistical and budget constraints, not all of these 500 plots can be added to our sampling design.

*=> The goal of this repository is to explore the best strategies to add new samples from 500 ENI to our initial set of sampled units.*


## Project description
### Structure
```
├── report/*                # A report summarising results
├── scripts/                # Core scripts to run
├── R/                      # Functions
├── data/
│   ├── config              # Configuration files
│   ├── raw_data            # Original datasets (hidden)
│   └── preprocessed_data   # Preprocessed datasets
├── outputs/
│   ├── results             # Data outputs of each HMSC run
│   └── other_subfolders    # Visualisations and plots
├── DESCRIPTION             # Standard DESCRIPTION file for R packages
├── renv.lock               # r environment parameters
├── .Rprofile               # Default file to activate R env 
└── README.md               # This file
```

## Getting Started
### Report
The full report is the best place to start (available [here](https://ll-mnhn.github.io/pseudo-STOC-OED/report/report.html)). It is generated automatically from our results using quarto. Comments and analyses were conducted manually.

### How to use scripts
1\. Clone this repository on your machine
```bash
cd /your/local/folder
git clone https://github.com/LL-mnhn/OED_biodicapt_simulation.git
```

> [!IMPORTANT]
> Make sure you have [git-lfs](https://git-lfs.com/) installed before cloning the repository. Some large files might need to get downloaded manually if it is not installed before.

2\. Install dependencies

Open the `OED_biodicapt_simulation` folder as a new session in [Rstudio](https://docs.posit.co/ide/user/) or [Positron](https://positron.posit.co/welcome.html). Use R 4.6.1 (version used during development) 

Install `renv` if not already installed on your machine. Then run:
```R
install.packages("renv")
renv::restore()
```

> [!NOTE]
> The `rnaturalearthhires` package might sometimes not install properly through `renv`. If you get an error try to install it manually with: 
> ```R
> install.packages("pak", repos = "https://cloud.r-project.org")
> pak::pkg_install("ropensci/rnaturalearthhires")
> ```


3\. Once your environment is ready, local scripts can be run. E.g.

```R
source("scripts/3-species-simulations.R")
```

### Usage Notes
Important data and figures are already saved in the `outputs` folder.

When running on a local machine: 
- `0-verify_datasets.R` and `1-pre_processing.R` can not be run (raw dataset are only accessible by authors).
- Other files in `scripts` should run without errors.

## Contact
For inquiries, please contact Loïc Lehnhoff (UMR CESCO - MNHN) at <loic.lehnhoff@mnhn.fr> 

*This work is supervised by Nicolas Parisey (UMR IGEPP - INRAE) and Karine Princé (UMR CESCO - MNHN) as part of the BIODICAPT project (see project [page](https://www.pepr-agroeconum.fr/les-projets-finances/socio-ecosysteme/laureat-aap/biodicapt)).*
