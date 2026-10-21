# OED_BIODICAPT_SIMULATION

## Overview
This repository explores **Optimal Experimental Designs (OEDs)** for the BIODICAPT project. It focuses on the optimisation of adaptive sampling designs to improve the predictions of SDM models and to study the effects of actionable variables within the BIODICAPT research project.

> [!TIP]
> For the report containing full context, methods and the analyses of results, see [**this document**](https://ll-mnhn.github.io/OED_biodicapt_simulation/report/report.html).

## Context
**BIODICAPT** is a French research project that aims at monitoring biodiversity of agricultural lands on a large (i.e. national) scale through the use of various recording devices.

The initiative collects data (field surveys, drone images and acoustic recordings) in agricultural fields. It is deployed in ~125 research plots, which are the initial sampled units. Through a partnership with another project (the 500 ENI network), the BIODICAPT project has access to ~500 more agricultural plots. However due to logistical and budget constraints, not all of these 500 plots can be added to our sampling design.

Therefore, the goal of this repository is to explore the best strategies to add new samples from 500 ENI to the initial set of sampled units.


## Project description
### Structure
```
├── report/*                # A report summarising results
├── scripts/                # Core scripts to run (indexed by execution order)
├── R/                      # Functions for scripts
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
The full report is the best place to start (available [here](https://ll-mnhn.github.io/pseudo-STOC-OED/report/report.html)). It is generated automatically from the results, using quarto. 

text, comments and analyses are written manually.

### How to use
1\. Clone this repository on your machine
```bash
cd /your/local/folder
git clone https://github.com/LL-mnhn/OED_biodicapt_simulation.git
```

> [!IMPORTANT]
> Make sure you have [`git-lfs`](https://git-lfs.com/) installed **BEFORE** cloning the repository, as some files (rasters mainly) are too large for base `git`.

2\. Install dependencies

Open the `OED_biodicapt_simulation` folder as a new session in either [Rstudio](https://docs.posit.co/ide/user/) or [Positron](https://positron.posit.co/welcome.html). Use R 4.6.1 (version used during development).

Install `renv` if not already installed on your machine. Then run:
```R
install.packages("renv")
source(.Rprofile)
renv::restore()
```

> [!NOTE]
> The `rnaturalearthhires` package sometimes does not install properly through `renv`. If you get an error try to install it manually with: 
> ```R
> install.packages("pak", repos = "https://cloud.r-project.org")
> pak::pkg_install("ropensci/rnaturalearthhires")
> ```


3\. Once your environment is ready, you can run any script of this repo. E.g.

```R
source("scripts/3-species-simulations.R")
```

### Usage Notes
Important data and figures are already saved and ready in the `outputs` folder.

When running on a local machine: 
- `0-verify_datasets.R` and `1-pre_processing.R` need the raw datasets (only accessible by authors for confidentiality reasons) and might not work properly.
- Other files in `scripts` should be safe to run without errors.


## Disclaimer
Some parts of this code were written/improved using [Anthropic](https://www.anthropic.com/)'s AI tool [claude.ai](https://claude.ai/). These parts were **systematically reviewed and verified**, no code generated directly by [claude.ai](https://claude.ai/) was left 'as is' in the scripts.

**On a more personal note**: even though LLMs are powerful tools, they should only be used as starting points. Like Academic articles, Wikipedia pages or even Reddit posts: **AI's answers are just a source of information whose reliability varies**. I believe it should be common practice to disclose their use wherever they were needed, just as with any other source.

## Acknowledgements
This work was supervised by Nicolas Parisey (UMR IGEPP - INRAE) and Karine Princé (UMR CESCO - MNHN) as part of the [BIODICAPT project](https://www.pepr-agroeconum.fr/les-projets-finances/socio-ecosysteme/laureat-aap/biodicapt) funded by the ["Agroécologie et Numérique" PEPR](https://www.pepr-agroeconum.fr/).

## Contact
For inquiries, please contact Loïc Lehnhoff (UMR CESCO - MNHN) at <loic.lehnhoff@mnhn.fr>. 
