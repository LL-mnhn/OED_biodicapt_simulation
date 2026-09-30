#!/bin/bash

#SBATCH --job-name=BIODICAPT
#SBATCH --partition=std
#SBATCH --mem=96G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3 # number of chains to run (NOT number of K-folds)
#SBATCH --time=00-48:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

## Load software modules
module load R/4.6.1/gcc-openmpi
cd "$SLURM_SUBMIT_DIR"

## Stop if there is an error
set -e

## Record start time
START_TIME=$(date '+%Y-%m-%d %H:%M:%S %Z')
START_SEC=$(date '+%s')

# ## Run script (STEP 3)
# rm -rf outputs/simulations/*
# Rscript -e 'source(".Rprofile"); source("scripts/3-species_simulations.R")'

# ## Run script (STEP 4)
# rm -rf outputs/results/*
# Rscript -e 'source(".Rprofile"); source("scripts/4-OED_sampling_strategies.R")'

## Run script (STEP 5)
rm -rf outputs/results_analysis/*
Rscript -e 'source(".Rprofile"); source("scripts/5-result_analysis.R")'


## Record stop time
STOP_TIME=$(date '+%Y-%m-%d %H:%M:%S %Z')
STOP_SEC=$(date '+%s')

## Compute total runtime
RUNTIME_SEC=$((STOP_SEC - START_SEC))
RUNTIME_FMT=$(printf '%02d:%02d:%02d' $((RUNTIME_SEC/3600)) $((RUNTIME_SEC%3600/60)) $((RUNTIME_SEC%60)))

## Save datetime of run
cat > outputs/Run.md <<EOF
# Run Summary

- **Start time:** ${START_TIME}
- **Stop time:** ${STOP_TIME}
- **Total runtime:** ${RUNTIME_FMT} (hh:mm:ss)
EOF

## Finish gracefully
exit 0