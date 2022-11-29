## Overview

Experimental and phylogenetic dynamics modeling pipeline. 

## How to use 
```

Whole pipeline runs with bash_scripts\expdyn_phydyn_v3.sh

Parameters:
INPUT_DIR=$1 # Absolute path input dir
OUT_DIR=$2 # Absolute path out dir
EQUATION_FILE=$3 # EQUATIONS and init min and max values for variables and parameters file absolute path
MSA_FASTA_OUT_FILE=$4 # MSA result in fasta format file absolute path
FASTA_INFO_FILE=$5 # File with information about all fasta sequences (sequence_name,time(as year.part_of_year),variable_name_in_ODE)
MCMC_ITERATIONS=$6 # Number of MCMC iterations
LOG_STEP=$7 # Size of MCMC logging step
MUTATION_RATE=$8 # Mutation rate for sequences in MCMC (default=1.1)
TARGET_FILE=$9 # Target trajectory file for Experimental dynamic
GLOB_ITERATIONS=$10 # Num of iterations of exp_dyn and phy_dyn (default=1)
FIRST_DYN=$11 # Set which dynamics will started first: Experimental dynamics (expdyn) or Phylodynamics (phydyn) (default=expdyn)
NUM_GOOD_RES_TR_L=$12 # The number of parameter sets that have a good residual and tree probability and have maximum variability
wR=$13 # Coeficient of Residual/min_residual ratio for weighted average sum (wR)

```

### Example
Example input files in /example_in directory.