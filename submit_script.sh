#!/bin/bash
#SBATCH --account=wetoowiroc
#SBATCH --output=logs/%j-%x.out
#SBATCH --time=8:59:00
#SBATCH --qos=high
#SBATCH --partition=standard
#SBATCH --licenses=gurobi@slurmdb:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50


# Load required modules
module load "julia/1.12.1"
module load gurobi


# Run simulation
julia --project=. --threads=50 jz_test_flowcancelling_losses_cats_v2.jl