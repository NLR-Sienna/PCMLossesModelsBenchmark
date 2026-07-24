#!/bin/bash
#SBATCH --account=wetoowiroc
#SBATCH --output=logs/%j-%x.out
#SBATCH --time=3:59:00
#SBATCH --partition=debug
#SBATCH --licenses=gurobi@slurmdb:1


# Load required modules
module load "julia/1.12.1"
module load gurobi


# Run simulation
julia --project=. --threads=36 5bus_acopf.jl