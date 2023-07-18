#!/bin/bash
#SBATCH --job-name=sweep_17L_1.5ϕx_1.5ϕy_-1V0_0V1_1.2J
#SBATCH --partition=owners,simes
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=nticea@stanford.edu
#SBATCH --output=/home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/stiffness/out/sweep_17L_1.5ϕx_1.5ϕy_-1V0_0V1_1.2J_output.txt
#SBATCH --error=/home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/stiffness/out/sweep_17L_1.5ϕx_1.5ϕy_-1V0_0V1_1.2J_error.txt
#SBATCH --open-mode=append

# load Julia module
ml julia

# multithreading
export JULIA_NUM_THREADS=$SLURM_CPUS_ON_NODE

# run the script
julia /home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/stiffness/collect_data.jl 17 1 0.6180339887498949 1.0e-8 0.4487989505128276 1.5 1.5 -1 0 1.2 1