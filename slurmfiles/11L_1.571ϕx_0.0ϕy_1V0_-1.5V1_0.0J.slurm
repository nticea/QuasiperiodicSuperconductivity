#!/bin/bash
#SBATCH --job-name=sweep_11L_1.571ϕx_0.0ϕy_1V0_-1.5V1_0.0J
#SBATCH --partition=owners
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=nticea@stanford.edu
#SBATCH --output=out/sweep_11L_1.571ϕx_0.0ϕy_1V0_-1.5V1_0.0J_output.txt
#SBATCH --error=out/sweep_11L_1.571ϕx_0.0ϕy_1V0_-1.5V1_0.0J_error.txt
#SBATCH --open-mode=append

# load Julia module
ml julia

# multithreading
export JULIA_NUM_THREADS=$SLURM_CPUS_ON_NODE

# run the script
julia /home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/stiffness/collect_data.jl 11 1 0.6180339887498949 1.0e-8 0.4487989505128276 1.5707963267948966 0.0 1 -1.5 0.0 1