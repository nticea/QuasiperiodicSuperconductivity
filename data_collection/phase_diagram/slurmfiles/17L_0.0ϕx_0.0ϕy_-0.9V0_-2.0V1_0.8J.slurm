#!/bin/bash
#SBATCH --job-name=phase_diagram_17L_0.0ϕx_0.0ϕy_-0.9V0_-2.0V1_0.8J
#SBATCH --partition=owners
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=nticea@stanford.edu
#SBATCH --output=/home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/phase_diagram/out/phase_diagram_17L_0.0ϕx_0.0ϕy_-0.9V0_-2.0V1_0.8J_output.txt
#SBATCH --error=/home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/phase_diagram/out/phase_diagram_17L_0.0ϕx_0.0ϕy_-0.9V0_-2.0V1_0.8J_error.txt
#SBATCH --open-mode=append

# load Julia module
ml julia

# multithreading
export JULIA_NUM_THREADS=$SLURM_CPUS_ON_NODE

# run the script
julia /home/users/nticea/code/QuasiperiodicSuperconductivity/data_collection/phase_diagram/collect_data.jl 17 1 0.6180339887498949 1.0e-8 0.4487989505128276 0 0 -0.9 -2.0 0.8 1