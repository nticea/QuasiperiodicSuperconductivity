include("model.jl")

function submit_job(params::ModelParams, filepath, job_prefix; nodes=1, ntasks=1, cpus_per_task=8, mem=256, partition="owners")
    mkpath("out")
    mkpath("slurmfiles")

    name = "$(params.L)L_$(round(params.ϕx, digits=3))ϕx_$(round(params.ϕy, digits=3))ϕy_$(params.V0)V0_$(params.V1)V1_$(params.J)J"

    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=48:00:00
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks=$ntasks
    #SBATCH --cpus-per-task=$cpus_per_task
    #SBATCH --mem=$(mem)G
    #SBATCH --mail-type=BEGIN,FAIL,END
    #SBATCH --mail-user=nticea@stanford.edu
    #SBATCH --output=out/$(job_prefix*"_"*name)_output.txt
    #SBATCH --error=out/$(job_prefix*"_"*name)_error.txt
    #SBATCH --open-mode=append

    # load Julia module
    ml julia

    # multithreading
    export JULIA_NUM_THREADS=\$SLURM_CPUS_ON_NODE

    # run the script
    julia $filepath $(params.L) $(params.t) $(params.Q) $(params.μ) $(params.θ) $(params.ϕx) $(params.ϕy) $(params.V0) $(params.V1) $(params.J) $(params.periodic)"""

    open("slurmfiles/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch slurmfiles/$(name).slurm`)
end