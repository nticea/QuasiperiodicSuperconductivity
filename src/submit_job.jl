include("model.jl")

function submit_job(params::ModelParams, filepath, dirpath, job_prefix; nodes=1, ntasks=1, cpus_per_task=8, mem=256, partition="owners,simes", kwargs="", time="24:00:00")
    outpath = joinpath(dirpath, "out")
    slurmpath = joinpath(dirpath, "slurmfiles")
    mkpath(outpath)
    mkpath(slurmpath)

    name = "$(params.ndims)D_$(params.L)L_$(round(params.ϕx, digits=3))ϕx_$(round(params.ϕy, digits=3))ϕy_$(round(params.ϕz, digits=3))ϕz_$(params.V0)V0_$(params.V1)V1_$(params.J)J"

    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks=$ntasks
    #SBATCH --cpus-per-task=$cpus_per_task
    #SBATCH --mem=$(mem)G
    #SBATCH --mail-type=BEGIN,FAIL,END
    #SBATCH --mail-user=nticea@stanford.edu
    #SBATCH --output=$outpath/$(job_prefix*"_"*name)_output.txt
    #SBATCH --error=$outpath/$(job_prefix*"_"*name)_error.txt
    #SBATCH --open-mode=append
    #SBATCH --sockets-per-node=2

    # load Julia module
    ml julia

    # multithreading
    export JULIA_NUM_THREADS=\$SLURM_CPUS_ON_NODE

    # run the script
    julia $filepath $(params.L) $(params.t) $(params.Q) $(params.μ) $(params.θ) $(params.ϕx) $(params.ϕy) $(params.ϕz) $(params.V0) $(params.V1) $(params.J) $(params.periodic) $(params.ndims) $(params.disorder) $kwargs"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end