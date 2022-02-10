#!/bin/bash

#SBATCH --partition=XXXX
#SBATCH --qos=XXXX
#SBATCH --nodes=4
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=10
#SBATCH --mem=40GB
#SBATCH --job-name=pjulia
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --account=XXXX
#SBATCH --mail-user=XXXX

ml julia

# Delete this line if you have already installed the packages
julia julia_dependencies.jl

julia abcmk_cli.jl rates --samples 661 --gamNeg -2000,-200 --gL 1,10 --gH 200,2000 --rho 0.001 --theta 0.001 --solutions 100000 --output /home/jmurga/rates.jld2 --dac 1,2,4,5,10,20,50,100,200,400,500,661,925,1000 --scheduler slurm --nthreads 40

julia abcmk_cli.jl /home/jmurga/parseData --analysisFolder /home/jmurga/tgp/

julia abcmk_cli.jl summaries --analysisFolder /home/jmurga/tgp/ --rates rates2.jld2 --samples 661 --replicas 100 --summstatSize 100000 --dac 2,4,5,10,20,50,200,661,925 --nthreads 8 --scheduler slurm

julia abcmk_cli.jl abcInference --analysisFolder /home/jmurga/tgp/ --replicas 100 --P 5 --S 9 --tol 0.001 --ABCreg /home/jmurga/ABCreg/src/reg --parallel false --nthreads 1
