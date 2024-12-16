#!/bin/bash
#SBATCH --job-name=setup_fe
#SBATCH --output=setup_fe.out        # STDOUT output file
#SBATCH --error=setup_fe.err         # STDERR output file
#SBATCH --partition=flex 
#SBATCH --nodes=1
#SBATCH -n 4
#SBATCH --export=ALL
#SBATCH --time=0-02:00:00


source ~/.bashrc
source /work2/10162/piskuliche/frontera/Software/FE-Workflow/FE-Workflow.bashrc

setup_fe
