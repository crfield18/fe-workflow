#!/bin/bash
#SBATCH --job-name="pr_ethane~methanol.slurm"
#SBATCH --output="pr_ethane~methanol.slurm.slurmout"
#SBATCH --error="pr_ethane~methanol.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00

lams=(0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000)
# check if AMBERHOME is set
if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
echo "running replica ti"
mpirun -np ${#lams[@]} ${EXE} -rem 3 -ng ${#lams[@]} -groupfile inputs/ti.groupfile
