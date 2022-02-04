#!/bin/bash
#SBATCH --job-name="pr_1h1q~1oiu.slurm"
#SBATCH --output="pr_1h1q~1oiu.slurm.slurmout"
#SBATCH --error="pr_1h1q~1oiu.slurm.slurmerr"
#SBATCH --partition=general-long-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=11
#SBATCH --gres=gpu:8
#SBATCH --time=3-00:00:00

lams=(0.00000000 0.10000000 0.20000000 0.30000000 0.40000000 0.50000000 0.60000000 0.70000000 0.80000000 0.90000000 1.00000000)
# check if AMBERHOME is set
if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
echo "running replica ti"
mpirun -np ${#lams[@]} ${EXE} -rem 3 -ng ${#lams[@]} -groupfile inputs/ti.groupfile
