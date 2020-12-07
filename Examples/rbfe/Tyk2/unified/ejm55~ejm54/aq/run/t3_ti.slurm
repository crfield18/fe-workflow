#!/bin/bash
#SBATCH --job-name="pr_ejm55~ejm54.slurm"
#SBATCH --output="pr_ejm55~ejm54.slurm.slurmout"
#SBATCH --error="pr_ejm55~ejm54.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00

lams=(0.00000000 0.09090909 0.18181818 0.27272727 0.36363636 0.45454545 0.54545455 0.63636364 0.72727273 0.81818182 0.90909091 1.00000000)
# check if AMBERHOME is set
if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
echo "running replica ti"
mpirun -np ${#lams[@]} ${EXE} -rem 3 -ng ${#lams[@]} -groupfile inputs/t3_ti.groupfile
