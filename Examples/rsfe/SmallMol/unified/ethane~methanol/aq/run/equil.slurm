#!/bin/bash
#SBATCH --job-name="eq_ethane~methanol.slurm"
#SBATCH --output="eq_ethane~methanol.slurm.slurmout"
#SBATCH --error="eq_ethane~methanol.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --gres=gpu:8
#SBATCH --share
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=24:00:00

lams=(0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000)

top=${PWD}
# check if AMBERHOME is set
if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=${AMBERHOME}/bin/pmemd.cuda
for lam in ${lams[@]};do
    init=${lam}_init
    mini=${lam}_mini
    heat=${lam}_heat


    echo "running ${lam} mini"
    ${EXE} -O -p ${top}/unisc.parm7 -c current/${init}.rst7 -i inputs/${mini}.mdin -o current/${mini}.mdout -r current/${mini}.rst7 -x current/${mini}.nc -inf current/${mini}.info -ref current/${init}.rst7 
        cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin current/${mini}.rst7
autoimage
trajout current/${mini}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/${mini}_centered.rst7 current/${mini}.rst7
        sleep 1


    echo "running ${lam} heat"
    ${EXE} -O -p ${top}/unisc.parm7 -c current/${mini}.rst7 -i inputs/${heat}.mdin -o current/${heat}.mdout -r current/${heat}.rst7 -x current/${heat}.nc -inf current/${heat}.info 
        cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin current/${heat}.rst7
autoimage
trajout current/${heat}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/${heat}_centered.rst7 current/${heat}.rst7
        sleep 1


done
