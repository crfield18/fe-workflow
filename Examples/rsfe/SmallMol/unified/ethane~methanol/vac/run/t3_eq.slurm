#!/bin/bash
#SBATCH --job-name="eq_ethane~methanol.slurm"
#SBATCH --output="eq_ethane~methanol.slurm.slurmout"
#SBATCH --error="eq_ethane~methanol.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=21
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00

top=${PWD}
lams=(0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000)
steps=(mini heat)
if [ ! -d t3 ];then mkdir t3; fi

for step in ${steps[@]}; do
        if [ "$step" == "mini" ]; then
                for lam in ${lams[@]};do
                        mini=${lam}_mini
                        init=${lam}_init
			# check if AMBERHOME is set
			if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                        export LAUNCH="srun"
                        export EXE=${AMBERHOME}/bin/pmemd.cuda

                        ${LAUNCH} ${EXE} -O -p ${top}/unisc.parm7 -c t3/${init}.rst7 -i inputs/${mini}.mdin -o t3/${mini}.mdout -r t3/${mini}.rst7 -x t3/${mini}.nc -inf t3/${lam}.info -ref t3/${init}.rst7 
                        cat <<EOF > center.in
parm ${top}/unisc.parm7
trajin t3/${mini}.rst7
autoimage
trajout t3/${mini}_centered.rst7
go
quit
EOF
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv t3/${mini}_centered.rst7 t3/${mini}.rst7
                done
        else
		# check if AMBERHOME is set
		if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                export LAUNCH="mpirun -np ${#lams[@]}"
                export EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
                export MV2_ENABLE_AFFINITY=0
                ${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/t3_${step}.groupfile

                for lam in ${lams[@]};do
                        cat <<EOF > center.in
parm ${top}/unisc.parm7
trajin t3/${lam}_${step}.rst7
autoimage
trajout t3/${lam}_${step}_centered.rst7
go
quit
EOF
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv t3/${lam}_${step}_centered.rst7 t3/${lam}_${step}.rst7
                done
        fi

done

