#!/bin/bash
#SBATCH --job-name="eq_ejm42~ejm55.slurm"
#SBATCH --output="eq_ejm42~ejm55.slurm.slurmout"
#SBATCH --error="eq_ejm42~ejm55.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00

top=${PWD}
lams=(0.00000000 0.09090909 0.18181818 0.27272727 0.36363636 0.45454545 0.54545455 0.63636364 0.72727273 0.81818182 0.90909091 1.00000000)
steps=(mini heat npt)
if [ ! -d t3 ];then mkdir t3; fi

for step in ${steps[@]}; do
        if [ "$step" == "mini" ]; then
                for lam in ${lams[@]};do
                        init=${lam}_init
                        mini=${lam}_mini
			# check if AMBERHOME is set
			if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                        export LAUNCH="srun"
                        export EXE=${AMBERHOME}/bin/pmemd.cuda

                        ${LAUNCH} ${EXE} -O -p ${top}/unisc.parm7 -c t3/${init}.rst7 -i inputs/${mini}.mdin -o t3/${mini}.mdout -r t3/${mini}.rst7 -x t3/${mini}.nc -inf t3/${mini}.info -ref t3/${init}.rst7 
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
        		if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv t3/${lam}_${step}_centered.rst7 t3/${lam}_${step}.rst7
                done
        fi

done

