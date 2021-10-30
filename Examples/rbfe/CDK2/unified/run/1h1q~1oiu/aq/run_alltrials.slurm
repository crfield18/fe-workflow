#!/bin/bash
#SBATCH --job-name="eq_1h1q~1oiu.slurm"
#SBATCH --output="eq_1h1q~1oiu.slurm.slurmout"
#SBATCH --error="eq_1h1q~1oiu.slurm.slurmerr"
#SBATCH --partition=general-long-gpu
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=11
#SBATCH --gres=gpu:8
#SBATCH --time=3-00:00:00

top=${PWD}
endstates=(0.00000000)
lams=(0.00000000 0.10000000 0.20000000 0.30000000 0.40000000 0.50000000 0.60000000 0.70000000 0.80000000 0.90000000 1.00000000)
twostate=true
eqstage=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqV eqP eqA minTI eqpre1P0TI eqpre2P0TI eqP0TI eqATI preTI)


# check if AMBERHOME is set
#if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

for trial in $(seq 1 1 3); do

	if [ ! -d t${trial} ];then mkdir t${trial}; fi

	count=-1; alllams=0
	for stage in ${eqstage[@]}; do
        	count=$((${count}+1))
        	lastcount=$((${count}-1))
        	if [ "${stage}" == "init" ] || [ "${stage}" == "eqpre1P0TI" ] || [ "${stage}" == "eqpre2P0TI" ] || [ "${stage}" == "eqP0TI" ]; then continue; fi
        	laststage=${eqstage[${lastcount}]}

        	if [ "${stage}" == "minTI" ];then alllams=1; fi

        	if [ ${alllams} -eq 0 ];then

                	# check if pmemd.cuda is present
                	if ! command -v ${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi

                	export LAUNCH="srun"
                	export EXE=${AMBERHOME}/bin/pmemd.cuda

                	lam=${endstates[0]}
                	echo "Running $stage for lambda ${lam}..."
                	${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
                	cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin t${trial}/${lam}_${stage}.rst7
autoimage
trajout t${trial}/${lam}_${stage}_centered.rst7
go
quit
EOF2
                	# check if cpptraj is present
                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                	cpptraj < center.in
                	sleep 1
                	mv t${trial}/${lam}_${stage}_centered.rst7 t${trial}/${lam}_${stage}.rst7

        	elif [ ${alllams} -eq 1 ] && [ "${stage}" == "minTI" ];then
                	# check if pmemd.cuda is present
                	if ! command -v ${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
                	export LAUNCH="srun"
                	export EXE=${AMBERHOME}/bin/pmemd.cuda
			for i in ${!lams[@]}; do
				lam=${lams[$i]}
				if [ "${i}" -eq 0 ]; then
					init=${endstates[0]}_eqA.rst7
				else
					init=${lams[$(($i-1))]}_eqP0TI.rst7
				fi

                        	echo "Running $stage for lambda ${lam}..."

				stage=minTI
				${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${init} -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${init}
                                sleep 1

				laststage=minTI; stage=eqpre1P0TI
                        	${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
				sleep 1

				laststage=eqpre1P0TI; stage=eqpre2P0TI
                        	${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
				sleep 1

				laststage=eqpre2P0TI; stage=eqP0TI
                        	${EXE} -O -p ${top}/unisc.parm7 -c t${trial}/${lam}_${laststage}.rst7 -i inputs/${lam}_${stage}.mdin -o t${trial}/${lam}_${stage}.mdout -r t${trial}/${lam}_${stage}.rst7 -ref t${trial}/${lam}_${laststage}.rst7
				sleep 1

                        	cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin t${trial}/${lam}_${stage}.rst7
autoimage
trajout t${trial}/${lam}_${stage}_centered.rst7
go
quit
EOF2
                         	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                         	cpptraj < center.in
                         	sleep 1
                         	mv t${trial}/${lam}_${stage}_centered.rst7 t${trial}/${lam}_${stage}.rst7
                	done
			laststage=eqP0TI
        	else
                	# check if pmemd.cuda.MPI is present
                	if ! command -v ${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                	export LAUNCH="mpirun -np ${#lams[@]}"
                	export EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
                	export MV2_ENABLE_AFFINITY=0
                	${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/t${trial}_${stage}.groupfile

                	for lam in ${lams[@]};do
                        	cat <<EOF2 > center.in
parm ${top}/unisc.parm7
trajin t${trial}/${lam}_${stage}.rst7
autoimage
trajout t${trial}/${lam}_${stage}_centered.rst7
go
quit
EOF2
                        	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        	cpptraj < center.in
                        	sleep 1
                        	mv t${trial}/${lam}_${stage}_centered.rst7 t${trial}/${lam}_${stage}.rst7
                	done
        	fi
	done

        # run production
        EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
        echo "running replica ti"
        mpirun -np ${#lams[@]} ${EXE} -rem 3 -ng ${#lams[@]} -groupfile inputs/t${trial}_ti.groupfile
done

