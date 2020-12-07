#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi

    ###################################################################
    local write_template=write_parallel_cpu_template
    #local write_template=write_serial_cpu_template
    #local write_template=write_serial_gpu_template
    ###################################################################
    local template=equilibration.slurm
    ${write_template} "${template}"

    #
    # do we already have stuff running in this directory?
    # if so, then append the current jobs to the end of the dependency chain
    #

    local lastid=$(squeue -h -o "%i %o" | grep ${PWD} | sed 's/ .*//' | sort -rn | head -n 1)
    if [ "${lastid}" != "" ]; then
       sbatch --dependency=afterany:${lastid} ${template}
    else
       sbatch ${template}
    fi

    
}


##############################################################################

write_parallel_cpu_template() {
    local fname="$1"
    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=main
#SBATCH --nodes=21
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=2-00:00:00

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="srun --mpi=pmi2"


export EXE=${AMBERHOME}/bin/pmemd.MPI
prefix="unisc"
lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )

EOF
write_parallel_template "${fname}"
}

write_parallel_template() {
    local fname="$1"
cat << 'EOF' >> ${fname}
if [ ! -d equilib ]; then
    mkdir equilib
fi
if [ ! -d current ]; then
    mkdir current
fi

${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_minimiz.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_heating.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press01.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press02.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press03.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press04.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press05.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press06.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press07.groupfile
${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_press08.groupfile

for lam in ${lams[@]}; do
    base=${prefix}_${lam}
    rest=${base}_restart
    heat=${base}_heating
    pre1=${base}_press01
    pre2=${base}_press02
    pre3=${base}_press03
    pre4=${base}_press04
    pre5=${base}_press05
    pre6=${base}_press06
    pre7=${base}_press07
    pre8=${base}_press08

    cp equilib/${pre8}.rst7 current/${rest}.rst7

    for tmpfile in current/${base}.mdout current/${base}.nc current/${base}.mdinfo current/${base}.logfile; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

EOF
}

##############################################################################
write_serial_cpu_template() {
    local fname="$1"
    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --array=1-21
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=2-00:00:00
##SBATCH --exclude=cuda[001-008]

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="srun --mpi=pmi2"


export EXE=${AMBERHOME}/bin/pmemd.MPI
prefix="unisc"
lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )

export parm=${prefix}.parm7

EOF
write_serial_template ${fname}
}

##############################################################################
write_serial_gpu_template() {
    local fname="$1"
    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --array=1-21
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mincpus=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#
##SBATCH --share
##SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --export=ALL
#SBATCH --time=2-00:00:00

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="srun --mpi=pmi2"


export EXE=${AMBERHOME}/bin/pmemd.cuda
prefix="unisc"
lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )

export parm=${prefix}.parm7

EOF
write_serial_template ${fname}
}

write_serial_template() {
    local fname="$1"
cat <<'EOF' >> ${fname} 
if [ ! -d equilib ]; then
    mkdir equilib
fi
if [ ! -d current ]; then
    mkdir current
fi


myidx=$(( ${SLURM_ARRAY_TASK_ID} - 1 ))
mylam=${lams[${myidx}]}

for lam in ${mylam}; do
    base=${prefix}_${lam}
    init=${base}_initial
    rest=${base}_restart
    mini=${base}_minimiz
    heat=${base}_heating
    pre1=${base}_press01
    pre2=${base}_press02
    pre3=${base}_press03
    pre4=${base}_press04
    pre5=${base}_press05
    pre6=${base}_press06
    pre7=${base}_press07
    pre8=${base}_press08

    ${LAUNCH} ${EXE} -O -c inputs/${init}.rst7 -p ${parm} -i inputs/${mini}.mdin -o equilib/${mini}.mdout -r equilib/${mini}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${mini}.rst7 -p ${parm} -i inputs/${heat}.mdin -o equilib/${heat}.mdout -r equilib/${heat}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${heat}.rst7 -p ${parm} -i inputs/${pre1}.mdin -o equilib/${pre1}.mdout -r equilib/${pre1}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre1}.rst7 -p ${parm} -i inputs/${pre2}.mdin -o equilib/${pre2}.mdout -r equilib/${pre2}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre1}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre2}.rst7 -p ${parm} -i inputs/${pre3}.mdin -o equilib/${pre3}.mdout -r equilib/${pre3}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre2}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre3}.rst7 -p ${parm} -i inputs/${pre4}.mdin -o equilib/${pre4}.mdout -r equilib/${pre4}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre3}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre4}.rst7 -p ${parm} -i inputs/${pre5}.mdin -o equilib/${pre5}.mdout -r equilib/${pre5}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre5}.rst7 -p ${parm} -i inputs/${pre6}.mdin -o equilib/${pre6}.mdout -r equilib/${pre6}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre6}.rst7 -p ${parm} -i inputs/${pre7}.mdin -o equilib/${pre7}.mdout -r equilib/${pre7}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7
    ${LAUNCH} ${EXE} -O -c equilib/${pre7}.rst7 -p ${parm} -i inputs/${pre8}.mdin -o equilib/${pre8}.mdout -r equilib/${pre8}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7


    cp equilib/${pre8}.rst7 current/${rest}.rst7

    for tmpfile in current/${base}.mdout current/${base}.nc current/${base}.mdinfo logfile*; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

EOF
}


# ---------------------------
# call to the main function
# ---------------------------

main

