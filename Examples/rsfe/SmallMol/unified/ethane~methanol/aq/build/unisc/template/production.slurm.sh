#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi


    ###################################################################

    local first_job_gets_to_step=1
    local final_job_gets_to_step=1
    local steps_per_job=1

    CHECKNAN="F"

    local write_template=write_parallel_cpu_template
    #local write_template=write_parallel_gpu_template
    # THE SERIAL-CASES CAN ONLY BE USED FOR NON-REMD SIMULATIONS
    #local write_template=write_serial_cpu_template
    #local write_template=write_serial_gpu_template

    ###################################################################


    local template=production.slurm
    
    local JARR=($( seq ${first_job_gets_to_step} ${steps_per_job} ${final_job_gets_to_step} ))
    if [ "${JARR[${#JARR[@]}-1]}" != "${final_job_gets_to_step}" ]; then JARR+=(${final_job_gets_to_step}); fi

    ${write_template} "${template}" "${CHECKNAN}" "${JARR[@]}"

    local lastid=$(squeue -h -o "%F %o" | grep ${PWD} | sed 's/ .*//' | sort -rn | head -n 1)
    if [ "${lastid}" != "" ]; then
       cmdstr="sbatch --dependency=afterany:${lastid} ${template}"
       line=$( sbatch --dependency=afterany:${lastid} ${template} )
    else
       cmdstr="sbatch ${template}"
       line=$( sbatch ${template} )
    fi
    echo "${cmdstr}"
    echo "${line}"
}



##############################################################################

write_parallel_cpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=main
#SBATCH --nodes=21
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
REVERSE=false
if [ "${REVERSE}" = true ]; then lams=( 1.00000000 0.95000000 0.90000000 0.85000000 0.80000000 0.75000000 0.70000000 0.65000000 0.60000000 0.55000000 0.50000000 0.45000000 0.40000000 0.35000000 0.30000000 0.25000000 0.20000000 0.15000000 0.10000000 0.05000000 0.00000000 ); fi

EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_parallel_body "${fname}"
}

##############################################################################
write_serial_cpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%1
EOF
    cat << 'EOF' >> ${fname}
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
REVERSE=false
if [ "${REVERSE}" = true ]; then lams=( 1.00000000 0.95000000 0.90000000 0.85000000 0.80000000 0.75000000 0.70000000 0.65000000 0.60000000 0.55000000 0.50000000 0.45000000 0.40000000 0.35000000 0.30000000 0.25000000 0.20000000 0.15000000 0.10000000 0.05000000 0.00000000 ); fi

EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_serial_body "${fname}"
}

##############################################################################

write_parallel_gpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=gpu
#
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=3
##SBATCH --mincpus=3
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=1
#
#SBATCH --share
##SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=2-00:00:00

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="srun --mpi=pmi2"


export EXE=${AMBERHOME}/bin/pmemd.cuda.MPI
prefix="unisc"
lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )
REVERSE=false
if [ "${REVERSE}" = true ]; then lams=( 1.00000000 0.95000000 0.90000000 0.85000000 0.80000000 0.75000000 0.70000000 0.65000000 0.60000000 0.55000000 0.50000000 0.45000000 0.40000000 0.35000000 0.30000000 0.25000000 0.20000000 0.15000000 0.10000000 0.05000000 0.00000000 ); fi

EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF
write_parallel_body "${fname}"
}

##############################################################################
write_serial_gpu_template() {
    local fname="$1"
    shift
    local CHECKNAN="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > ${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}%1
EOF
    cat << 'EOF' >> ${fname}
#SBATCH --partition=gpu
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --mincpus=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#
#SBATCH --share
##SBATCH -C gpu_shared # set this on xstream.stanford.xsede.org INSTEAD of --share
#
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=2-00:00:00

export MV2_ENABLE_AFFINITY=0
source ${AMBERHOME}/amber.sh
export LAUNCH="srun --mpi=pmi2"

export EXE=${AMBERHOME}/bin/pmemd.cuda
prefix="unisc"
lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )
REVERSE=false
if [ "${REVERSE}" = true ]; then lams=( 1.00000000 0.95000000 0.90000000 0.85000000 0.80000000 0.75000000 0.70000000 0.65000000 0.60000000 0.55000000 0.50000000 0.45000000 0.40000000 0.35000000 0.30000000 0.25000000 0.20000000 0.15000000 0.10000000 0.05000000 0.00000000 ); fi

EOF
   cat << EOF >> ${fname}
SARR=(${sarr[@]})
CHECKNAN="${CHECKNAN}"
EOF
   cat << 'EOF' >> ${fname}
#laststep=${SARR[$(((${SLURM_ARRAY_TASK_ID}-1)))]}
laststep=${SLURM_ARRAY_TASK_ID}
EOF

write_serial_body "${fname}"
}

##############################################################################
write_parallel_body() {
   local fname="$1"
   cat << 'EOF' >> ${fname}
#
# is the input file asking for a replica-exchange simulation?
#

numexchg=$( grep numexchg inputs/${prefix}_${lams[0]}_initial.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#


if [ ! -d current ]; then
   mkdir current
fi


#
# each step consists of running all windows for some length of time
#

for istep in $(seq ${laststep}); do

   istep=$(printf "%06i" ${istep})

   #
   # if we've already run this step, then skip to the next step
   #

   if [ -d production/${istep} ]; then
      continue
   fi

   #
   # this infinite-loop will exit if all simulations run correctly
   #

   while true; do

      #
      # run all simulations at the same time using a groupfile.
      # there are different mdin and rst7 files depending on if we are
      # starting or restarting the simulation
      #

      if [ -e "current/${prefix}_${lams[0]}_restart.rst7" ]; then

         ${LAUNCH} ${EXE} -rem ${rem} -ng ${#lams[@]} -groupfile inputs/${prefix}_restart.groupfile

      else 

         echo "Missing current/${prefix}_${lams[0]}_restart.rst7"
         echo "You should first run equilibration.slurm.sh"
         exit 1

#         if [ ! -d equilib ]; then
#            mkdir equilib
#         fi
#         ${LAUNCH} ${EXE} -ng ${#lams[@]} -groupfile inputs/${prefix}_minimiz.groupfile
#         for lam in ${lams[@]}; do
#             initial=${prefix}_${lam}_initial
#             minimiz=${prefix}_${lam}_minimiz
#             if [ -L inputs/${initial}.rst7 ]; then
#                  rm -f inputs/${initial}.rst7
#             fi
#             cp equilib/${minimiz}.rst7 inputs/${initial}.rst7
#         done
#         ${LAUNCH} ${EXE} -rem ${rem} -ng ${#lams[@]} -groupfile inputs/${prefix}_initial.groupfile

      fi

      ok=1
      for lambda in ${lams[@]}; do
           base="${prefix}_${lambda}"
           if [ ! -e current/${base}.mdout ]; then
              ok=0
           fi
      done
      if [ "${ok}" == "0" ]; then
         echo "There was an unrecoverable error. A mdout file was not produced"
         exit 1
      else
	  incomplete=($(grep -L '5.  TIMINGS' current/${prefix}_*.mdout))
	  if [ "${#incomplete[@]}" -gt 0 ]; then
	      echo "The jobs did not complete for an unknown reason"
	      exit 1
	  fi
      fi

      #
      # if there aren't any not-a-numbers in the mdout files,
      # then the simulations must have run correctly,
      # so copy the output files to a subdirectory and exit
      # the infinite loop
      #

      if [ "$(grep -l NaN current/*.mdout)" == "" -o "${CHECKNAN}" == "F" ]; then

         mkdir -p production/${istep}
         for lambda in ${lams[@]}; do
            base="${prefix}_${lambda}"

            #
            # save the coordinates to restart the next simulation
            #

            cp current/${base}.rst7 current/${base}_restart.rst7

            # 
            # now copy all the outputs to a subdirectory
            #

            mv current/${base}.nc current/${base}.rst7 current/${base}.mdout current/${base}.mdinfo production/${istep}/
            for dumpave in current/*.dumpave; do
               if [ -e "${dumpave}" ]; then
                  mv ${dumpave} production/${istep}/$(basename ${dumpave})
               fi
            done
            if [ -e "current/${base}.logfile" ]; then
               mv current/${base}.logfile production/${istep}/
            fi
         done
         if [ -e "rem.log" ]; then
            mv rem.log production/${istep}/
         fi

         # cd production/${istep}
         # ln -s ../../*.parm7 .
         # cd ../../

         break # exit the infinite loop

     fi
  done
done

EOF
    
}

##############################################################################
write_serial_body() {
local fname="$1"
cat << 'EOF' >> ${fname}

#
# is the input file asking for a replica-exchange simulation?
#

numexchg=$( grep numexchg inputs/${prefix}_${lams[0]}_initial.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#

if [ "${rem}" != "0" ]; then
   echo "CANNOT RUN REPLICA EXCHANGE IN SERIAL-MODE"
   exit 1
fi


if [ ! -d current ]; then
    mkdir current
fi


#
# each step consists of running all windows for some length of time
#

for istep in $(seq ${laststep}); do

   istep=$(printf "%06i" ${istep})

   #
   # if we've already run this step, then skip to the next step
   #

   if [ -d production/${istep} ]; then
      continue
   fi


   #
   # run a simulation for each window
   #

   for lambda in ${lams[@]}; do
       base="${prefix}_${lambda}"
       rest="${base}_restart"
       init="${base}_initial"
       mini="${base}_minimiz"
       parm="${prefix}.parm7"


       #
       # do we need to run this window?
       #

       needed=1
       if [ -e "current/${base}.mdout" ]; then
          if [ "$(grep -LE 'Run   done at|Final Performance Info' current/${base}.mdout)" == "" ]; then
             needed=0
          fi
       fi

       #
       # if we need this window, then run it until it until it terminates normally
       #

       if [ "${needed}" == "1" ]; then
          while true; do

             #
             # are we starting or restarting a simulation?
             # we need to use different mdin and rst7 files depending on the answer
             #

             if [ -e "current/${rest}.rst7" ]; then

                ${LAUNCH} ${EXE} -O -c current/${rest}.rst7 -p ${parm} -i inputs/${rest}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo

             else 

                echo "Missing current/${rest}.rst7"
                echo "You should first run equilibration.slurm.sh"
                exit 1

#                if [ ! -d equilib ]; then
#                    mkdir equilib
#                fi
#                ${LAUNCH} ${EXE} -O -c inputs/${init}.rst7 -p ${parm} -i inputs/${mini}.mdin -o current/${base}.mdout -r equilib/${mini}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile
#                if [ -L inputs/${init}.rst7 ]; then
#                    rm -f inputs/${init}.rst7
#                fi
#                cp equilib/${mini}.rst7 inputs/${init}.rst7
#                ${LAUNCH} ${EXE} -O -c inputs/${init}.rst7 -p ${parm} -i inputs/${init}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo

             fi

             #
             # if we don't get not-a-numbers, then it must've worked, so leave this loop
             #

             if [ ! -e current/${base}.mdout ]; then
                echo "There was an unrecoverable error. A mdout file was not produced"
                exit 1
             else
                incomplete=($(grep -L '5.  TIMINGS' current/${base}.mdout))
                if [ "${#incomplete[@]}" -gt 0 ]; then
                   echo "The jobs did not complete for an unknown reason"
                   exit 1
                fi
             fi 

             if [ "$(grep -l NaN current/${base}.mdout)" == "" -o "${CHECKNAN}" == "F" ]; then
                break
             fi
          done
       fi
   done

   #
   # we've run all the windows for this step, so store the outputs in a subdirectory
   #

   mkdir -p production/${istep}
   for lambda in ${lams[@]}; do
      base="${prefix}_${lambda}"

      #
      # save these coordinates to restart the next step
      # 

      cp current/${base}.rst7 current/${base}_restart.rst7

      # 
      # copy the outputs
      #

      mv current/${base}.nc current/${base}.rst7 current/${base}.mdout current/${base}.mdinfo production/${istep}/
      for dumpave in current/*.dumpave; do
          if [ -e "${dumpave}" ]; then
              mv ${dumpave} production/${istep}/$(basename ${dumpave})
          fi
      done
      if [ -e "current/${base}.logfile" ]; then
         mv current/${base}.logfile production/${istep}/
      fi

   done
   if [ -e "rem.log" ]; then
      mv rem.log production/${istep}/
   fi

   if [ -e logfile ]; then
      rm -f logfile
   fi

   # cd production/${istep}
   # ln -s ../../*.parm7 .
   # cd ../../

done

EOF
    
}


# ---------------------------
# call to the main function
# ---------------------------

main

