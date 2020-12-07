#!/bin/bash

main() {

    if [ -z ${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi

    local base="unisc" 
    local lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )

    #
    # if we don't have any completed production, then there's nothing to analyze
    #

    if [ ! -d production ]; then 
       echo "Nothing to analyze because you don't have a production directory yet"
       exit
    fi

    #
    # the analysis will use a python utility to transform the data in the
    # re-analyzed mdout files to a series of .dat files
    #

    if [ ! -e "production/stdti_step2dats.py" ]; then
       write_python_script
    fi

    #
    # the production dir contains many subdirs that are 0-padded integers.
    # what is the largest number that we can find?
    #

    local last_step=$(for f in $(ls production | grep -v '\.' | tail -n 1); do bc -l <<< $f; done)

    #
    # if we couldn't find a valid subdirectory, then exit now
    #

    if [ "${last_step}" == "" ]; then 
       echo "Nothing to analyze because I couldn't find a valid subdirectory name in production/"
       exit
    fi


    #
    # collect the existing results, if any
    #

    echo "Collecting existing results before submitting new analysis..."
    python2.7 analysis2results.py

    echo ""
    echo ""
    echo "Searching for un-analyzed production directories..."

    local script="analysis.slurm"
    local jarr=()
    for step in $(seq ${last_step}); do

       #
       # this is the zero-padded name
       # 

       local step_name=$(printf "%06i" ${step})

       #
       # if the subdir does not exist, then skip it
       #

       if [ ! -d "production/${step_name}" ]; then
          echo "skipping ${step_name} because dir does not exist"
          continue
       fi


       #
       # do we actually have to analyze this subdir?
       #


       #
       # if it doesn't have mdout files, then I think I already deleted it to save disk space
       #

       local ok=1
       local cnt=0
       for lam in ${lams[@]}; do
           if [ -e "production/${step_name}/${base}_${lam}.mdout" ]; then
               cnt=$(( ${cnt} + 1 ))
           else
               ok=0
           fi
       done
       if [ "${cnt}" == "0" ]; then
          echo "skipping ${step_name} because it was probably already deleted from this machine"
          continue
       elif [ "${ok}" == "0" ]; then
          echo "skipping ${step_name} because it has some, but not all, of the mdout files (something wrong here?)"
          continue
       fi


       ok=1
       for lam in ${lams[@]}; do
           if [ ! -e "production/${step_name}/${base}_${lam}.nc" ]; then
               ok=0
               echo "Missing production/${step_name}/${base}_${lam}.nc"
           fi
       done
       if [ "${ok}" == "0" ]; then
          echo "skipping ${step_name} because it is missing 1-or-more nc files"
          continue
       fi


       #
       # do we have the mbar trace file generated from the analysis
       #

       local testfile=$(printf "production/${step_name}/efep_%.8f_%.8f.dat" 1. 1.)
       if [ -e ${testfile} ]; then
          echo "skipping ${step_name} because it has already been analyzed"
          continue
       fi


       cd production/${step_name}
       python2.7 ../stdti_step2dats.py ${base}_*.mdout
       cd ../../


    done


    python2.7 analysis2results.py


}



##############################################################################


write_template() {
    local fname="$1"
    shift
    local sarr=("$@")
    local nsarr=${#sarr[@]}
    local lsarr=$((${nsarr}-1))
    local sline=${sarr[0]}
    for istep in $(seq ${lsarr}); do sline="${sline},${sarr[${istep}]}"; done

    cat << EOF > production/${fname}
#!/bin/bash
#SBATCH --job-name="${fname}"
#SBATCH --output="${fname}.slurmout"
#SBATCH --error="${fname}.slurmerr"
#SBATCH --array=${sline}
EOF
    cat << 'EOF' >> production/${fname}
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


export EXE="${AMBERHOME}/bin/pmemd.MPI"

prefix="unisc"
lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )
EOF
    cat << 'EOF' >> production/${fname}
istep=${SLURM_ARRAY_TASK_ID}


#
# is the input file asking for a replica-exchange simulation?
#
numexchg=$( grep numexchg ../inputs/${prefix}_${lams[0]}_analyze.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\1/' )
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
   echo "CANNOT RUN REPLICA EXCHANGE IN ANALYSIS-MODE"
   exit 1
fi

#
# each step consists of running all windows for some length of time
#

istep=$(printf "%06i" ${istep})

#
# if don't have the subdirectory, then something's gone wrong
#

if [ ! -d ${istep} ]; then
   exit
fi

cd ${istep}

#
# run analysis for each window
#

for lambda in ${lams[@]}; do
    base="${prefix}_${lambda}"
    rest="${base}_restart"
    init="${base}_initial"
    anal="${base}_analyze"

    ${LAUNCH} ${EXE} -O -c ${base}.rst7 -p ../../${prefix}.parm7 -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -y ${base}.nc

    for tmpfile in ${anal}.rst7 ${anal}.nc ${anal}.mdinfo logfile; do
       if [ -e "${tmpfile}" ]; then
          rm -f "${tmpfile}"
       fi
    done
done

python2.7 ../stdti_step2dats.py *_analyze.mdout

ok=1
for lambda in ${lams[@]}; do
    if [ ! -e "dvdl_${lambda}.dat" ]; then
      ok=0
    fi
done
if [ "${ok}" == "1" ]; then
    for lambda in ${lams[@]}; do
       base="${prefix}_${lambda}"
       anal="${base}_analyze"
       if [ -e "${anal}.mdout" ]; then
          rm -f "${anal}.mdout"
       fi
    done
fi

EOF
    
}




write_python_script() {

cat << 'EOF' > production/stdti_step2dats.py
#!/usr/bin/env python2.7
import sys,os

def extract_traditional_ti( fname, write=False ):
    import os
    from collections import defaultdict as ddict

    fh = open(fname,"r")
    if not fh:
        raise Exception("Could not open %s\n"%(fname))



    numexchg=0
    nstlim=None
    ntpr=None
    dt=None
    irest=0
    for line in fh:
        cmdstr,sepstr,comstr = line.partition("!")
        if "ntpr" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "ntpr":
                    ntpr = int( cols[icol+1] )
                    break
        if "dt" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "dt":
                    dt = float( cols[icol+1] )
                    break
        if "numexchg" in cmdstr:
            cols = cmdstr.replace("=","").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "numexchg":
                    numexchg = int( cols[icol+1] )
                    break
        if "nstlim" in cmdstr:
            cols = cmdstr.replace("="," ").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "nstlim":
                    nstlim = int( cols[icol+1] )
                    break
        if "irest" in cmdstr:
            cols = cmdstr.replace("="," ").replace(",","").strip().split()
            for icol in range(len(cols)-1):
                if cols[icol] == "irest":
                    irest = int( cols[icol+1] )
                    break

    if ntpr is None:
        raise Exception("Could not determine ntpr from %s"%(fname))

    if dt is None:
        raise Exception("Could not determine dt from %s"%(fname))

    if nstlim is None:
        raise Exception("Could not determine nstlim from %s"%(fname))

    if numexchg < 1:
        numexchg = 1

    dt = dt
    nstep_per_sim = nstlim * numexchg
    nframe_per_sim = nstep_per_sim / ntpr

    if nstep_per_sim % ntpr != 0:
        print "num md steps per simulation is not a multiple of ntpr. Unclear how the simulation time works"

    t_per_frame = dt * ntpr
    t_per_sim = t_per_frame * nframe_per_sim


    fh = open(fname,"r")

    
    efeps = []
    dvdls = []
    efep = []
    reading_region_1 = False

    lam = None
    nlam = 0
    
    for line in fh:
        if "A V E R A G E S" in line:
            break
        if "MBAR Energy analysis:" in line:
            efep = []
        if "clambda" in line:
            if lam is None:
                cols = line.replace("="," ").replace(","," ").split()
                for i in range(len(cols)):
                    if cols[i] == "clambda":
                        lam = float(cols[i+1])
                        break
        elif "Energy at " in line:
            #print line
            val = line.strip().split()[-1]
            if "****" in val:
                val = 10000.00
                #if len(efep) > 0:
                #   if efep[-1] < 0:
                #       val = -val
            else:
                val = float(val)
            efep.append( val )
        elif "TI region  1" in line:
            reading_region_1 = True
            dvdl = 0
        elif "| TI region  2" in line:
            #print line
            reading_region_1 = False
            #print dvdl
            dvdls.append( dvdl )
            if len( efep ) > 0:
                efeps.append( efep )
                nlam = len(efep)
        elif "TI region " in line:
            reading_region_1 = False

        if "DV/DL  =" in line and reading_region_1:
            #print line
            cols = line.strip().split()
            #print cols
            dvdl = float( cols[-1] )
            #dvdls.append( float( cols[-1] ) )
            #if len( efep ) > 0:
            #    efeps.append( efep )
            #    nlam = len(efep)
    if write:
        lams = [ float(i) / ( nlam-1. ) for i in range(nlam) ]
        for l in lams:
            if abs(l-lam) < 0.001:
                lam = l
                break
        head, tail = os.path.split(fname)
        dvdl_fname = os.path.join( head, "dvdl_%.8f.dat"%( lam ) )

        if irest == 0:
           dvdls=dvdls[1:]

        fh = file(dvdl_fname,"w")
        for i in range(len(dvdls)):
            fh.write("%.4f %18.6f\n"%((i+1)*t_per_frame,dvdls[i]))
        fh.close()
        for ilam,plam in enumerate(lams):
            efep_fname = os.path.join( head, "efep_%.8f_%.8f.dat"%( lam, plam ) )
            fh = file(efep_fname,"w")
            for i in range(len(efeps)):
                fh.write("%.4f %18.6f\n"%((i+1)*t_per_frame,efeps[i][ilam]))
            fh.close()

    return dvdls,efeps


for arg in sys.argv[1:]:
    if os.path.isfile( arg ):
        if ".mdout" in arg:
            extract_traditional_ti( arg, write=True )
        else:
            print "File does not end in .mdout: %s"%(arg)
    else:
        print "File not found: %s"%(arg)

EOF

chmod u+x production/stdti_step2dats.py

}





# ---------------------------
# call to the main function
# ---------------------------

main

