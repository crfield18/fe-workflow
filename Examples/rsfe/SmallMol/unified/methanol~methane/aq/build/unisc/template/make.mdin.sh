#!/bin/bash

lams=( 0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000 )

prefix="unisc"

mdin=TEMPLATE.mdin



#
# if you set cpu_equil=1, then gti options are removed from equilibration mdin files
#
cpu_equil=0

#
# is this softcore?
#
ifsc="1"

#
# is this parameter interpolated?
#
piti="0"

#
# is this a parameter-interpolated softcore calculation?
#

piscti="0"
if [ "${piscti}" == "1" ]; then
    comment="\!"
fi

noshake=$(grep noshakemask TEMPLATE.mdin | sed -e "s/^.*= *//" -e "s/'//g")


#
# is the input file asking for a replica-exchange simulation?
#

numexchg=$( grep numexchg ${mdin} | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\1/' )
if [ "${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if this is replica exchange AND we are doing piscti, then do acyclic exchange loops
#

acyclic=0
if [ "${piscti}" == "1" -a "${rem}" == "3" ]; then
   acyclic=1
fi

    

if [ ! -d inputs ]; then
   mkdir inputs
fi

for lam in ${lams[@]}; do   
    base=${prefix}_${lam}
    init=inputs/${base}_initial
    rest=inputs/${base}_restart
    anal=inputs/${base}_analyze
    mini=inputs/${base}_minimiz
    heat=inputs/${base}_heating
    pre1=inputs/${base}_press01
    pre2=inputs/${base}_press02
    pre3=inputs/${base}_press03
    pre4=inputs/${base}_press04
    pre5=inputs/${base}_press05
    pre6=inputs/${base}_press06
    pre7=inputs/${base}_press07
    pre8=inputs/${base}_press08


    sed -e "s/CLAMBDA/${lam}/" \
        -e "s/irest *= *[0-9]*/irest = 1/" \
        -e "s/ntx *= *[0-9]*/ntx = 5/" \
        -e "s/acyclic_remd *= *[0-9]*/acyclic_remd = ${acyclic}/" \
        -e "s/clambda/${comment}clambda/" \
        -e "s/ifsc/${comment}ifsc/" \
        -e "s/icfe/${comment}icfe/" \
        -e "s/timask/${comment}timask/" \
        -e "s/scmask/${comment}scmask/" \
        -e "s/crgmask/${comment}crgmask/" \
         ${mdin} > ${rest}.mdin

    sed -e "s/irest *= *[0-9]*/irest = 0/" \
        -e "s/ntx *= *[0-9]*/ntx = 1/" \
         ${rest}.mdin > ${init}.mdin

    sed -e "s/imin *= *[0-9]*/imin = 1/" \
	-e "s/maxcyc *= *[0-9]*/maxcyc = 1000/" \
	-e "s/ntc *= *[0-9]*/ntc = 1/" \
	-e "s/iwrap *= *[0-9]*/iwrap = 0/" \
        -e "s/ntwx *= *[0-9]*/ntwx = 500/" \
	-e "s/ntwr *= *[0-9]*/ntwr = 500/" \
	-e "s/ntpr *= *[0-9]*/ntpr = 500/" \
	-e "s/ntr *= *[0-9]*/ntr = 1/" \
	-e "s/restraintmask *=.*/restraintmask = '!:WAT \& !@H= \& !${noshake}'/" \
	-e "s/restraint_wt *= *[0-9\.]*/restraint_wt = 100./" \
        -e "s/acyclic_remd *= *[0-9]*/acyclic_remd = 0/" \
        -e 's/gremd_acyc/!gremd_acyc/' \
	-e "s/numexchg *= *[0-9]*/numexchg = 0/" \
        -e "s/ifmbar *= *[0-9]*/ifmbar = 0/" \
        ${init}.mdin > ${mini}.mdin

    if [ "${cpu_equil}" == "1" ]; then
       sed -i -e 's/gti/!gti/' ${mini}.mdin
    fi

    sed -e "s/imin *= *[0-9]*/imin = 0/" \
	-e "s/maxcyc *= *[0-9]*/maxcyc = 0/" \
	-e "s/nstlim *= *[0-9]*/nstlim = 100000/" \
	-e "s/tempi *= *[0-9\.]*/tempi = 100./" \
	-e "s/temp0 *= *[0-9\.]*/temp0 = 298./" \
	-e "s/ntc *= *[0-9]*/ntc = 2/" \
        -e "s/ntwx *= *[0-9]*/ntwx = 10000/" \
	-e "s/ntwr *= *[0-9]*/ntwr = 1000/" \
	-e "s/ntpr *= *[0-9]*/ntpr = 1000/" \
	-e "s/nmropt *= *[0-9]*/nmropt = 1/" \
	${mini}.mdin > ${heat}.mdin

    
    cat <<EOF >> ${heat}.mdin

 &wt
     TYPE="TEMP0", istep1=0, istep2=100000,
     value1=100., value2=298.,
 /
 &wt
     TYPE="END",
 /

EOF
    
    if [ "${rem}" == "0" ]; then

    sed -e "s/imin *= *[0-9]*/imin = 0/" \
	-e "s/maxcyc *= *[0-9]*/maxcyc = 0/" \
	-e "s/nstlim *= *[0-9]*/nstlim = 2000/" \
	-e "s/irest *= *[0-9]*/irest = 1/" \
	-e "s/ntx *= *[0-9]*/ntx = 5/" \
	-e "s/tempi *= *[0-9\.]*/tempi = 298./" \
	-e "s/temp0 *= *[0-9\.]*/temp0 = 298./" \
	-e "s/ntc *= *[0-9]*/ntc = 2/" \
        -e "s/ntwx *= *[0-9]*/ntwx = 10000/" \
	-e "s/ntwr *= *[0-9]*/ntwr = 1000/" \
	-e "s/ntpr *= *[0-9]*/ntpr = 1000/" \
	-e "s/ntb *= *[0-9]*/ntb = 2/" \
	-e "s/ntp *= *[0-9]*/ntp = 1/" \
	${mini}.mdin > ${pre1}.mdin

    else


    sed -e "s/imin *= *[0-9]*/imin = 0/" \
	-e "s/maxcyc *= *[0-9]*/maxcyc = 0/" \
	-e "s/nstlim *= *[0-9]*/nstlim = 2000/" \
	-e "s/irest *= *[0-9]*/irest = 1/" \
	-e "s/ntx *= *[0-9]*/ntx = 5/" \
	-e "s/tempi *= *[0-9\.]*/tempi = 298./" \
	-e "s/temp0 *= *[0-9\.]*/temp0 = 298./" \
	-e "s/ntc *= *[0-9]*/ntc = 2/" \
        -e "s/ntwx *= *[0-9]*/ntwx = 10000/" \
	-e "s/ntwr *= *[0-9]*/ntwr = 1000/" \
	-e "s/ntpr *= *[0-9]*/ntpr = 1000/" \
	${mini}.mdin > ${pre1}.mdin

    fi

    sed -e "s/nstlim *= *[0-9\.]*/nstlim = 18000/" \
	${pre1}.mdin > ${pre2}.mdin

    sed -e "s/nstlim *= *[0-9\.]*/nstlim = 25000/" \
	${pre1}.mdin > ${pre3}.mdin

    sed -e "s/nstlim *= *[0-9\.]*/nstlim = 55000/" \
	${pre1}.mdin > ${pre4}.mdin

    sed -e "s/nstlim *= *[0-9\.]*/nstlim = 100000/" \
        -e "s/restraint_wt *= *[0-9\.]*/restraint_wt = 10./" \
	-e "s/irest *= *[0-9]*/irest = 1/" \
	-e "s/ntx *= *[0-9]*/ntx = 5/" \
	${pre4}.mdin > ${pre5}.mdin

    sed -e "s/restraintmask *=.*/restraintmask = '@CA,N,C \& !${noshake}'/" \
	${pre5}.mdin > ${pre6}.mdin

    sed -e "s/restraint_wt *= *[0-9\.]*/restraint_wt = 1./" \
	${pre6}.mdin > ${pre7}.mdin

    sed -e "s/restraint_wt *= *[0-9\.]*/restraint_wt = 0.1/" \
	${pre7}.mdin > ${pre8}.mdin

    sed -e "s/${comment}clambda/clambda/" \
        -e "s/${comment}ifsc/ifsc/" \
        -e "s/${comment}icfe/icfe/" \
        -e "s/${comment}timask/timask/" \
        -e "s/${comment}scmask/scmask/" \
        -e "s/${comment}crgmask/crgmask/" \
        -e "s/imin *= *[0-9]*/imin = 6/" \
        -e "s/ntpr *= *[0-9]*/ntpr = 1/" \
        -e "s/ntwx *= *[0-9]*/ntwx = 0/" \
        -e "s/ntwr *= *[0-9]*/ntwr = 0/" \
        -e "s/icfe *= *[0-9]*/icfe = 1/" \
        -e "s/ifsc *= *[0-9]*/ifsc = ${ifsc}/" \
        -e "s/numexchg *= *[0-9]*/numexchg = 0/" \
        -e "s/acyclic_remd/\!acyclic_remd/" \
        -e "s/chngmask *= *[0-9]*/chngmask = 1/" \
         ${rest}.mdin > ${anal}.mdin
done


truncate -s0 inputs/${prefix}_initial.groupfile
truncate -s0 inputs/${prefix}_restart.groupfile
truncate -s0 inputs/${prefix}_analyze.groupfile
truncate -s0 inputs/${prefix}_minimiz.groupfile
truncate -s0 inputs/${prefix}_heating.groupfile
truncate -s0 inputs/${prefix}_press01.groupfile
truncate -s0 inputs/${prefix}_press02.groupfile
truncate -s0 inputs/${prefix}_press03.groupfile
truncate -s0 inputs/${prefix}_press04.groupfile
truncate -s0 inputs/${prefix}_press05.groupfile
truncate -s0 inputs/${prefix}_press06.groupfile
truncate -s0 inputs/${prefix}_press07.groupfile
truncate -s0 inputs/${prefix}_press08.groupfile

for lam in ${lams[@]}; do
    base=${prefix}_${lam}
    init=${base}_initial
    rest=${base}_restart
    anal=${base}_analyze
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

    parm=${prefix}.parm7
    if [ "${piti}" == "1" ]; then
       parm=${base}.parm7
    fi

    echo "-O -c inputs/${init}.rst7 -p ${parm} -i inputs/${mini}.mdin -o equilib/${mini}.mdout -r equilib/${mini}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7" >> inputs/${prefix}_minimiz.groupfile
    echo "-O -c equilib/${mini}.rst7 -p ${parm} -i inputs/${heat}.mdin -o equilib/${heat}.mdout -r equilib/${heat}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7" >> inputs/${prefix}_heating.groupfile
    echo "-O -c equilib/${heat}.rst7 -p ${parm} -i inputs/${pre1}.mdin -o equilib/${pre1}.mdout -r equilib/${pre1}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref inputs/${init}.rst7" >> inputs/${prefix}_press01.groupfile
    echo "-O -c equilib/${pre1}.rst7 -p ${parm} -i inputs/${pre2}.mdin -o equilib/${pre2}.mdout -r equilib/${pre2}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre1}.rst7" >> inputs/${prefix}_press02.groupfile
    echo "-O -c equilib/${pre2}.rst7 -p ${parm} -i inputs/${pre3}.mdin -o equilib/${pre3}.mdout -r equilib/${pre3}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre2}.rst7" >> inputs/${prefix}_press03.groupfile
    echo "-O -c equilib/${pre3}.rst7 -p ${parm} -i inputs/${pre4}.mdin -o equilib/${pre4}.mdout -r equilib/${pre4}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre3}.rst7" >> inputs/${prefix}_press04.groupfile
    echo "-O -c equilib/${pre4}.rst7 -p ${parm} -i inputs/${pre5}.mdin -o equilib/${pre5}.mdout -r equilib/${pre5}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press05.groupfile
    echo "-O -c equilib/${pre5}.rst7 -p ${parm} -i inputs/${pre6}.mdin -o equilib/${pre6}.mdout -r equilib/${pre6}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press06.groupfile
    echo "-O -c equilib/${pre6}.rst7 -p ${parm} -i inputs/${pre7}.mdin -o equilib/${pre7}.mdout -r equilib/${pre7}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press07.groupfile
    echo "-O -c equilib/${pre7}.rst7 -p ${parm} -i inputs/${pre8}.mdin -o equilib/${pre8}.mdout -r equilib/${pre8}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile -ref equilib/${pre4}.rst7" >> inputs/${prefix}_press08.groupfile


    
    echo "-O -c inputs/${init}.rst7 -p ${parm} -i inputs/${init}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile" >> inputs/${prefix}_initial.groupfile
    echo "-O -c current/${rest}.rst7 -p ${parm} -i inputs/${rest}.mdin -o current/${base}.mdout -r current/${base}.rst7 -x current/${base}.nc -inf current/${base}.mdinfo -l current/${base}.logfile" >> inputs/${prefix}_restart.groupfile
    echo "-O -c ${base}.rst7 -p ../../${parm} -i ../../inputs/${anal}.mdin -o ${anal}.mdout -r ${anal}.rst7 -x ${anal}.nc -inf ${anal}.mdinfo -l ${anal}.logfile -y ${base}.nc " >> inputs/${prefix}_analyze.groupfile

    if [ ! -e "inputs/${init}.rst7" ]; then
       if [ -e "${base}.rst7" ]; then
          cd inputs
          ln -s ../${base}.rst7 ${init}.rst7
          cd ../
       elif [ -e "${prefix}.rst7" ]; then
          cd inputs
          ln -s ../${prefix}.rst7 ${init}.rst7
          cd ../
       fi
    fi
done
