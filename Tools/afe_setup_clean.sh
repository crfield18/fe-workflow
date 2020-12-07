#!/bin/bash

##########################################
path=`pwd`
###########################################

##########################################
# functions

##########################################
# read input file
function read_input {

while read line; do
        varlist=(system translist nlambda protocol mapping ntrials cutoff mincyc nstlimnvt nstlimnpt repex nstlimti numexchgti hmr scalpha scbeta gti_add_sc gti_scale_beta gti_cut gti_cut_sc_on gti_cut_sc_off gti_lam_sch gti_ele_sc gti_vdw_sc gti_cut_sc gti_ele_exp gti_vdw_exp stage ticalc partition nnodes ngpus wallclock) 
	IFS=$'\t| |=' read -ra args <<< $line
        if [[ "${args[0]}" =~ ^#.* ]]; then continue; fi
        keyword=${args[0]}; value=${args[1]}
        if [ "${args[0]}" == "translist" ]; then
                IFS='[=()]' read -ra args <<< $line
                keyword=${args[0]}; value=${args[2]}
        fi
        for var in ${varlist[@]}; do
                declare -n arr="$var"
                if [ "$var" == "$keyword" ]; then
                        if [ "$var" == "translist" ]; then
                                arr=($value)
                        else
                                arr=$value
                        fi
                fi
        done


done < $1
}
##########################################


##########################################
# parse input from user
function parse_input {
        if [ ! -z "$1" ]; then
                if [ "$1" == "-h" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ] ; then
                        cat << EOFN > input.template
system=cdk2                     # system
translist=(1h1q~1h1s 1h1q~1h1r 1h1q~1oiy 1h1q~1oi9 1oi9~1h1s 1oi9~1oiy 1h1q~1oiu 1oiu~1h1r)               # list of transformations format:ligandA~ligandB represents transformation ligandA->ligandB
nlambda=4                       # number of lambda windows
protocol=unified                # unified protocol for TI
mapping=MCSS                    # "MCSS/manual/checked" manual if atom mapping between state A and state B needs to be manually inspected. set to 'default' if manual inspection is not needed, or 'checked' is manual inspection has already been done after running the script with mapping=manual
ntrials=3

cutoff=8                        # non-bonded cutoff
mincyc=2000                     # max minimization cycles in lambda window equilibration
nstlimnvt=298000                # length of 0->298K NVT heating of lambda windows
nstlimnpt=500000                # length of NPT equilibration of lambda windows. Ingored if ticalc is set to 'abfe'
repex=true                      # true/false corresponding to use of replica exchange in TI simulations
nstlimti=5000                   # length of TI simulations
numexchgti=1000                 # number of exchanges in replica exchange TI simulations. if repex=false, numexchgti is ignored
hmr=false                       # "true/false" If hmr=true, dt (timestep) is set to 4fs
scalpha=0.5                     # scalpha
scbeta=1.0                      # scbeta
gti_add_sc=5                         # gti_add_sc
gti_scale_beta=1                       # gti_scale_beta
gti_cut=1			# gti_cut
gti_cut_sc_on=6                       # gti_cut_sc_on
gti_cut_sc_off=8                      # gti_cut_sc_off
gti_lam_sch=1                     # gti_lam_sch
gti_ele_sc=1                      # gti_ele_sc
gti_vdw_sc=1                      # gti_vdw_sc
gti_cut_sc=2                      # gti_cut_sc
gti_ele_exp=2                     # gti_ele_exp
gti_vdw_exp=2                     # gti_vdw_exp

ticalc=rbfe                     # "rbfe -> relative binding free energy/rsfe -> relative solvation free energy"
stage=setup                     # "setup/run-equil/check-equil/run-TI/check-TI"

# job submission related
partition=v100                  # name of specific partition on HPC. Use "null" is not relevant
nnodes=1                        # number of nodes to be used for each transformation
ngpus=4                         # number of gpus/node to be used for each transformation
wallclock=24:00:00              # wallclock for individual jobs


EOFN
                        echo "Script expects a file named \"input\" in working directory. Check input.template for details."
                        exit 0
                fi
        else
                if [ ! -f input ]; then echo "Script expects a file named \"input\" in working directory. Run script with -h/-help to generate template input file" && exit 0; fi
        fi

        # read input file
        read_input input


        # check if AMBERHOME is set
        if [ -z "${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        # check if parmed is present
        if ! command -v parmed &> /dev/null;  then echo "parmed is missing." && exit 0; fi

        # check if input directories and files are present
        if [ ! -d initial/${system}/${dir1} ] || [ ! -d initial/${system}/${dir2} ];
                then echo "initial/${system}/${dir1} or initial/${system}/${dir2} folder(s) missing" && exit 0
        else
                for s in ${dir1} ${dir2}; do
                        cd initial/${system}/$s
                                for i in "${!translist[@]}";do
                                        stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                                        if [ ! -f ${stA}_${s}.parm7 ] || [ ! -f ${stA}_${s}.rst7 ] || [ ! -f ${stB}_${s}.parm7 ] || [ ! -f ${stB}_${s}.rst7 ]; then
                                                echo "one or more of ${stA}/${stB} parm/rst files are missing in initial/${system}/${s}"
                                        fi
                                done
                        cd $path
                done
        fi

        # check input file parameters
        if [ "${protocol}" != "unified" ]; then echo "Script currently supports only \"unified\" protocol" && exit 0; fi
        if [ "${mapping}"  != "MCSS" ] && [ "${mapping}" != "manual" ] && [ "${mapping}" != "checked" ];then  echo "\"mapping\" should be set to \"MCSS\", \"manual\", or \"checked\"" && exit 0; fi
        if [ "${repex}"    != "true" ] && [ "${repex}"   != "false"  ]; then echo "\"repex\" should be set to \"true\" or \"false\"" && exit 0; fi
        if [ "${hmr}"      != "true" ] && [ "${hmr}"     != "false"  ]; then echo "\"hmr\" should be set to \"true\" or \"false\""   && exit 0; fi
	if [ "${gti_add_sc}"    -lt 1 ] || [ "${gti_add_sc}"       -gt 5 ]; then echo "Acceptable values for \"gti_add_sc\" are 1, 2(Recommended), 3, 4, and 5" && exit 0; fi
        if [ "${gti_lam_sch}" != 0 ] && [ "${gti_lam_sch}" != 1 ]; then echo "Acceptable values of \"gti_lam_sch\" are 0(default) and 1" && exit 0; fi
        if [ "${gti_lam_sch}" -eq 1 ]; then
                if [ "${gti_ele_sc}" != 1 ] || [ "${gti_vdw_sc}" != 1 ]; then echo "If gti_lam_sch is set to 1, both gti_ele_sc and gti_vdw_sc should be set to 1" && exit 0; fi
        fi
        if [ "${gti_scale_beta}" != 0 ] && [ "${gti_scale_beta}" != 1 ]; then echo "Acceptable values of \"gti_scale_beta\" are 0(default) and 1" && exit 0; fi
        if [ "${gti_scale_beta}" -eq 1 ]; then
                if [ -z "${gti_ele_exp}" ] || [ -z "${gti_vdw_exp}"  ] || [ -z "$scalpha" ] || [ -z "$scbeta" ]; then echo "If gti_scale_beta is set to 1, gti_ele_exp gti_vdw_exp scalpha scbeta all needs to be defined. Recommended values are gti_lam_sch = 1,  gti_ele_sc = 1, gti_vdw_sc = 1, gti_ele_exp = 2, gti_vdw_exp = 2,  gti_scale_beta = 1, scalpha = 0.5, scbeta = 1.0" && exit 0; fi
        fi
        if [ "${gti_cut}" != 0 ] && [ "${gti_cut}" != 1 ]; then echo "Acceptable values of \"gti_cut\" are 0(default) and 1" && exit 0; fi
        if [ "${gti_cut_sc}" != 0 ] && [ "${gti_cut_sc}" != 1 ] && [ "${gti_cut_sc}" != 2 ]; then echo "Acceptable values of \"gti_cut_sc\" are 0(default), 1, 2" && exit 0; fi
        if [ "${gti_cut_sc}" -eq 1 ] || [ "${gti_cut_sc}" -eq 2 ]; then
                if [ -z "${gti_cut_sc_on}" ] || [ -z "${gti_cut_sc_off}" ]; then echo "if \"gti_cut_sc_on\" is set to 1 or 2, gti_cut_sc_on and gti_cut_sc_off must be defined." && exit 0; fi
        fi
        if [ "${cutoff}" -lt "${gti_cut_sc_on}" ] || [ "${cutoff}" -lt "${gti_cut_sc_off}" ] || [ "${gti_cut_sc_off}" -lt "${gti_cut_sc_on}" ]; then echo "Should be \"cutoff\" >= \"gti_cut_sc_off\" > \"gti_cut_sc_on\" " && exit 0; fi



}

##########################################

##########################################
# generate equally spaced lambda windows
function gen_lambdas {
        cat <<EOF > gen_lambda.py
#!/usr/bin/env python

import numpy as np
import sys

nlam = sys.argv[1]; nlam = float (nlam)
a = np.linspace(0, 1, num = nlam, dtype = float)
for x in range(len(a)):
    print("{:.8f}".format(a[x])),

EOF

local nlambda=$1
local lams=$(python2.7 gen_lambda.py ${nlambda})
echo $lams
}

###########################################
#write TEMPLATE
#for relative binding free energies -> [ protein-lig complex - lig in water ]
#################################
function writetemplate_rbfe
{
        #echo "$#"
        local args=$*; args=($args)
        local varlist=(nlambda CUTOFF MINCYC NSTLIMNVT NSTLIMNPT REPEX NSTLIMTI NUMEXCHGTI TIMASK1 TIMASK2 SCMASK1 SCMASK2 SCALPHA SCBETA GTISC GTIBETA GTICUT GTISCON GTISCOFF GTILAMSCH GTISCELE GTISCVDW GTISCCUT GTIEXPELE GTIEXPVDW stA stB)
        local i=0
        for a in "${varlist[@]}"; do
                declare -n arr="$a"
                arr=${args[$i]}
                i=$(($i+1))
        done

lams=($(gen_lambdas $nlambda))
cat <<EOFN >TEMPLATE.sh
#!/bin/bash
top=\`pwd\`

lams=(${lams[@]})

#generate initial configurations for lambda windows
mkdir -p current inputs
truncate -s0 inputs/heat.groupfile inputs/npt.groupfile inputs/ti.groupfile
for lam in \${lams[@]};do
        cp unisc/template/unisc.rst7  current/\${lam}_init.rst7
        cp unisc/template/unisc.parm7 unisc.parm7

        init=\${lam}_init
        mini=\${lam}_mini
        heat=\${lam}_heat
        npt=\${lam}_npt
        ti=\${lam}_ti

        cat<<EOF>inputs/\${lam}_mini.mdin
&cntrl
imin 		= 1
maxcyc 		= ${MINCYC}
drms 		= 0.1 
dx0 		= 0.005 
ntmin 		= 2
irest 		= 0
ntx 		= 1
ntxo 		= 1
ntc 		= 1  
ntf 		= 1
tishake 	= 1
ntwx 		= 500
ntwr 		= 500
ntpr 		= 500
cut 		= ${CUTOFF}
iwrap 		= 1
ig 		= -1

ntb 		= 1 

ifsc 		= 1
icfe 		= 1

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
EOF


        cat<<EOF>inputs/\${lam}_heat.mdin
&cntrl
imin 		= 0
nstlim 		= ${NSTLIMNVT}
dt 		= 0.001
irest 		= 0
ntx 		= 1
ntxo 		= 1
ntc 		= 1 
ntf 		= 1 
tishake 	= 1 
ntwx 		= 10000
ntwr 		= 5000
ntpr 		= 1000
cut 		= ${CUTOFF}
iwrap 		= 1

ntb 		= 1 
tempi 		= 5.
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 	

ifsc 		= 1
icfe 		= 1

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
 &ewald
 /
&wt
type = 'TEMP0',
istep1 = 0, istep2 = 298000,
value1 = 0, value2 = 298,
/
&wt type = 'END'
/
EOF

        cat<<EOF>inputs/\${lam}_npt.mdin
&cntrl
imin 		= 0
nstlim 		= ${NSTLIMNPT}
dt 		= 0.001
irest 		= 1
ntx 		= 5
ntxo 		= 1
ntc 		= 1  
ntf 		= 1 
tishake 	= 1 
ntwx 		= 10000
ntwr 		= 5000
ntpr 		= 1000
cut 		= ${CUTOFF}
iwrap 		= 1

ntb 		= 2 
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1
ntp 		= 1
barostat 	= 2
pres0 		= 1.01325
taup 		= 5.0

ifsc 		= 1
icfe 		= 1

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
 &ewald
 /
EOF

        cat<< EOF >inputs/\${lam}_ti.mdin
&cntrl
imin 		= 0
nstlim 		= ${NSTLIMTI}
numexchg 	= ${NUMEXCHGTI}
dt 		= 0.001
irest 		= 1
ntx 		= 5
ntxo 		= 1
ntc 		= 1
ntf 		= 1
tishake 	= 1 
ntwx 		= 10000
ntwr 		= 5000
ntpr 		= 1000
cut 		= ${CUTOFF}
iwrap 		= 1

ntb 		= 2 
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 
ntp 		= 1
barostat 	= 2
pres0 		= 1.01325
taup 		= 5.0

ifsc 		= 1
icfe 		= 1

ifmbar 		= 1
bar_intervall   = 1
mbar_states 	= \${#lams[@]}
mbar_lambda 	= \${lams[@]}

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}

gremd_acyc 	= 1 
/
 &ewald
 /
EOF

        cat<< EOF >inputs/\${lam}_analyze.mdin
&cntrl
imin 		= 6
nstlim 		= ${NSTLIMTI}
numexchg 	= ${NUMEXCHGTI}
dt 		= 0.001
irest 		= 1
ntx 		= 5
ntxo 		= 1
ntc 		= 1 
ntf 		= 1
tishake 	= 1
ntwx 		= 0
ntwr 		= 0
ntpr 		= 1
cut 		= ${CUTOFF}
iwrap 		= 0

ntb 		= 1 
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 
ntp 		= 1
barostat 	= 2
pres0 		= 1.01325
taup	 	= 1.0

ifsc 		= 1
icfe 		= 1

ifmbar 		= 1
mbar_states 	= \${#lams[@]}
mbar_lambda 	= \${lams[@]}
clambda 	= \${lam}

timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
 &ewald
 /

EOF

        cat<<EOF>> inputs/heat.groupfile
-O -p unisc.parm7 -c current/\${mini}.rst7 -i inputs/\${heat}.mdin -o current/\${heat}.mdout -r current/\${heat}.rst7 -x current/\${heat}.nc -inf current/\${heat}.info 
EOF
        cat<<EOF>> inputs/npt.groupfile
-O -p unisc.parm7 -c current/\${heat}.rst7 -i inputs/\${npt}.mdin -o current/\${npt}.mdout -r current/\${npt}.rst7 -x current/\${npt}.nc -inf current/\${npt}.info 
EOF
        cat<<EOF>> inputs/ti.groupfile
-O -p unisc.parm7 -c current/\${npt}.rst7 -i inputs/\${ti}.mdin -o current/\${ti}.mdout -r current/\${ti}.rst7 -x current/\${ti}.nc -inf current/\${lam}.info 
EOF

done

#submit independent jobs
cat<<EOF > equil.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --share
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=${wallclock}

lams=(\${lams[@]})

top=\\\${PWD}
# check if AMBERHOME is set
if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi
# check if cpptraj is present
if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi


EXE=\\\${AMBERHOME}/bin/pmemd.cuda
for lam in \\\${lams[@]};do
    	init=\\\${lam}_init
    	mini=\\\${lam}_mini
    	heat=\\\${lam}_heat
    	npt=\\\${lam}_npt

    	echo "running \\\${lam} mini"
    	\\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${init}.rst7 -i inputs/\\\${mini}.mdin -o current/\\\${mini}.mdout -r current/\\\${mini}.rst7 -x current/\\\${mini}.nc -inf current/\\\${mini}.info -ref current/\\\${init}.rst7 
	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${mini}.rst7
autoimage
trajout current/\\\${mini}_centered.rst7
go
quit
EOF2
       	# check if cpptraj is present
       	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
       	cpptraj < center.in
       	sleep 1
       	mv current/\\\${mini}_centered.rst7 current/\\\${mini}.rst7
       	sleep 1
	
    	echo "running \\\${lam} heat"
    	\\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${mini}.rst7 -i inputs/\\\${heat}.mdin -o current/\\\${heat}.mdout -r current/\\\${heat}.rst7 -x current/\\\${heat}.nc -inf current/\\\${heat}.info 
        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${heat}.rst7
autoimage
trajout current/\\\${heat}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/\\\${heat}_centered.rst7 current/\\\${heat}.rst7
        sleep 1

    	echo "running \\\${lam} npt"
    	\\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${heat}.rst7 -i inputs/\\\${npt}.mdin -o current/\\\${npt}.mdout -r current/\\\${npt}.rst7 -x current/\\\${npt}.nc -inf current/\\\${press}.info 
        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${npt}.rst7
autoimage
trajout current/\\\${npt}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/\\\${npt}_centered.rst7 current/\\\${npt}.rst7
        sleep 1


done
EOF

# submit group-ed jobs
cat<<EOF > equilgroup.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --time=${wallclock}

top=\\\${PWD}
lams=(\${lams[@]})
steps=(mini heat npt)
if [ ! -d current ];then mkdir current; fi

for step in \\\${steps[@]}; do
        if [ "\\\$step" == "mini" ]; then
                for lam in \\\${lams[@]};do
                        init=\\\${lam}_init
                        mini=\\\${lam}_mini
			# check if AMBERHOME is set
			if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                        export LAUNCH="srun"
                        export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

                        \\\${LAUNCH} \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${init}.rst7 -i inputs/\\\${mini}.mdin -o current/\\\${mini}.mdout -r current/\\\${mini}.rst7 -x current/\\\${mini}.nc -inf current/\\\${mini}.info -ref current/\\\${init}.rst7 
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${mini}.rst7
autoimage
trajout current/\\\${mini}_centered.rst7
go
quit
EOF2
		        # check if cpptraj is present
        		if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${mini}_centered.rst7 current/\\\${mini}.rst7
                done
        else
		# check if AMBERHOME is set
		if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                export LAUNCH="mpirun -np \\\${#lams[@]}"
                export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                export MV2_ENABLE_AFFINITY=0
                \\\${LAUNCH} \\\${EXE} -ng \\\${#lams[@]} -groupfile inputs/current_\\\${step}.groupfile

                for lam in \\\${lams[@]};do
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${step}.rst7
autoimage
trajout current/\\\${lam}_\\\${step}_centered.rst7
go
quit
EOF2
        		if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${lam}_\\\${step}_centered.rst7 current/\\\${lam}_\\\${step}.rst7
                done
        fi

done

EOF


truncate -s0 prod.slurm
cat<<EOF > prod.slurm
#!/bin/bash
#SBATCH --job-name="pr_${stA}~${stB}.slurm"
#SBATCH --output="pr_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="pr_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --time=${wallclock}

lams=(\${lams[@]})
# check if AMBERHOME is set
if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
echo "running replica ti"
mpirun -np \\\${#lams[@]} \\\${EXE} -rem 3 -ng \\\${#lams[@]} -groupfile inputs/ti.groupfile
EOF

if [ "\${REPEX}" == "false" ]; then
        sed -i -e 's/numexchg/!numexchg/g' -e 's/gremd_acyc = 1/!gremd_acyc = 1/g' inputs/\\\${lam}_ti.mdin inputs/\\\${lam}_analyze.mdin
        sed -i 's/ -rem 3//g' inputs/ti.groupfile
fi

EOFN
}

#################################


###########################################
#write TEMPLATE
#for relative solvation free energies -> [ lig in water - lig in vacuum ]
###########################################

function writetemplate_rsfe
{
        local args=$*; args=($args)
        local varlist=(nlambda CUTOFF MINCYC NSTLIMNVT NSTLIMNPT REPEX NSTLIMTI NUMEXCHGTI TIMASK1 TIMASK2 SCMASK1 SCMASK2 SCALPHA SCBETA GTISC GTIBETA GTICUT GTISCON GTISCOFF GTILAMSCH GTISCELE GTISCVDW GTISCCUT GTIEXPELE GTIEXPVDW stA stB)
        local i=0
        for a in "${varlist[@]}"; do
                declare -n arr="$a"
                arr=${args[$i]}
                i=$(($i+1))
        done

lams=($(gen_lambdas $nlambda))
cat <<EOFN >TEMPLATE.sh
#!/bin/bash
top=\`pwd\`

lams=(${lams[@]})

#generate initial configurations for lambda windows
mkdir -p current inputs
truncate -s0 inputs/heat.groupfile inputs/npt.groupfile inputs/ti.groupfile
for lam in \${lams[@]};do
        cp unisc/template/unisc.rst7  current/\${lam}_init.rst7
        cp unisc/template/unisc.parm7 unisc.parm7

        init=\${lam}_init
        mini=\${lam}_mini
        heat=\${lam}_heat
        ti=\${lam}_ti

        cat<<EOF>inputs/\${lam}_mini.mdin
&cntrl
imin 		= 1
maxcyc 		= ${MINCYC}
drms 		= 0.1
dx0 		= 0.005
ntmin 		= 2
irest 		= 0
ntx 		= 1
ntxo 		= 1
ntc 		= 1 
ntf 		= 1 
tishake 	= 1 
ntwx 		= 500
ntwr 		= 500
ntpr 		= 500
cut 		= ${CUTOFF}
iwrap 		= 1
ig 		= -1

ntb             = 1 

ifsc 		= 1
icfe 		= 1

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
EOF


        cat<<EOF>inputs/\${lam}_heat.mdin
&cntrl
imin 		= 0
nstlim 		= ${NSTLIMNVT}
dt 		= 0.001
irest 		= 0
ntx 		= 1
ntxo 		= 1
ntc 		= 1 
ntf 		= 1 
tishake 	= 1 
ntwx 		= 10000
ntwr 		= 5000
ntpr 		= 1000
cut 		= ${CUTOFF}
iwrap 		= 1

ntb		= 1 
tempi 		= 5.
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 
ntp 		= 0 

clambda 	= \${lam}
ifsc 		= 1
icfe 		= 1

timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
 &ewald
 /
&wt
type = 'TEMP0',
istep1 = 0, istep2 = 298000,
value1 = 0, value2 = 298,
/
&wt type = 'END'
/
EOF


        cat<< EOF >inputs/\${lam}_ti.mdin
&cntrl
imin 		= 0
nstlim 		= ${NSTLIMTI}
numexchg 	= ${NUMEXCHGTI}
dt 		= 0.001
irest 		= 1
ntx 		= 5
ntxo 		= 1
ntc 		= 1
ntf 		= 1
tishake 	= 1 
ntwx 		= 10000
ntwr 		= 5000
ntpr 		= 1000
cut 		= ${CUTOFF}
iwrap 		= 1



ntb 		= 1 
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 
ntp 		= 0

ifsc 		= 1
icfe 		= 1

ifmbar 		= 1
bar_intervall 	= 1
mbar_states 	= \${#lams[@]}
mbar_lambda 	= \${lams[@]}

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
gremd_acyc 	= 1 
/
 &ewald
 /
EOF

        cat<< EOF >inputs/\${lam}_analyze.mdin
&cntrl
imin 		= 6
nstlim 		= ${NSTLIMTI}
numexchg 	= ${NUMEXCHGTI}
dt 		= 0.001
irest 		= 1
ntx 		= 5
ntxo 		= 1
ntc 		= 1
vlimit		= 20
ntf 		= 1
tishake 	= 1 
ntwx 		= 0
ntwr 		= 0
ntpr 		= 1
cut 		= ${CUTOFF}
iwrap 		= 0

ntb 		= 1
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 
ntp 		= 0
klambda 	= 1 

ifsc 		= 1
icfe 		= 1

ifmbar 		= 1
bar_intervall 	= 1
mbar_states 	= \${#lams[@]}
mbar_lambda 	= \${lams[@]}

clambda 	= \${lam}
timask1 	= ${TIMASK1}
timask2 	= ${TIMASK2}
crgmask 	= ""
scmask1 	= ${SCMASK1}
scmask2 	= ${SCMASK2}
scalpha 	= ${SCALPHA}
scbeta 		= ${SCBETA}

gti_cut 	= ${GTICUT}
gti_output 	= 1
gti_add_sc 	= ${GTISC}
gti_scale_beta	= ${GTIBETA}
gti_cut_sc_on 	= ${GTISCON}
gti_cut_sc_off 	= ${GTISCOFF}
gti_lam_sch 	= ${GTILAMSCH}
gti_ele_sc 	= ${GTISCELE}
gti_vdw_sc 	= ${GTISCVDW}
gti_cut_sc 	= ${GTISCCUT}
gti_ele_exp 	= ${GTIEXPELE}
gti_vdw_exp 	= ${GTIEXPVDW}
/
 &ewald
 /

EOF

        cat<<EOF>> inputs/heat.groupfile
-O -p unisc.parm7 -c current/\${mini}.rst7 -i inputs/\${heat}.mdin -o current/\${heat}.mdout -r current/\${heat}.rst7 -x current/\${heat}.nc -inf current/\${heat}.info 
EOF

        cat<<EOF>> inputs/ti.groupfile
-O -p unisc.parm7 -c current/\${heat}.rst7 -i inputs/\${ti}.mdin -o current/\${ti}.mdout -r current/\${ti}.rst7 -x current/\${ti}.nc -inf current/\${ti}.info 
EOF

done

#submit independent jobs
cat<<EOF > equil.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --share
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=${wallclock}

lams=(\${lams[@]})

top=\\\${PWD}
# check if AMBERHOME is set
if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=\\\${AMBERHOME}/bin/pmemd.cuda
for lam in \\\${lams[@]};do
    init=\\\${lam}_init
    mini=\\\${lam}_mini
    heat=\\\${lam}_heat


    echo "running \\\${lam} mini"
    \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${init}.rst7 -i inputs/\\\${mini}.mdin -o current/\\\${mini}.mdout -r current/\\\${mini}.rst7 -x current/\\\${mini}.nc -inf current/\\\${mini}.info -ref current/\\\${init}.rst7 
        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${mini}.rst7
autoimage
trajout current/\\\${mini}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/\\\${mini}_centered.rst7 current/\\\${mini}.rst7
        sleep 1


    echo "running \\\${lam} heat"
    \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${mini}.rst7 -i inputs/\\\${heat}.mdin -o current/\\\${heat}.mdout -r current/\\\${heat}.rst7 -x current/\\\${heat}.nc -inf current/\\\${heat}.info 
        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${heat}.rst7
autoimage
trajout current/\\\${heat}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/\\\${heat}_centered.rst7 current/\\\${heat}.rst7
        sleep 1


done
EOF

# submit group-ed jobs
cat<<EOF > equilgroup.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --time=${wallclock}

top=\\\${PWD}
lams=(\${lams[@]})
steps=(mini heat)
if [ ! -d current ];then mkdir current; fi

for step in \\\${steps[@]}; do
        if [ "\\\$step" == "mini" ]; then
                for lam in \\\${lams[@]};do
                        mini=\\\${lam}_mini
                        init=\\\${lam}_init
			# check if AMBERHOME is set
			if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                        export LAUNCH="srun"
                        export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

                        \\\${LAUNCH} \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${init}.rst7 -i inputs/\\\${mini}.mdin -o current/\\\${mini}.mdout -r current/\\\${mini}.rst7 -x current/\\\${mini}.nc -inf current/\\\${lam}.info -ref current/\\\${init}.rst7 
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${mini}.rst7
autoimage
trajout current/\\\${mini}_centered.rst7
go
quit
EOF2
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${mini}_centered.rst7 current/\\\${mini}.rst7
                done
        else
		# check if AMBERHOME is set
		if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                export LAUNCH="mpirun -np \\\${#lams[@]}"
                export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                export MV2_ENABLE_AFFINITY=0
                \\\${LAUNCH} \\\${EXE} -ng \\\${#lams[@]} -groupfile inputs/current_\\\${step}.groupfile

                for lam in \\\${lams[@]};do
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${step}.rst7
autoimage
trajout current/\\\${lam}_\\\${step}_centered.rst7
go
quit
EOF2
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${lam}_\\\${step}_centered.rst7 current/\\\${lam}_\\\${step}.rst7
                done
        fi

done

EOF


truncate -s0 prod.slurm
cat<<EOF > prod.slurm
#!/bin/bash
#SBATCH --job-name="pr_${stA}~${stB}.slurm"
#SBATCH --output="pr_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="pr_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --time=${wallclock}

lams=(\${lams[@]})
# check if AMBERHOME is set
if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
echo "running replica ti"
mpirun -np \\\${#lams[@]} \\\${EXE} -rem 3 -ng \\\${#lams[@]} -groupfile inputs/ti.groupfile
EOF

if [ "\${REPEX}" == "false" ]; then
        sed -i -e 's/numexchg/!numexchg/g' -e 's/gremd_acyc = 1/!gremd_acyc = 1/g' inputs/\\\${lam}_ti.mdin inputs/\\\${lam}_analyze.mdin
        sed -i 's/ -rem 3//g' inputs/ti.groupfile
fi

EOFN
}


###########################################

#################################
function write_analysis2results {
	local lamlist=("$@")
	cat << EOFN > analysis2results.py
#!/usr/bin/env python2.7
import os
from collections import defaultdict as ddict
prefix="unisc"
${lamlist[@]}
merge_gaps=False

mdin = "TEMPLATE.mdin"
fh = file(mdin,"r")
numexchg=0
nstlim=None
ntwx=None
dt=None
for line in fh:
    cmdstr,sepstr,comstr = line.partition("!")
    if "ntpr" in cmdstr:
        cols = cmdstr.replace("="," ").replace(",","").strip().split()
        for icol in range(len(cols)-1):
            if cols[icol] == "ntpr":
                ntwx = int( cols[icol+1] )
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
        cols = cmdstr.replace("=","").replace(",","").strip().split()
        for icol in range(len(cols)-1):
            if cols[icol] == "nstlim":
                nstlim = int( cols[icol+1] )
                break

if ntwx is None:
    raise Exception("Could not determine ntwx from %s"%(mdin))

if dt is None:
    raise Exception("Could not determine dt from %s"%(mdin))

if nstlim is None:
    raise Exception("Could not determine nstlim from %s"%(mdin))

if numexchg < 1:
    numexchg = 1

dt = dt
nstep_per_sim = nstlim * numexchg
nframe_per_sim = nstep_per_sim / ntwx

if nstep_per_sim % ntwx != 0:
    print "num md steps per simulation is not a multiple of ntwx. Unclear how the simulation time works"

t_per_frame = dt * ntwx
t_per_sim = t_per_frame * nframe_per_sim

dvdl_data = ddict( lambda: ddict( float ) )
efep_data = ddict( lambda: ddict( lambda: ddict( float ) ) )

last_read_sim=0
missing_dirs=[]
for isim in range(1,100001):

    dirstr = "production/%06i"%(isim)
    if not os.path.isdir(dirstr):
        missing_dirs.append(dirstr)
        continue

    if last_read_sim != (isim-1):
        for d in missing_dirs:
            print "TIME GAP! Missing directory: %s"%(d)
        missing_dirs=[]

    t0 = (isim-1) * t_per_sim + t_per_frame

    missing_files=False
    error_msgs=[]
    error=False
    data = ddict(list)
    for lam in lams:
        nframe = 0
        dat = "%s/dvdl_%s.dat"%(dirstr,lam)
        if not os.path.isfile(dat):
            error=True
            missing_files=True
        else:
            fh = file(dat,"r")
            for line in fh:
                cols = line.strip().split()
                if len(cols) == 2:
                    nframe += 1
                    data[lam].append( cols[-1] )
        if nframe != nframe_per_sim and nframe-1 != nframe_per_sim:
            msg="%s expected %i frames, but found %i"%(dat,nframe_per_sim,nframe)
            error_msgs.append(msg)
            error=True
    if not error:
        for iframe in range(nframe_per_sim):
            t = t0 + iframe * t_per_frame
            for lam in lams:
                dvdl_data[t][lam] = data[lam][iframe]

    if missing_files and len(error_msgs) == len(lams):
        print "%s doesn't appear to have been analyzed yet"%(dirstr)
    else:
        for msg in error_msgs:
            print msg
    if len(error_msgs) > 0:
        missing_dirs.append(dirstr)
        continue



    missing_files=False
    error_msgs=[]
    error=False
    data = ddict(lambda: ddict(list))
    for tlam in lams:
        for plam in lams:
            nframe = 0
            dat = "%s/efep_%s_%s.dat"%(dirstr,tlam,plam)
            if not os.path.isfile(dat):
                error=True
                missing_files=True
            else:
                fh = file(dat,"r")
                for line in fh:
                    cols = line.strip().split()
                    if len(cols) == 2:
                        nframe += 1
                        data[tlam][plam].append( cols[-1] )
            if nframe != nframe_per_sim and nframe-1 != nframe_per_sim:
                msg="%s expected %i frames, but found %i"%(dat,nframe_per_sim,nframe)
                error_msgs.append(msg)
                error=True
    if not error:
        for iframe in range(nframe_per_sim):
            t = t0 + iframe * t_per_frame
            for tlam in lams:
                for plam in lams:
                    efep_data[t][tlam][plam] = data[tlam][plam][iframe]


    for msg in error_msgs:
        print msg
    if len(error_msgs) > 0:
        missing_dirs.append(dirstr)
        continue


if not os.path.exists("results/data"):
    os.makedirs("results/data")

ts=[ t for t in dvdl_data ]
if len( ts ) > 0:
    for lam in lams:
        dat = "results/data/dvdl_%s.dat"%(lam)
        fh = file(dat,"w")
        for i,t in enumerate(sorted(dvdl_data)):
            time=t
            if merge_gaps:
                time = (i+1)*t_per_frame
            fh.write("%12.1f %s\n"%(t,dvdl_data[t][lam]))
        fh.close()

ts=[ t for t in dvdl_data ]
if len( ts ) > 0:
    for tlam in lams:
        for plam in lams:
            dat = "results/data/efep_%s_%s.dat"%(tlam,plam)
            fh = file(dat,"w")
            for i,t in enumerate(sorted(efep_data)):
                time=t
                if merge_gaps:
                    time = (i+1)*t_per_frame
                fh.write("%12.1f %s\n"%(t,efep_data[t][tlam][plam]))
            fh.close()
EOFN

}


#################################


#################################
function write_analysis_slurm {
	cat << EOFN > analyze.slurm.sh
#!/bin/bash

main() {

    if [ -z \${AMBERHOME+x} ]; then
       echo "AMBERHOME is unset. Please set AMBERHOME and try again."
       exit 1
    fi

    local base="unisc"
    local lams=( ${lams[@]} )

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

    local last_step=\$(for f in \$(ls production | grep -v '\.' | tail -n 1); do bc -l <<< \$f; done)

    #
    # if we couldn't find a valid subdirectory, then exit now
    #

    if [ "\${last_step}" == "" ]; then
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
    for step in \$(seq \${last_step}); do

       #
       # this is the zero-padded name
       #

       local step_name=\$(printf "%06i" \${step})

       #
       # if the subdir does not exist, then skip it
       #

       if [ ! -d "production/\${step_name}" ]; then
          echo "skipping \${step_name} because dir does not exist"
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
       for lam in \${lams[@]}; do
           if [ -e "production/\${step_name}/\${base}_\${lam}.mdout" ]; then
               cnt=\$(( \${cnt} + 1 ))
           else
               ok=0
           fi
       done
       if [ "\${cnt}" == "0" ]; then
          echo "skipping \${step_name} because it was probably already deleted from this machine"
          continue
       elif [ "\${ok}" == "0" ]; then
          echo "skipping \${step_name} because it has some, but not all, of the mdout files (something wrong here?)"
          continue
       fi


       ok=1
       for lam in \${lams[@]}; do
           if [ ! -e "production/\${step_name}/\${base}_${lam}.nc" ]; then
               ok=0
               echo "Missing production/\${step_name}/\${base}_${lam}.nc"
           fi
       done
       if [ "\${ok}" == "0" ]; then
          echo "skipping \${step_name} because it is missing 1-or-more nc files"
          continue
       fi


       #
       # do we have the mbar trace file generated from the analysis
       #

       local testfile=\$(printf "production/\${step_name}/efep_%.8f_%.8f.dat" 1. 1.)
       if [ -e \${testfile} ]; then
          echo "skipping \${step_name} because it has already been analyzed"
          continue
       fi


       cd production/\${step_name}
       python2.7 ../stdti_step2dats.py \${base}_*.mdout
       cd ../../


    done


    python2.7 analysis2results.py


}



##############################################################################


write_template() {
    local fname="\$1"
    shift
    local sarr=("\$@")
    local nsarr=\${#sarr[@]}
    local lsarr=\$((${nsarr}-1))
    local sline=\${sarr[0]}
    for istep in \$(seq \${lsarr}); do sline="\${sline},\${sarr[${istep}]}"; done

    cat << EOF > production/\${fname}
#!/bin/bash
#SBATCH --job-name="\${fname}"
#SBATCH --output="\${fname}.slurmout"
#SBATCH --error="\${fname}.slurmerr"
#SBATCH --array=\${sline}
EOF
    cat << 'EOF' >> production/\${fname}
#SBATCH --partition=main
#SBATCH --nodes=${nnodes}
#SBATCH --ntasks-per-node=25
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=1-00:00:00
##SBATCH --exclude=cuda[001-008]

export MV2_ENABLE_AFFINITY=0
source \${AMBERHOME}/amber.sh
export LAUNCH="srun --mpi=pmi2"


export EXE="\${AMBERHOME}/bin/pmemd.MPI"

prefix="unisc"
lams=(${lams[@]})
EOF
    cat << 'EOF' >> production/\${fname}
istep=\${SLURM_ARRAY_TASK_ID}


#
# is the input file asking for a replica-exchange simulation?
#
numexchg=\$( grep numexchg ../inputs/\${prefix}_\${lams[0]}_analyze.mdin | sed -e 's/\!.*//' -e 's/ *//g' -e 's/.*numexchg\=\([0-9]*\)/\1/' )
if [ "\${numexchg}" == "" ]; then
   numexchg="0"
fi
rem=0
if [ "\${numexchg}" != "0" ]; then
   rem="3"
fi

#
# if rem > 0, then it IS asking for replica exchange
#

if [ "\${rem}" != "0" ]; then
   echo "CANNOT RUN REPLICA EXCHANGE IN ANALYSIS-MODE"
   exit 1
fi

#
# each step consists of running all windows for some length of time
#

istep=\$(printf "%06i" \${istep})

#
# if don't have the subdirectory, then something's gone wrong
#

if [ ! -d \${istep} ]; then
   exit
fi

cd \${istep}

#
# run analysis for each window
#

for lambda in \${lams[@]}; do
    base="\${prefix}_\${lambda}"
    rest="\${base}_restart"
    init="\${base}_initial"
    anal="\${base}_analyze"

    \${LAUNCH} \${EXE} -O -c \${base}.rst7 -p ../../\${prefix}.parm7 -i ../../inputs/\${anal}.mdin -o \${anal}.mdout -r \${anal}.rst7 -x \${anal}.nc -inf \${anal}.mdinfo -y \${base}.nc

    for tmpfile in \${anal}.rst7 \${anal}.nc \${anal}.mdinfo logfile; do
       if [ -e "\${tmpfile}" ]; then
          rm -f "\${tmpfile}"
       fi
    done
done

python2.7 ../stdti_step2dats.py *_analyze.mdout

ok=1
for lambda in \${lams[@]}; do
    if [ ! -e "dvdl_\${lambda}.dat" ]; then
      ok=0
    fi
done
if [ "\${ok}" == "1" ]; then
    for lambda in \${lams[@]}; do
       base="\${prefix}_\${lambda}"
       anal="\${base}_analyze"
       if [ -e "\${anal}.mdout" ]; then
          rm -f "\${anal}.mdout"
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

EOFN
}
#################################


#################################
function writelatex {
cat << EOFN > makelatex.py
#!/usr/bin/env python2.7

if __name__ == "__main__":

    import matplotlib
    matplotlib.use("Agg")

    from ndmbar.tianalysis import DataLoc as dloc
    from ndmbar.tianalysis import UsePkaUnits
    from ndmbar.tianalysis import make_latex_document
    from collections import defaultdict as ddict
    import os
    import glob

    #UsePkaUnits()

    tequil=0
    tmax=1.e+10
    odir="latex"

    D = ddict( lambda: ddict( lambda: ddict( list ) ) )

    scs=[]
    for itrial in range(10):
        trial="t%i"%(itrial+1)
        d = "com/%s/results/data/"%(trial)
        if os.path.exists(d):
            if len(glob.glob(d+"dvdl_*.dat")) > 0:
                scs.append( dloc(trial,d,"",tequil,tmax) )

    if len(scs) > 0:
        D["lig"]["bio"]["uni"]  = scs

    scs=[]
    for itrial in range(10):
        trial="t%i"%(itrial+1)
        d = "lig/%s/results/data/"%(trial)
        if os.path.exists(d):
            if len(glob.glob(d+"dvdl_*.dat")) > 0:
                scs.append( dloc(trial,d,"",tequil,tmax) )

    if len(scs) > 0:
        D["lig"]["ref"]["uni"]  = scs

    make_latex_document( odir, D, methods=["TI","TI3","BAR","MBAR"] )

EOFN
}

#################################



##################################
# Main program begins

#read and parse input data
parse_input $1
####

#################################
##BEGIN STAGE=setup
if [ "$stage" == "setup" ]; then
#################################

#check if initial structures are present
if [ "${ticalc}" == "rbfe" ]; then dir1=com; dir2=aq; else dir1=aq; dir2=vac; fi

if [ ! -d initial/${system}/${dir1} ] || [ ! -d initial/${system}/${dir2} ]; 
	then echo "initial/${system}/${dir1} or initial/${system}/${dir2} folder(s) missing" && exit 0
else
	for s in ${dir1} ${dir2}; do
		cd initial/${system}/$s
			for i in "${!translist[@]}";do
        			stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
				if [ ! -f ${stA}_${s}.parm7 ] || [ ! -f ${stA}_${s}.rst7 ] || [ ! -f ${stB}_${s}.parm7 ] || [ ! -f ${stB}_${s}.rst7 ]; then
					echo "one or more of ${stA}/${stB} parm/rst files are missing in initial/${system}/${s}"
				fi
			done
		cd $path
	done
fi
#####

#setup TI files for UNIFIED protocol
if [ "${protocol}" == "unified" ]; then
	lams=($(gen_lambdas $nlambda))
	for i in "${!translist[@]}";do
		stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		# generate stateB mol2 file
		cat << EOF > genmol2.in
parm initial/${system}/aq/${stB}_aq.parm7
trajin initial/${system}/aq/${stB}_aq.rst7
strip :WAT,K+,Na+,Cl-
change parmindex 0 resname from :L1 to ${stB}
trajout ${stB}.mol2
go
quit
EOF
		cpptraj < genmol2.in > output
		sleep 1
		for s in ${dir1} ${dir2}; do
			if [ "${mapping}" != "checked" ]; then
				parmutils-timutate.py -p initial/${system}/$s/${stA}_${s}.parm7 -c initial/${system}/$s/${stA}_${s}.rst7 --target ":L1" --mol2 ${stB}.mol2 --uniti --nlambda ${nlambda} >> output 2>&1
				sleep 1
				mkdir -p        ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build
                                mv ticopy.*     ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/
                                cp ${stB}.mol2  ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/

			fi

			if [ "${mapping}" == "manual" ]; then
				cat <<EOF
######################################################################################################################
1. Double check ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/ticopy.map.txt
2. Double check timask1/2 scmask1/2 definitions in ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/TEMPLATE.mdin
######################################################################################################################
EOF
				continue
			else
				if [ ! -d ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build ]; then
					echo "${path}/${system}/${protocol}/${stA}~${stB}/${s}/build does not exist"
			       		continue	       
				else
					cd ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build
						sh ticopy.sh >> output 2>&1
						sleep 1

						timask1=$(grep molmask2 ticopy.sh|awk -F " " '{print $2}'|awk -F "=" '{print $2}')
						scmask1=$(grep molmask2 ticopy.sh|awk -F " " '{print $3}'|awk -F "=" '{print $2}')
						timask2=$(grep molmask2 ticopy.sh|awk -F " " '{print $5}'|awk -F "=" '{print $2}')
						scmask2=$(grep molmask2 ticopy.sh|awk -F " " '{print $6}'|awk -F "=" '{print $2}')
						
						if [ "${ticalc}" == "rbfe" ]; then
							writetemplate_rbfe $nlambda $cutoff $mincyc $nstlimnvt $nstlimnpt $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp $stA $stB
						else
							writetemplate_rsfe $nlambda $cutoff $mincyc $nstlimnvt $nstlimnpt $repex $nstlimti $numexchgti $timask1 $timask2 $scmask1 $scmask2 $scalpha $scbeta $gti_add_sc $gti_scale_beta $gti_cut $gti_cut_sc_on $gti_cut_sc_off $gti_lam_sch $gti_ele_sc $gti_vdw_sc $gti_cut_sc $gti_ele_exp $gti_vdw_exp $stA $stB
						fi
						##########################################
						# setup H-mass repartitioning
						if [ "${hmr}" == "true" ]; then
							if [ -f hmr.parm7 ] || [ -f hmr.rst7 ]; then rm -rf hmr.parm7 hmr.rst7; fi
                					cat <<EOF > hmr.in
HMassRepartition
outparm hmr.parm7 hmr.rst7
EOF
                					parmed -i hmr.in -p unisc/template/unisc.parm7 -c unisc/template/unisc.rst7 >> output 2>&1
                					sleep 1
							mv hmr.parm7 unisc/template/unisc.parm7
							mv hmr.rst7  unisc/template/unisc.rst7
						fi
						##########################################

						sh TEMPLATE.sh
						sleep 1
						# if hmr=true set timestep to 4fs
						if [ "${hmr}" == "true" ]; then
							sed -i '/dt.*.=.*.*/c\dt              = 0.004' ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/inputs/*_ti.mdin
						fi
						mkdir -p ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run
						rm -rf ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run/current ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run/inputs
						cd ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run
							mv -f ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/unisc.parm7	.
							mv -f ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/current 		.
							mv -f ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/inputs 		.
							mv -f ${path}/${system}/${protocol}/${stA}~${stB}/${s}/build/*slurm 		.

							truncate -s0 sub-eq.sh sub-ti.sh
							for(( t=1;t<=${ntrials};t++));do
								mkdir -p t${t}
								cp current/*_init.rst7 t${t}/
								sed "s/current/t${t}/g" inputs/heat.groupfile  > inputs/t${t}_heat.groupfile
                        					sed "s/current/t${t}/g" inputs/npt.groupfile > inputs/t${t}_npt.groupfile
                        					sed "s/current/t${t}/g" inputs/ti.groupfile    > inputs/t${t}_ti.groupfile

								sed 	-e "s/current/t${t}/g" \
									-e "s/EOF2/EOF/g" equilgroup.slurm > t${t}_eq.slurm
                        					sed    	-e "s/current/t${t}/g" \
                               						-e "s/ti.groupfile/t${t}_ti.groupfile/g" prod.slurm > t${t}_ti.slurm

								echo "sbatch t${t}_eq.slurm" >> sub-eq.sh
								echo "sbatch t${t}_ti.slurm" >> sub-ti.sh
							done
						cd $path	
					cd $path
				fi
			fi
		done
		echo "Done with ${stA}~${stB}..."
	done
fi

#################################
##END STAGE=setup
fi
#################################



#################################
##BEGIN STAGE=run-equil
if [ "$stage" == "run-equil" ]; then
#################################
	for i in "${!translist[@]}";do
		stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		for s in com aq; do
			cd ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run
				echo "Submitting TI equil jobs in ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run"
				sh sub-eq.sh
				sleep 1
                        cd ${path}
		done
	done
#################################
##END STAGE=run-equil
fi
#################################



#################################
##BEGIN STAGE=check-equil
if [ "$stage" == "check-equil" ]; then
#################################
	lams=($(gen_lambdas $nlambda))
	for i in "${!translist[@]}";do
		stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		for s in com aq; do
			cd ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run
				for(( t=1;t<=${ntrials};t++));do
					for file in heat npt; do
						for i in "${!lams[@]}";do
							if ! grep -Eq 'EPtot      = **************|NaN' t${t}/${lams[$i]}_${file}.mdout; then
                                                		if grep -Eq 'Final Performance Info' t${t}/${lams[$i]}_${file}.mdout; then continue; fi
                                      	  		fi
                                        		echo "check ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run/t${t}/${lams[$i]}_${file}.mdout"
						done
					done
				done
			cd ${path}
		done
	done
#################################
##END STAGE=check-equil
fi
#################################

#################################
##BEGIN STAGE=run-TI
if [ "$stage" == "run-TI" ]; then
#################################
        for i in "${!translist[@]}";do
                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                for s in com aq; do
                        cd ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run
                                echo "Submitting TI production jobs in ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run"
                                sh sub-ti.sh
                                sleep 1
                        cd ${path}
                done
        done
#################################
##END STAGE=run-TI
fi
#################################


#################################
##BEGIN STAGE=check-TI
if [ "$stage" == "check-TI" ]; then
#################################
        lams=($(gen_lambdas $nlambda))
        for i in "${!translist[@]}";do
                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                for s in com aq; do
			if [ ! -d ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run ]; then echo "Folder ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run missing" && continue; fi
                        cd ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run
                                for(( t=1;t<=${ntrials};t++));do
                                      	 for i in "${!lams[@]}";do
                                                	if ! grep -Eq 'EPtot      = **************|NaN' t${t}/${lams[$i]}_ti.mdout; then
                                               		if grep -Eq 'Final Performance Info' t${t}/${lams[$i]}_ti.mdout; then continue; fi
                                                	fi
                                                	echo "check ${path}/${system}/${protocol}/${stA}~${stB}/${s}/run/t${t}/${lams[$i]}_ti.mdout"
                                        	done
                                done
                        cd ${path}
                done
        done
#################################
##END STAGE=check-TI
fi
#################################


