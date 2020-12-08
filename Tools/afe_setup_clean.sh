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
        varlist=(system translist path_to_input nlambda protocol mapping ntrials cutoff mincyc nstlimnvt nstlimnpt repex nstlimti numexchgti hmr scalpha scbeta gti_add_sc gti_scale_beta gti_cut gti_cut_sc_on gti_cut_sc_off gti_lam_sch gti_ele_sc gti_vdw_sc gti_cut_sc gti_ele_exp gti_vdw_exp stage ticalc partition nnodes ngpus wallclock) 
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
system=Tyk2                     # system
translist=(ejm42~ejm54 ejm42~ejm55 ejm55~ejm54)
path_to_input=../initial	# path to folder containing input configuration files 
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
        if [ ! -d ${path_to_input}/${system}/${dir1} ] || [ ! -d ${path_to_input}/${system}/${dir2} ];
                then echo "${path_to_input}/${system}/${dir1} or ${path_to_input}/${system}/${dir2} folder(s) missing" && exit 0
        else
                for s in ${dir1} ${dir2}; do
                        cd ${path_to_input}/${system}/$s
                                for i in "${!translist[@]}";do
                                        stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                                        if [ ! -f ${stA}_${s}.parm7 ] || [ ! -f ${stA}_${s}.rst7 ] || [ ! -f ${stB}_${s}.parm7 ] || [ ! -f ${stB}_${s}.rst7 ]; then
                                                echo "one or more of ${stA}/${stB} parm/rst files are missing in ${path_to_input}/${system}/${s}"
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
	if [ "${ticalc}" == "rbfe" ]; then dir1=com; dir2=aq; else dir1=aq; dir2=vac; fi



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
rm -rf gen_lambda.py
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


###########################################
# Main program begins

#read and parse input data
parse_input $1
####

#################################
##BEGIN STAGE=setup
if [ "$stage" == "setup" ]; then
#################################


#setup TI files for UNIFIED protocol
if [ "${protocol}" == "unified" ]; then
	lams=($(gen_lambdas $nlambda))
	for i in "${!translist[@]}";do
		stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		# generate stateB mol2 file
		cat << EOF > genmol2.in
parm ${path_to_input}/${system}/aq/${stB}_aq.parm7
trajin ${path_to_input}/${system}/aq/${stB}_aq.rst7
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
				python3 parmutils-timutate.py -p ${path_to_input}/${system}/$s/${stA}_${s}.parm7 -c ${path_to_input}/${system}/$s/${stA}_${s}.rst7 --target ":L1" --mol2 ${stB}.mol2 --uniti --nlambda ${nlambda} >> output 2>&1
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
						sed -i "s,parmutils-tigen.py,python3 ${path}/parmutils-tigen.py,g" ticopy.sh
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
	rm -rf ${path}/*mol2 ${path}/genmol2.in
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
		for s in ${dir1} ${dir2}; do
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
		for s in ${dir1} ${dir2}; do
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
                for s in ${dir1} ${dir2}; do
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
                for s in ${dir1} ${dir2}; do
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


