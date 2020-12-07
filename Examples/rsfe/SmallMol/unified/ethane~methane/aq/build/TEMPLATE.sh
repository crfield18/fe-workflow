#!/bin/bash
top=`pwd`

lams=(0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000)

#generate initial configurations for lambda windows
mkdir -p current inputs
truncate -s0 inputs/heat.groupfile inputs/npt.groupfile inputs/ti.groupfile
for lam in ${lams[@]};do
        cp unisc/template/unisc.rst7  current/${lam}_init.rst7
        cp unisc/template/unisc.parm7 unisc.parm7

        init=${lam}_init
        mini=${lam}_mini
        heat=${lam}_heat
        ti=${lam}_ti

        cat<<EOF>inputs/${lam}_mini.mdin
&cntrl
imin 		= 1
maxcyc 		= 2000
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
cut 		= 8
iwrap 		= 1
ig 		= -1

ntb             = 1 

ifsc 		= 1
icfe 		= 1

clambda 	= ${lam}
timask1 	= "@1-8"
timask2 	= "@9-13"
crgmask 	= ""
scmask1 	= "@6-8"
scmask2 	= ""
scalpha 	= 0.5
scbeta 		= 1.0

gti_cut 	= 1
gti_output 	= 1
gti_add_sc 	= 5
gti_scale_beta	= 1
gti_cut_sc_on 	= 6
gti_cut_sc_off 	= 8
gti_lam_sch 	= 1
gti_ele_sc 	= 1
gti_vdw_sc 	= 1
gti_cut_sc 	= 2
gti_ele_exp 	= 2
gti_vdw_exp 	= 2
/
EOF


        cat<<EOF>inputs/${lam}_heat.mdin
&cntrl
imin 		= 0
nstlim 		= 298000
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
cut 		= 8
iwrap 		= 1

ntb		= 1 
tempi 		= 5.
temp0 		= 298.
ntt 		= 3
gamma_ln 	= 2.
tautp 		= 1 
ntp 		= 0 

clambda 	= ${lam}
ifsc 		= 1
icfe 		= 1

timask1 	= "@1-8"
timask2 	= "@9-13"
crgmask 	= ""
scmask1 	= "@6-8"
scmask2 	= ""
scalpha 	= 0.5
scbeta 		= 1.0

gti_cut 	= 1
gti_output 	= 1
gti_add_sc 	= 5
gti_scale_beta	= 1
gti_cut_sc_on 	= 6
gti_cut_sc_off 	= 8
gti_lam_sch 	= 1
gti_ele_sc 	= 1
gti_vdw_sc 	= 1
gti_cut_sc 	= 2
gti_ele_exp 	= 2
gti_vdw_exp 	= 2
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


        cat<< EOF >inputs/${lam}_ti.mdin
&cntrl
imin 		= 0
nstlim 		= 1000
numexchg 	= 500
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
cut 		= 8
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
mbar_states 	= ${#lams[@]}
mbar_lambda 	= ${lams[@]}

clambda 	= ${lam}
timask1 	= "@1-8"
timask2 	= "@9-13"
crgmask 	= ""
scmask1 	= "@6-8"
scmask2 	= ""
scalpha 	= 0.5
scbeta 		= 1.0

gti_cut 	= 1
gti_output	= 1
gti_add_sc 	= 5
gti_scale_beta	= 1
gti_cut_sc_on 	= 6
gti_cut_sc_off 	= 8
gti_lam_sch 	= 1
gti_ele_sc 	= 1
gti_vdw_sc 	= 1
gti_cut_sc 	= 2
gti_ele_exp 	= 2
gti_vdw_exp 	= 2
gremd_acyc 	= 1 
/
 &ewald
 /
EOF

        cat<< EOF >inputs/${lam}_analyze.mdin
&cntrl
imin 		= 6
nstlim 		= 1000
numexchg 	= 500
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
cut 		= 8
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
mbar_states 	= ${#lams[@]}
mbar_lambda 	= ${lams[@]}

clambda 	= ${lam}
timask1 	= "@1-8"
timask2 	= "@9-13"
crgmask 	= ""
scmask1 	= "@6-8"
scmask2 	= ""
scalpha 	= 0.5
scbeta 		= 1.0

gti_cut 	= 1
gti_output 	= 1
gti_add_sc 	= 5
gti_scale_beta	= 1
gti_cut_sc_on 	= 6
gti_cut_sc_off 	= 8
gti_lam_sch 	= 1
gti_ele_sc 	= 1
gti_vdw_sc 	= 1
gti_cut_sc 	= 2
gti_ele_exp 	= 2
gti_vdw_exp 	= 2
/
 &ewald
 /

EOF

        cat<<EOF>> inputs/heat.groupfile
-O -p unisc.parm7 -c current/${mini}.rst7 -i inputs/${heat}.mdin -o current/${heat}.mdout -r current/${heat}.rst7 -x current/${heat}.nc -inf current/${heat}.info 
EOF

        cat<<EOF>> inputs/ti.groupfile
-O -p unisc.parm7 -c current/${heat}.rst7 -i inputs/${ti}.mdin -o current/${ti}.mdout -r current/${ti}.rst7 -x current/${ti}.nc -inf current/${ti}.info 
EOF

done

#submit independent jobs
cat<<EOF > equil.slurm
#!/bin/bash
#SBATCH --job-name="eq_ethane~methane.slurm"
#SBATCH --output="eq_ethane~methane.slurm.slurmout"
#SBATCH --error="eq_ethane~methane.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${#lams[@]}
#SBATCH --gres=gpu:8
#SBATCH --share
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=24:00:00

lams=(${lams[@]})

top=\${PWD}
# check if AMBERHOME is set
if [ -z "\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=\${AMBERHOME}/bin/pmemd.cuda
for lam in \${lams[@]};do
    init=\${lam}_init
    mini=\${lam}_mini
    heat=\${lam}_heat


    echo "running \${lam} mini"
    \${EXE} -O -p \${top}/unisc.parm7 -c current/\${init}.rst7 -i inputs/\${mini}.mdin -o current/\${mini}.mdout -r current/\${mini}.rst7 -x current/\${mini}.nc -inf current/\${mini}.info -ref current/\${init}.rst7 
        cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin current/\${mini}.rst7
autoimage
trajout current/\${mini}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/\${mini}_centered.rst7 current/\${mini}.rst7
        sleep 1


    echo "running \${lam} heat"
    \${EXE} -O -p \${top}/unisc.parm7 -c current/\${mini}.rst7 -i inputs/\${heat}.mdin -o current/\${heat}.mdout -r current/\${heat}.rst7 -x current/\${heat}.nc -inf current/\${heat}.info 
        cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin current/\${heat}.rst7
autoimage
trajout current/\${heat}_centered.rst7
go
quit
EOF2
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        cpptraj < center.in
        sleep 1
        mv current/\${heat}_centered.rst7 current/\${heat}.rst7
        sleep 1


done
EOF

# submit group-ed jobs
cat<<EOF > equilgroup.slurm
#!/bin/bash
#SBATCH --job-name="eq_ethane~methane.slurm"
#SBATCH --output="eq_ethane~methane.slurm.slurmout"
#SBATCH --error="eq_ethane~methane.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${#lams[@]}
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00

top=\${PWD}
lams=(${lams[@]})
steps=(mini heat)
if [ ! -d current ];then mkdir current; fi

for step in \${steps[@]}; do
        if [ "\$step" == "mini" ]; then
                for lam in \${lams[@]};do
                        mini=\${lam}_mini
                        init=\${lam}_init
			# check if AMBERHOME is set
			if [ -z "\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                        export LAUNCH="srun"
                        export EXE=\${AMBERHOME}/bin/pmemd.cuda

                        \${LAUNCH} \${EXE} -O -p \${top}/unisc.parm7 -c current/\${init}.rst7 -i inputs/\${mini}.mdin -o current/\${mini}.mdout -r current/\${mini}.rst7 -x current/\${mini}.nc -inf current/\${lam}.info -ref current/\${init}.rst7 
                        cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin current/\${mini}.rst7
autoimage
trajout current/\${mini}_centered.rst7
go
quit
EOF2
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\${mini}_centered.rst7 current/\${mini}.rst7
                done
        else
		# check if AMBERHOME is set
		if [ -z "\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

                export LAUNCH="mpirun -np \${#lams[@]}"
                export EXE=\${AMBERHOME}/bin/pmemd.cuda.MPI
                export MV2_ENABLE_AFFINITY=0
                \${LAUNCH} \${EXE} -ng \${#lams[@]} -groupfile inputs/current_\${step}.groupfile

                for lam in \${lams[@]};do
                        cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin current/\${lam}_\${step}.rst7
autoimage
trajout current/\${lam}_\${step}_centered.rst7
go
quit
EOF2
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\${lam}_\${step}_centered.rst7 current/\${lam}_\${step}.rst7
                done
        fi

done

EOF


truncate -s0 prod.slurm
cat<<EOF > prod.slurm
#!/bin/bash
#SBATCH --job-name="pr_ethane~methane.slurm"
#SBATCH --output="pr_ethane~methane.slurm.slurmout"
#SBATCH --error="pr_ethane~methane.slurm.slurmerr"
#SBATCH --partition=null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${#lams[@]}
#SBATCH --gres=gpu:8
#SBATCH --time=24:00:00

lams=(${lams[@]})
# check if AMBERHOME is set
if [ -z "\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

EXE=\${AMBERHOME}/bin/pmemd.cuda.MPI
echo "running replica ti"
mpirun -np \${#lams[@]} \${EXE} -rem 3 -ng \${#lams[@]} -groupfile inputs/ti.groupfile
EOF

if [ "${REPEX}" == "false" ]; then
        sed -i -e 's/numexchg/!numexchg/g' -e 's/gremd_acyc = 1/!gremd_acyc = 1/g' inputs/\${lam}_ti.mdin inputs/\${lam}_analyze.mdin
        sed -i 's/ -rem 3//g' inputs/ti.groupfile
fi

