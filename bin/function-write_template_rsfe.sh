#write TEMPLATE
#for relative solution free energies -> [ lig in water - lig in vacuum ]
#################################
function writetemplate_rsfe
{
        #echo "$#"
        local args=$*; args=($args)
        local varlist=(nlambda CUTOFF REPEX NSTLIMTI NUMEXCHGTI TIMASK1 TIMASK2 SCMASK1 SCMASK2 NOSHAKEMASK SCALPHA SCBETA GTISC GTIBETA GTICUT GTISCON GTISCOFF GTILAMSCH GTISCELE GTISCVDW GTISCCUT GTIEXPELE GTIEXPVDW stA stB twostate)
        local i=0
        for a in "${varlist[@]}"; do
                declare -n arr="$a"
                arr=${args[$i]}
                i=$(($i+1))
        done

	lams=($(gen_lambdas $nlambda))

	eqstagelist=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqV eqP eqA minTI eqpre1P0TI eqpre2P0TI eqP0TI eqATI preTI)
	endstates=(0.00000000)
	parmbase=unisc
	rstbase=(stateA stateA)


	cat <<EOFN >TEMPLATE.sh
#!/bin/bash
top=\`pwd\`

parmbase=${parmbase}
rstbase=(${rstbase[*]})
endstates=(${endstates[*]})
twostate=${twostate}
lams=(${lams[@]})

#generate initial configurations for lambda windows
mkdir -p current inputs
count=0
for lam in \${endstates[@]};do
        cp \${rstbase[\$count]}.rst7  current/\${lam}_init.rst7
        #cp \${parmbase}.parm7 unisc.parm7

        init=\${lam}_init
        min1=\${lam}_min1
        min2=\${lam}_min2
	eqpre1P0=\${lam}_eqpre1P0
	eqpre2P0=\${lam}_eqpre2P0
        eqP0=\${lam}_eqP0
        eqNTP4=\${lam}_eqNTP4 
        eqV=\${lam}_eqV
        eqP=\${lam}_eqP
        eqA=\${lam}_eqA
        ti=\${lam}_ti

        cat<<EOF>inputs/\${min1}.mdin
Minimization with Cartesian restraints for the solute
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5

ntr             = 1 
restraintmask 	= '!:WAT,Cl-,K+,Na+ & !@H='
restraint_wt	= 5 

ifsc            = 1
icfe            = 1

noshakemask	= ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}

gti_syn_mass	= 0
/
EOF
	cat<<EOF>inputs/\${min2}.mdin
Minimization of the entire molecular system
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5


ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}

gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/\${eqpre1P0}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ntr             = 1
restraint_wt    = 5
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}
clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/\${eqpre2P0}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ntr             = 1
restraint_wt    = 5
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/\${eqP0}.mdin
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1 
ntf             = 1 
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp		= 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2     
barostat	= 2

ig		= -1

ntr		= 1
restraint_wt	= 5
restraintmask 	= '!:WAT,Cl-,K+,Na+ & !@H='

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass	= 0
/
EOF
        cat<<EOF>inputs/\${eqNTP4}.mdin.template
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 4
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0

&end
 &ewald
   target_a=ABOX,
   target_b=BBOX,
   target_c=CBOX,
   target_n=500,

/
EOF
        cat<<EOF>inputs/\${eqV}.mdin
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 1
ntp		= 0
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2

ig              = -1

ntr             = 1
restraint_wt    = 5
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF

        cat<<EOF>inputs/\${eqP}.mdin
&cntrl
imin            = 0
nstlim          = 500000 
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 1 
ntf             = 1 
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2 
ntp		= 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat	= 2

ig              = -1

ntr             = 1
restraint_wt    = 5
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF

        cat<<EOF>inputs/\${eqA}.mdin
&cntrl
imin            = 0
nstlim          = 2000000
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 1 
ntf             = 1 
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 1
tempi           = 300.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2     

ig              = -1

ntr             = 1
restraint_wt    = 5
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='

nmropt		= 1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
&wt type='TEMP0', istep1=0,istep2=50000,value1=300.,
           value2=600.,    /
&wt type='TEMP0', istep1=50001, istep2=150000, value1=600.0,
           value2=600.0,     /
&wt type='TEMP0', istep1=150001, istep2=200000, value1=600.0,
           value2=300.0,     /
&wt type='END'  /

EOF

count=\$((\$count+1))
done

if [ \${twostate} == true ]; then

	truncate -s0 inputs/eqpre1P0.groupfile inputs/eqpre2P0.groupfile inputs/eqP0.groupfile inputs/eqNTP4.groupfile inputs/eqV.groupfile inputs/eqP.groupfile inputs/eqA.groupfile inputs/eqProt2.groupfile inputs/eqProt1.groupfile inputs/eqProt05.groupfile inputs/eqProt025.groupfile inputs/eqProt01.groupfile inputs/eqProt0.groupfile
	for lam in \${endstates[@]};do

        	init=\${lam}_init
        	min1=\${lam}_min1
        	min2=\${lam}_min2
		eqpre1P0=\${lam}_eqpre1P0
		eqpre2P0=\${lam}_eqpre2P0
        	eqP0=\${lam}_eqP0
        	eqNTP4=\${lam}_eqNTP4
        	eqV=\${lam}_eqV
        	eqP=\${lam}_eqP
        	eqA=\${lam}_eqA

                cat<<EOF>> inputs/eqpre1P0.groupfile
-O -p unisc.parm7 -c current/\${min2}.rst7 -i inputs/\${eqpre1P0}.mdin -o current/\${eqpre1P0}.mdout -r current/\${eqpre1P0}.rst7 -x current/\${eqpre1P0}.nc -ref current/\${min2}.rst7
EOF
                cat<<EOF>> inputs/eqpre2P0.groupfile
-O -p unisc.parm7 -c current/\${eqpre1P0}.rst7 -i inputs/\${eqpre2P0}.mdin -o current/\${eqpre2P0}.mdout -r current/\${eqpre2P0}.rst7 -x current/\${eqpre2P0}.nc -ref current/\${eqpre1P0}.rst7
EOF
        	cat<<EOF>> inputs/eqP0.groupfile
-O -p unisc.parm7 -c current/\${min2}.rst7 -i inputs/\${eqP0}.mdin -o current/\${eqP0}.mdout -r current/\${eqP0}.rst7 -x current/\${eqP0}.nc -ref current/\${min2}.rst7
EOF
        	cat<<EOF>> inputs/eqNTP4.groupfile
-O -p unisc.parm7 -c current/\${eqP0}.rst7 -i inputs/\${eqNTP4}.mdin -o current/\${eqNTP4}.mdout -r current/\${eqNTP4}.rst7 -x current/\${eqNTP4}.nc -ref current/\${eqP0}.rst7
EOF
        	cat<<EOF>> inputs/eqV.groupfile
-O -p unisc.parm7 -c current/\${eqNTP4}.rst7 -i inputs/\${eqV}.mdin -o current/\${eqV}.mdout -r current/\${eqV}.rst7 -x current/\${eqV}.nc -ref current/\${eqNTP4}.rst7
EOF
        	cat<<EOF>> inputs/eqP.groupfile
-O -p unisc.parm7 -c current/\${eqV}.rst7 -i inputs/\${eqP}.mdin -o current/\${eqP}.mdout -r current/\${eqP}.rst7 -x current/\${eqP}.nc -ref current/\${eqV}.rst7
EOF
        	cat<<EOF>> inputs/eqA.groupfile
-O -p unisc.parm7 -c current/\${eqP}.rst7 -i inputs/\${eqA}.mdin -o current/\${eqA}.mdout -r current/\${eqA}.rst7 -x current/\${eqA}.nc -ref current/\${eqP}.rst7
EOF

	done
fi


for lam in \${lams[@]}; do

        init2=\${lam}_init2
        minTI=\${lam}_minTI
        eqpre1P0TI=\${lam}_eqpre1P0TI
        eqpre2P0TI=\${lam}_eqpre2P0TI
        eqP0TI=\${lam}_eqP0TI
        eqATI=\${lam}_eqATI
        preTI=\${lam}_preTI
        ti=\${lam}_ti
        anal=\${lam}_analyze

	cat << EOF > inputs/\${minTI}.mdin
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5


ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0

/
EOF
        cat<<EOF>inputs/\${eqpre1P0TI}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/\${eqpre2P0TI}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/\${eqP0TI}.mdin
&cntrl
imin            = 0
nstlim          = 100000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF
        cat << EOF > inputs/\${eqATI}.mdin
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
&wt
 TYPE='TEMP0',
 ISTEP1=0, ISTEP2=298000,
 VALUE1=0.0, VALUE2=298.0,
&wt
 TYPE='TEMP0',
 ISTEP1=298001, ISTEP2=500000,
 VALUE1=298.0, VALUE2=298.0,
/
&wt TYPE='END' /
EOF
	cat << EOF > inputs/\${preTI}.mdin
&cntrl
imin            = 0
nstlim          = 2000000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ifsc            = 1
icfe            = 1

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
gti_syn_mass    = 0
/
EOF
        cat<< EOF >inputs/\${ti}.mdin
&cntrl
imin            = 0
nstlim          = ${NSTLIMTI}
numexchg        = ${NUMEXCHGTI}
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = ${CUTOFF}
iwrap           = 1

ntb             = 2
temp0           = 298.
ntt             = 3
gamma_ln        = 2.
tautp           = 1
ntp             = 1
barostat        = 2
pres0           = 1.01325
taup            = 5.0

ifsc            = 1
icfe            = 1

ifmbar          = 1
bar_intervall   = 1
mbar_states     = \${#lams[@]}
EOF
ilam=0
for l in \${lams[@]}; do
     ilam=\$(( \${ilam} + 1 ))
     echo "mbar_lambda(\${ilam})     = \${l}" >> inputs/\${ti}.mdin
done
        cat<< EOF >> inputs/\${ti}.mdin

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}
timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}

gremd_acyc      = 1
/
 &ewald
 /
EOF

        cat<< EOF >inputs/\${lam}_analyze.mdin
&cntrl
imin            = 6
nstlim          = ${NSTLIMTI}
numexchg        = ${NUMEXCHGTI}
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 1
ntf             = 1
ntwx            = 0
ntwr            = 0
ntpr            = 1
cut             = ${CUTOFF}
iwrap           = 0

ntb             = 2
temp0           = 298.
ntt             = 3
gamma_ln        = 2.
tautp           = 1
ntp             = 1
barostat        = 2
pres0           = 1.01325
taup            = 5.0

ifsc            = 1
icfe            = 1

ifmbar          = 1
bar_intervall   = 1
mbar_states     = \${#lams[@]}
EOF
ilam=0
for l in \${lams[@]}; do
     ilam=\$(( \${ilam} + 1 ))
     echo "mbar_lambda(\${ilam})     = \${l}" >> inputs/\${ti}.mdin
done
        cat<< EOF >> inputs/\${ti}.mdin

noshakemask     = ${NOSHAKEMASK}

clambda         = \${lam}

timask1         = ${TIMASK1}
timask2         = ${TIMASK2}
crgmask         = ""
scmask1         = ${SCMASK1}
scmask2         = ${SCMASK2}
scalpha         = ${SCALPHA}
scbeta          = ${SCBETA}

gti_cut         = ${GTICUT}
gti_output      = 1
gti_add_sc      = ${GTISC}
gti_scale_beta  = ${GTIBETA}
gti_cut_sc_on   = ${GTISCON}
gti_cut_sc_off  = ${GTISCOFF}
gti_lam_sch     = ${GTILAMSCH}
gti_ele_sc      = ${GTISCELE}
gti_vdw_sc      = ${GTISCVDW}
gti_cut_sc      = ${GTISCCUT}
gti_ele_exp     = ${GTIEXPELE}
gti_vdw_exp     = ${GTIEXPVDW}
/
 &ewald
 /

EOF
done

truncate -s0 inputs/eqATI.groupfile inputs/preTI.groupfile inputs/ti.groupfile
for lam in \${lams[@]};do
	minTI=\${lam}_minTI
        eqpre1P0TI=\${lam}_eqpre1P0TI
        eqpre2P0TI=\${lam}_eqpre2P0TI
        eqP0TI=\${lam}_eqP0TI
	eqATI=\${lam}_eqATI
	preTI=\${lam}_preTI
	ti=\${lam}_ti

	cat<<EOF>> inputs/eqATI.groupfile
-O -p unisc.parm7 -c current/\${eqP0TI}.rst7 -i inputs/\${eqATI}.mdin -o current/\${eqATI}.mdout -r current/\${eqATI}.rst7 -x current/\${eqATI}.nc -ref current/\${minTI}.rst7
EOF
	cat<<EOF>> inputs/preTI.groupfile
-O -p unisc.parm7 -c current/\${eqATI}.rst7 -i inputs/\${preTI}.mdin -o current/\${preTI}.mdout -r current/\${preTI}.rst7 -x current/\${preTI}.nc -ref current/\${eqATI}.rst7
EOF
	cat<<EOF>> inputs/ti.groupfile
-O -p unisc.parm7 -c current/\${preTI}.rst7 -i inputs/\${ti}.mdin -o current/\${ti}.mdout -r current/\${ti}.rst7 -x current/\${ti}.nc -ref current/\${preTI}.rst7
EOF

done

# slurm to submit all trials
###############
###############

cat<<EOF > run_alltrials.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${ntrials}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --share
#SBATCH --gres-flags=enforce-binding
#SBATCH --export=ALL
#SBATCH --time=${wallclock}

endstates=(\${endstates[@]})
lams=(\${lams[@]})
twostate=\${twostate}
eqstage=(${eqstagelist[*]})

top=\\\${PWD}
# check if AMBERHOME is set
if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi
# check if cpptraj is present
if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi

EXE=\\\${AMBERHOME}/bin/pmemd.cuda

for trial in \\\$(seq 1 1 ${ntrials}); do

	count=-1; alllams=0
	for stage in \\\${eqstage[@]}; do
        	count=\\\$((\\\$count+1))
        	lastcount=\\\$((\\\$count-1))
		if [ "\\\${stage}" == "init" ] || [ "\\\${stage}" == "eqpre1P0TI" ] || [ "\\\${stage}" == "eqpre2P0TI" ] || [ "\\\${stage}" == "eqP0TI" ]; then continue; fi
        	laststage=\\\${eqstage[\\\${lastcount}]}


        	if [ "\\\${stage}" == "minTI" ];then
                	for lam in \\\${lams[@]}; do
                        	cp t\\\${trial}/\\\${endstates[0]}_\\\${laststage}.rst7 t\\\${trial}/\\\${lam}_init2.rst7
                	done
                	laststage="init2"
                	alllams=1
        	fi


        	if [ \\\${alllams} -eq 0 ];then
                	lam=\${endstates[0]}
                	echo "Running \\\$stage for lambda \\\${lam}..."
                	\\\${EXE} -O -p \\\${top}/unisc.parm7 -c t\\\${trial}/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o t\\\${trial}/\\\${lam}_\\\${stage}.mdout -r t\\\${trial}/\\\${lam}_\\\${stage}.rst7 -ref t\\\${trial}/\\\${lam}_\\\${laststage}.rst7
                	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin t\\\${trial}/\\\${lam}_\\\${stage}.rst7
autoimage
trajout t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                	# check if cpptraj is present
                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                	cpptraj < center.in
                	sleep 1
                	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7

		else
			if [ "\\\${stage}" == "minTI" ];then laststage="init2"; fi
			# check if pmemd.cuda is present
			if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
			export LAUNCH="srun"
			export EXE=\\\${AMBERHOME}/bin/pmemd.cuda
			for lam in \\\${lams[@]}; do
				echo "Running \\\$stage for lambda \\\${lam}..."
                                \\\${EXE} -O -p \\\${top}/unisc.parm7 -c t\\\${trial}/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o t\\\${trial}/\\\${lam}_\\\${stage}.mdout -r t\\\${trial}/\\\${lam}_\\\${stage}.rst7 -ref t\\\${trial}/\\\${lam}_\\\${laststage}.rst7
                                sleep 1
				cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin t\\\${trial}/\\\${lam}_\\\${stage}.rst7
autoimage
trajout t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                cpptraj < center.in
                                sleep 1
                                mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                        done

        	fi

	done

	cat << EOFP > extract.py
#!/usr/bin/env python2.7

def OpenParm( fname, xyz=None ):
    import parmed
    from parmed.constants import IFBOX
    if ".mol2" in fname:
        param = parmed.load_file( fname, structure=True )
        #help(param)
    else:
        param = parmed.load_file( fname,xyz=xyz )
        if xyz is not None:
            if ".rst7" in xyz:
                param.load_rst7(xyz)
    if param.box is not None:
        if abs(param.box[3]-109.471219)<1.e-4 and \
           abs(param.box[4]-109.471219)<1.e-4 and \
           abs(param.box[5]-109.471219)<1.e-4:
            param.parm_data["POINTERS"][IFBOX]=2
            param.pointers["IFBOX"]=2
    return param

def CopyParm( parm ):
    import copy
    try:
        parm.remake_parm()
    except:
        pass
    p = copy.copy( parm )
    p.coordinates = copy.copy( parm.coordinates )
    p.box = copy.copy( parm.box )
    try:
        p.hasbox = copy.copy( parm.hasbox )
    except:
        p.hasbox = False
    return p

def Strip( parm, mask ):
    p = CopyParm( parm )
    p.strip( "%s"%(mask) )
    return p

def Extract( parm, mask ):
    return Strip( parm, "!(%s)"%(mask) )

def SaveParmRst( param, fname ):
    from parmed.constants import IFBOX
    for a in param.atoms:
        param.parm_data["CHARGE"][ a.idx ] = a.charge
    if param.box is not None:
       if abs(param.box[3]-109.471219)<1.e-4 and \
          abs(param.box[4]-109.471219)<1.e-4 and \
          abs(param.box[5]-109.471219)<1.e-4:
           param.parm_data["POINTERS"][IFBOX]=2
           param.pointers["IFBOX"]=2
    try:
        param.save( "{}.parm7".format(fname), overwrite=True )
        #param.save( fname, overwrite=True )
    except:
        param.save( "{}.parm7".format(fname) )
        #param.save( fname )
    rst = parmed.amber.Rst7(natom=len(param.atoms),title="BLAH")
    rst.coordinates = param.coordinates
    rst.box = [param.box[0], param.box[1], param.box[2], param.box[3], param.box[4], param.box[5]]
    rst.write( "{}.rst7".format(fname) )

if __name__ == "__main__":

    import argparse
    import parmed
    import re
    import sys

    parser = argparse.ArgumentParser \\
    ( formatter_class=argparse.RawDescriptionHelpFormatter,
      description="Extracts Amber parameters from a parm7 file" )

    parser.add_argument("-p","--parm",
                        help="Amber parm7 file",
                        type=str,
                        required=True)

    parser.add_argument("-c","--crd",
                        help="Amber rst7 file",
                        type=str,
                        required=True)

    parser.add_argument("-m","--mask",
                        help="Amber selection mask. Only residues within this mask will be extracted",
                        type=str,
                        default="@*",
                        required=False )

    parser.add_argument("-o","--output",
                        help="Output file basename if -n/--name is specified",
                        type=str,
                        required=False)

    args = parser.parse_args()

    p    = OpenParm(args.parm,xyz=args.crd)

    if args.mask is not None:
        if len(args.mask) > 0:
            q    = Extract(p,args.mask)
            SaveParmRst(q, args.output)

EOFP
	chmod a+x extract.py

	mkdir -p ../vac/t\\\${trial} ../vac/inputs
	for lam in \\\${lams[@]};do
        	./extract.py -p unisc.parm7 -c t\\\${trial}/\\\${lam}_preTI.rst7 -m '!:WAT,Na+,K+,Cl-' -o ../vac/t\\\${trial}/\\\${lam}_init
        	sed -e 's/nstlim.*/nstlim          = 500000/g' -e 's/restraint_wt.*/restraint_wt    = 0/g' -e 's/nmropt.*/nmropt          = 0/g' -e '59,65d' -e "s/clambda.*/clambda         = ${lam}/g" -e 's/irest.*/irest           = 0/g' -e 's/ntx.*/ntx             = 1/g' inputs/0.00000000_eqA.mdin > ../vac/inputs/\\\${lam}_preTI.mdin
        	sed -e 's/ntb.*/ntb             = 1/g' -e '/barostat.*/d' -e '/ntp.*/d' -e '/pres0.*/d' -e '/taup.*/d' -e '/numexchg/d' -e '/gremd_acyc/d' -e 's/nstlim.*/nstlim          = 1000000/g' inputs/\\\${lam}_ti.mdin > ../vac/inputs/\\\${lam}_ti.mdin

	done

	cd ../vac
        	cd t\\\${trial}
                	mv 0.00000000_init.parm7 ../unisc.parm7
                	rm *.parm7
        	cd ../

        	stage=preTI; laststage=init
        	EXE=\\\${AMBERHOME}/bin/pmemd.cuda
        	LAUNCH="srun --exclusive -N 1 -n 1 -c 1 --gres=gpu:1"
        	for lam in \\\${lams[@]};do
                	\\\${LAUNCH} \\\${EXE} -O -p unisc.parm7 -c t\\\${trial}/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o t\\\${trial}/\\\${lam}_\\\${stage}.mdout -r t\\\${trial}/\\\${lam}_\\\${stage}.rst7 -ref t\\\${trial}/\\\${lam}_\\\${laststage}.rst7 &
        	done
        	wait
	cd ../

	for dir in aq vac; do
		cd \\\${dir}
			stage=ti; laststage=preTI
			for lam in \\\${lams[@]};do
				\\\${LAUNCH} \\\${EXE} -O -p unisc.parm7 -c t\\\${trial}/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o t\\\${trial}/\\\${lam}_\\\${stage}.mdout -r t\\\${trial}/\\\${lam}_\\\${stage}.rst7 -ref t\\\${trial}/\\\${lam}_\\\${laststage}.rst7 &
			done
		cd ../
	done
	wait

done

EOF

EOFN

}

#################################

