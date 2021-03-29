#write TEMPLATE
#for relative binding free energies -> [ protein-lig complex - lig in water ]
#################################
function writetemplate_rbfe
{
        #echo "$#"
        local args=$*; args=($args)
        local varlist=(nlambda CUTOFF REPEX NSTLIMTI NUMEXCHGTI TIMASK1 TIMASK2 SCMASK1 SCMASK2 NOSHAKEMASK SCALPHA SCBETA GTISC GTIBETA GTICUT GTISCON GTISCOFF GTILAMSCH GTISCELE GTISCVDW GTISCCUT GTIEXPELE GTIEXPVDW stA stB s twostate)
        local i=0
        for a in "${varlist[@]}"; do
                declare -n arr="$a"
                arr=${args[$i]}
                i=$(($i+1))
        done

	lams=($(gen_lambdas $nlambda))
	if [ "${twostate}" == "true" ]; then
		endstates=(0.00000000 1.00000000)
		parmbase=merged-lam0
		rstbase=(merged-lam0 merged-lam1)
		if [ "$s" == "com" ]; then
	                eqstagelist=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqNTP4 eqV eqP eqA eqProt2 eqProt1 eqProt05 eqProt025 eqProt01 eqProt0 minTI eqATI preTI)
        	else
                	eqstagelist=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqNTP4 eqV eqP eqA minTI eqATI preTI)
        	fi

	else
		endstates=(0.00000000)
		parmbase=unisc
		rstbase=(unisc unisc)
		if [ "$s" == "com" ]; then
    	        	eqstagelist=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqV eqP eqA eqProt2 eqProt1 eqProt05 eqProt025 eqProt01 eqProt0 minTI eqATI preTI)
        	else
                	eqstagelist=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqV eqP eqA minTI eqATI preTI)
        	fi

	fi
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
        cp unisc/template/\${rstbase[\$count]}.rst7  current/\${lam}_init.rst7
        cp unisc/template/\${parmbase}.parm7 unisc.parm7

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
        eqProt2=\${lam}_eqProt2
        eqProt1=\${lam}_eqProt1
        eqProt05=\${lam}_eqProt05
        eqProt025=\${lam}_eqProt025
        eqProt01=\${lam}_eqProt01
        eqProt0=\${lam}_eqProt0
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

        cat<<EOF>inputs/\${eqProt2}.mdin
&cntrl
imin            = 0
nstlim          = 200000
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
ntp             = 1
tempi           = 5.
temp0           = 300.
ntt             = 3
gamma_ln        = 2.
tautp           = 2
barostat        = 2

ig              = -1

ntr             = 1
restraint_wt    = 2
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
        cat<<EOF>inputs/\${minProt2}.mdin
Minimization with Cartesian restraints for the solute
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5

ntr             = 1
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='
restraint_wt    = 2

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
        cat<<EOF>inputs/\${eqProt1}.mdin
&cntrl
imin            = 0
nstlim          = 200000
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
restraint_wt    = 1 
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
        cat<<EOF>inputs/\${minProt1}.mdin
Minimization with Cartesian restraints for the solute
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5

ntr             = 1
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='
restraint_wt    = 1

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
        cat<<EOF>inputs/\${eqProt05}.mdin
&cntrl
imin            = 0
nstlim          = 200000
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
restraint_wt    = 0.5
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
        cat<<EOF>inputs/\${minProt05}.mdin
Minimization with Cartesian restraints for the solute
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5

ntr             = 1
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='
restraint_wt    = 0.5

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
        cat<<EOF>inputs/\${eqProt025}.mdin
&cntrl
imin            = 0
nstlim          = 200000 
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
restraint_wt    = 0.25
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
        cat<<EOF>inputs/\${minProt025}.mdin
Minimization with Cartesian restraints for the solute
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5

ntr             = 1
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='
restraint_wt    = 0.25

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
        cat<<EOF>inputs/\${eqProt01}.mdin
&cntrl
imin            = 0
nstlim          = 200000 
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
restraint_wt    = 0.1
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
        cat<<EOF>inputs/\${minProt01}.mdin
Minimization with Cartesian restraints for the solute
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5

ntr             = 1
restraintmask   = '!:WAT,Cl-,K+,Na+ & !@H='
restraint_wt    = 0.1

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
        cat<<EOF>inputs/\${eqProt0}.mdin
&cntrl
imin            = 0
nstlim          = 200000 
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
        	eqProt2=\${lam}_eqProt2
       	 	eqProt1=\${lam}_eqProt1
        	eqProt05=\${lam}_eqProt05
        	eqProt025=\${lam}_eqProt025
        	eqProt01=\${lam}_eqProt01
        	eqProt0=\${lam}_eqProt0

        	cat<<EOF>> inputs/eqpre1P0.groupfile
-O -p unisc.parm7 -c current/\${min2}.rst7 -i inputs/\${eqpre1P0}.mdin -o current/\${eqpre1P0}.mdout -r current/\${eqpre1P0}.rst7 -x current/\${eqpre1P0}.nc -ref current/\${min2}.rst7
EOF
        	cat<<EOF>> inputs/eqpre2P0.groupfile
-O -p unisc.parm7 -c current/\${eqpre1P0}.rst7 -i inputs/\${eqpre2P0}.mdin -o current/\${eqpre2P0}.mdout -r current/\${eqpre2P0}.rst7 -x current/\${eqpre2P0}.nc -ref current/\${eqpre1P0}.rst7
EOF
        	cat<<EOF>> inputs/eqP0.groupfile
-O -p unisc.parm7 -c current/\${eqpre2P0}.rst7 -i inputs/\${eqP0}.mdin -o current/\${eqP0}.mdout -r current/\${eqP0}.rst7 -x current/\${eqP0}.nc -ref current/\${eqpre2P0}.rst7
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
        	cat<<EOF>> inputs/eqProt2.groupfile
-O -p unisc.parm7 -c current/\${eqA}.rst7 -i inputs/\${eqProt2}.mdin -o current/\${eqProt2}.mdout -r current/\${eqProt2}.rst7 -x current/\${eqProt2}.nc -ref current/\${eqA}.rst7
EOF
        	cat<<EOF>> inputs/eqProt1.groupfile
-O -p unisc.parm7 -c current/\${eqProt2}.rst7 -i inputs/\${eqProt1}.mdin -o current/\${eqProt1}.mdout -r current/\${eqProt1}.rst7 -x current/\${eqProt1}.nc -ref current/\${eqProt2}.rst7
EOF
        	cat<<EOF>> inputs/eqProt05.groupfile
-O -p unisc.parm7 -c current/\${eqProt1}.rst7 -i inputs/\${eqProt05}.mdin -o current/\${eqProt05}.mdout -r current/\${eqProt05}.rst7 -x current/\${eqProt05}.nc -ref current/\${eqProt1}.rst7
EOF
        	cat<<EOF>> inputs/eqProt025.groupfile
-O -p unisc.parm7 -c current/\${eqProt05}.rst7 -i inputs/\${eqProt025}.mdin -o current/\${eqProt025}.mdout -r current/\${eqProt025}.rst7 -x current/\${eqProt025}.nc -ref current/\${eqProt05}.rst7
EOF
        	cat<<EOF>> inputs/eqProt01.groupfile
-O -p unisc.parm7 -c current/\${eqProt025}.rst7 -i inputs/\${eqProt01}.mdin -o current/\${eqProt01}.mdout -r current/\${eqProt01}.rst7 -x current/\${eqProt01}.nc -ref current/\${eqProt025}.rst7
EOF
        	cat<<EOF>> inputs/eqProt0.groupfile
-O -p unisc.parm7 -c current/\${eqProt01}.rst7 -i inputs/\${eqProt0}.mdin -o current/\${eqProt0}.mdout -r current/\${eqProt0}.rst7 -x current/\${eqProt0}.nc -ref current/\${eqProt01}.rst7
EOF

	done
fi


for lam in \${lams[@]}; do

        init2=\${lam}_init2
        minTI=\${lam}_minTI
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
mbar_lambda     = \${lams[@]}

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
mbar_lambda     = \${lams[@]}

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
	eqATI=\${lam}_eqATI
	preTI=\${lam}_preTI
	ti=\${lam}_ti

	cat<<EOF>> inputs/eqATI.groupfile
-O -p unisc.parm7 -c current/\${minTI}.rst7 -i inputs/\${eqATI}.mdin -o current/\${eqATI}.mdout -r current/\${eqATI}.rst7 -x current/\${eqATI}.nc -ref current/\${minTI}.rst7
EOF
	cat<<EOF>> inputs/preTI.groupfile
-O -p unisc.parm7 -c current/\${eqATI}.rst7 -i inputs/\${preTI}.mdin -o current/\${preTI}.mdout -r current/\${preTI}.rst7 -x current/\${preTI}.nc -ref current/\${eqATI}.rst7
EOF
	cat<<EOF>> inputs/ti.groupfile
-O -p unisc.parm7 -c current/\${preTI}.rst7 -i inputs/\${ti}.mdin -o current/\${ti}.mdout -r current/\${ti}.rst7 -x current/\${ti}.nc -ref current/\${preTI}.rst7
EOF

done

# slurm to submit all trials together
#############
#############

if [ "\${twostate}" != "true" ]; then
        cat<<EOF > equil_alltrials.slurm
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

	if [ ! -d "t\\\${trial}" ]; then mkdir -p t\\\${trial}; fi

	count=-1; alllams=0
	for stage in \\\${eqstage[@]}; do
        	count=\\\$((\\\$count+1))
        	lastcount=\\\$((\\\$count-1))
        	if [ "\\\${stage}" == init ]; then continue; fi
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

                	for lam in \\\${lams[@]};do
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
                	done

        	fi

	done

done
EOF

else
        cat<<EOF > equil_alltrials.slurm
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

	if [ ! -d "t\\\${trial}" ]; then mkdir -p t\\\${trial}; fi

	count=-1; alllams=0
	for stage in \\\${eqstage[@]}; do
        	count=\\\$((\\\$count+1))
        	lastcount=\\\$((\\\$count-1))
        	if [ "\\\${stage}" == init ]; then continue; fi
        	laststage=\\\${eqstage[\\\${lastcount}]}

        	if [ "\\\${stage}" == "minTI" ];then
	
        	        for lam in \\\${lams[@]}; do
                	        if [ \\\$(bc -l <<< "\\\${lam} <= 0.5") -eq 1 ]; then
                        	        cp t\\\${trial}/\\\${endstates[0]}_\\\${laststage}.rst7 t\\\${trial}/\\\${lam}_init2.rst7
                        	else
                                	cp t\\\${trial}/\\\${endstates[1]}_\\\${laststage}.rst7 t\\\${trial}/\\\${lam}_init2.rst7
                        	fi
                	done
                	laststage="init2"
                	alllams=1

        	fi

        	if [ \\\${alllams} -eq 0 ];then

                	for lam in \\\${endstates[@]};do
                        	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                       		# check if pmemd.cuda is present
                        	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
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
                	done

        	else

                	for lam in \\\${lams[@]};do
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
                	done

        	fi

        	if [ "\\\${stage}" == "eqP0" ]; then
                	for lam in \\\${endstates[@]};do
                        	box=(\\\$(tail -1 \\\${lam}_\\\${laststage}.rst7))
                        	a=\\\$(awk "BEGIN {print ( \\\$a + \\\${box[0]} ) }")
                        	b=\\\$(awk "BEGIN {print ( \\\$b + \\\${box[1]} ) }")
                        	c=\\\$(awk "BEGIN {print ( \\\$c + \\\${box[2]} ) }")
                	done
                	a=\\\$(awk "BEGIN {print ( \\\$a / \\\${#endstates[@]} ) }")
                	b=\\\$(awk "BEGIN {print ( \\\$b / \\\${#endstates[@]} ) }")
                	c=\\\$(awk "BEGIN {print ( \\\$c / \\\${#endstates[@]} ) }")

                	a=\\\$(printf "%8.7f" \\\$a); b=\\\$(printf "%8.7f" \\\$b); c=\\\$(printf "%8.7f" \\\$c)
                	for lam in \\\${endstates[@]};do
                        	sed -e "s/ABOX/\\\${a}/g" -e "s/BBOX/\\\${b}/g" -e "s/CBOX/\\\${c}/g" inputs/\\\${lam}_eqNTP4.mdin.template > inputs/\\\${lam}_eqNTP4.mdin
                	done
                	sleep 1
        	fi
	done

done
EOF

fi


#############
#############

#submit independent jobs

if [ "\${twostate}" != "true" ]; then
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

count=-1; alllams=0
for stage in \\\${eqstage[@]}; do
        count=\\\$((\\\$count+1))
        lastcount=\\\$((\\\$count-1))
        if [ "\\\${stage}" == init ]; then continue; fi
        laststage=\\\${eqstage[\\\${lastcount}]}


        if [ "\\\${stage}" == "minTI" ];then
                for lam in \\\${lams[@]}; do
                	cp current/\\\${endstates[0]}_\\\${laststage}.rst7 current/\\\${lam}_init2.rst7
                done
                laststage="init2"
		alllams=1
        fi

	if [ \\\${alllams} -eq 0 ];then
		lam=\${endstates[0]}
                echo "Running \\\$stage for lambda \\\${lam}..."
                \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                # check if cpptraj is present
                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                cpptraj < center.in
                sleep 1
                mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7

        else

                for lam in \\\${lams[@]};do
                        echo "Running \\\$stage for lambda \\\${lam}..."
                        \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
                done

        fi

done
EOF

else
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

count=-1; alllams=0
for stage in \\\${eqstage[@]}; do
        count=\\\$((\\\$count+1))
        lastcount=\\\$((\\\$count-1))
        if [ "\\\${stage}" == init ]; then continue; fi
        laststage=\\\${eqstage[\\\${lastcount}]}

	if [ "\\\${stage}" == "minTI" ];then

		for lam in \\\${lams[@]}; do
                	if [ \\\$(bc -l <<< "\\\${lam} <= 0.5") -eq 1 ]; then
                        	cp current/\\\${endstates[0]}_\\\${laststage}.rst7 current/\\\${lam}_init2.rst7
                	else
                        	cp current/\\\${endstates[1]}_\\\${laststage}.rst7 current/\\\${lam}_init2.rst7
                	fi
		done
		laststage="init2"
		alllams=1

	fi

	if [ \\\${alllams} -eq 0 ];then

		for lam in \\\${endstates[@]};do
			source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                        # check if pmemd.cuda is present
                        if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
        		echo "Running \\\$stage for lambda \\\${lam}..."
			\\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
			cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
			# check if cpptraj is present
	        	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
        		cpptraj < center.in
        		sleep 1
			mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
		done

	else

		for lam in \\\${lams[@]};do
                        echo "Running \\\$stage for lambda \\\${lam}..."
                        \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                        # check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
                done

	fi
	
	if [ "\\\${stage}" == "eqP0" ]; then
		for lam in \\\${endstates[@]};do
			box=(\\\$(tail -1 \\\${lam}_\\\${laststage}.rst7))
		        a=\\\$(awk "BEGIN {print ( \\\$a + \\\${box[0]} ) }")
        		b=\\\$(awk "BEGIN {print ( \\\$b + \\\${box[1]} ) }")
        		c=\\\$(awk "BEGIN {print ( \\\$c + \\\${box[2]} ) }")
		done
		a=\\\$(awk "BEGIN {print ( \\\$a / \\\${#endstates[@]} ) }")
		b=\\\$(awk "BEGIN {print ( \\\$b / \\\${#endstates[@]} ) }")
		c=\\\$(awk "BEGIN {print ( \\\$c / \\\${#endstates[@]} ) }")

		a=\\\$(printf "%8.7f" \\\$a); b=\\\$(printf "%8.7f" \\\$b); c=\\\$(printf "%8.7f" \\\$c)
                for lam in \\\${endstates[@]};do
                        sed -e "s/ABOX/\\\${a}/g" -e "s/BBOX/\\\${b}/g" -e "s/CBOX/\\\${c}/g" inputs/\\\${lam}_eqNTP4.mdin.template > inputs/\\\${lam}_eqNTP4.mdin
                done
                sleep 1
	fi
		
done
EOF

fi

# slurm to submit all trials together
##########
##########
if [ "\${twostate}" != true ]; then
# submit group-ed jobs
        cat<<EOF > equilgroup_alltrials.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${ntrials}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --time=${wallclock}

top=\\\${PWD}
endstates=(\${endstates[@]})
lams=(\${lams[@]})
twostate=\${twostate}
eqstage=(${eqstagelist[*]})


# check if AMBERHOME is set
#if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

for trial in \\\$(seq 1 1 ${ntrials}); do

	if [ ! -d t\\\${trial} ];then mkdir t\\\${trial}; fi

	count=-1; alllams=0
	for stage in \\\${eqstage[@]}; do
        	count=\\\$((\\\${count}+1))
        	lastcount=\\\$((\\\${count}-1))
        	if [ "\\\${stage}" == "init" ]; then continue; fi
        	laststage=\\\${eqstage[\\\${lastcount}]}

        	if [ "\\\${stage}" == "minTI" ];then
                	for lam in \\\${lams[@]}; do
                        	cp t\\\${trial}/\\\${endstates[0]}_\\\${laststage}.rst7 t\\\${trial}/\\\${lam}_init2.rst7
                	done
                	laststage="init2"
                	alllams=1
        	fi

        	if [ \\\${alllams} -eq 0 ];then

                	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi

                	export LAUNCH="srun"
                	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

                	lam=\\\${endstates[0]}
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
                	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                	# check if cpptraj is present
                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                	cpptraj < center.in
                	sleep 1
                	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7

        	elif [ \\\${alllams} -eq 1 ] && [ "\\\${stage}" == "minTI" ];then
                	laststage="init2"
                	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
                	export LAUNCH="srun"
                	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda
                	for lam in \\\${lams[@]}; do
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
                         	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                         	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                         	cpptraj < center.in
                         	sleep 1
                         	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                	done
        	else
                	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda.MPI is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                	export LAUNCH="mpirun -np \\\${#lams[@]}"
                	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                	export MV2_ENABLE_AFFINITY=0
                	\\\${LAUNCH} \\\${EXE} -ng \\\${#lams[@]} -groupfile inputs/t\\\${trial}_\\\${stage}.groupfile

                	for lam in \\\${lams[@]};do
                        	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin t\\\${trial}/\\\${lam}_\\\${stage}.rst7
autoimage
trajout t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                        	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                        	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        	cpptraj < center.in
                        	sleep 1
                        	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                	done
        	fi
	done

done

EOF

else
# submit group-ed jobs
        cat<<EOF > equilgroup_alltrials.slurm
#!/bin/bash
#SBATCH --job-name="eq_${stA}~${stB}.slurm"
#SBATCH --output="eq_${stA}~${stB}.slurm.slurmout"
#SBATCH --error="eq_${stA}~${stB}.slurm.slurmerr"
#SBATCH --partition=${partition}
#SBATCH --nodes=${ntrials}
#SBATCH --ntasks-per-node=\${#lams[@]}
#SBATCH --gres=gpu:${ngpus}
#SBATCH --time=${wallclock}

top=\\\${PWD}
endstates=(\${endstates[@]})
lams=(\${lams[@]})
twostate=\${twostate}
eqstage=(${eqstagelist[*]})


# check if AMBERHOME is set
#if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

for trial in \\\$(seq 1 1 ${ntrials}); do

	if [ ! -d t\\\${trial} ];then mkdir t\\\${trial}; fi

	count=-1; alllams=0
	for stage in \\\${eqstage[@]}; do
        	count=\\\$((\\\${count}+1))
        	lastcount=\\\$((\\\${count}-1))
        	if [ "\\\${stage}" == "init" ]; then continue; fi
        	laststage=\\\${eqstage[\\\${lastcount}]}

        	if [ "\\\$stage" == "minTI" ]; then
                	for lam in \\\${lams[@]};do
                        	if [ \\\$(bc -l <<< "\\\${lam} <= 0.5") -eq 1 ]; then
                                	cp t\\\${trial}/\\\${endstates[0]}_\\\${laststage}.rst7 t\\\${trial}/\\\${lam}_init2.rst7
                        	else
                                	cp t\\\${trial}/\\\${endstates[1]}_\\\${laststage}.rst7 t\\\${trial}/\\\${lam}_init2.rst7
                        	fi
                	done
                	laststage="init2"
                	alllams=1
        	fi

        	if [ \\\${alllams} -eq 0 ];then

                	if [ "\\\$stage" == "min1" ] || [ "\\\$stage" == "min2" ] || [ "\\\$stage" == "minProt2" ] || [ "\\\$stage" == "minProt1" ] || [ "\\\$stage" == "minProt05" ] || [ "\\\$stage" == "minProt025" ] || [ "\\\$stage" == "minProt01" ]; then
                        	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                        	# check if pmemd.cuda is present
                        	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi

                        	export LAUNCH="srun"
                        	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

                        	for lam in \\\${endstates[@]};do
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
                                	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                                	# check if cpptraj is present
                                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                	cpptraj < center.in
                                	sleep 1
                                	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                        	done

                	else
                        	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                        	# check if pmemd.cuda.MPI is present
                        	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                        	export LAUNCH="mpirun -np \\\${#endstates[@]}"
                        	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                        	export MV2_ENABLE_AFFINITY=0
                        	\\\${LAUNCH} \\\${EXE} -ng \\\${#endstates[@]} -groupfile inputs/t\\\${trial}_\\\${stage}.groupfile

                        	for lam in \\\${endstates[@]};do
                                	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin t\\\${trial}/\\\${lam}_\\\${stage}.rst7
autoimage
trajout t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                                	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                	cpptraj < center.in
                                	sleep 1
                                	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                        	done
                	fi

        	elif [ "\\\${alllams}" == 1 ] && [ "\\\$stage" == "minTI" ]; then

                	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
                	export LAUNCH="srun"
                	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

                	for lam in \\\${lams[@]};do
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
                        	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                        	# check if cpptraj is present
                        	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        	cpptraj < center.in
                        	sleep 1
                        	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                	done
        	else
                	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda.MPI is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                	export LAUNCH="mpirun -np \\\${#lams[@]}"
                	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                	export MV2_ENABLE_AFFINITY=0
                	\\\${LAUNCH} \\\${EXE} -ng \\\${#lams[@]} -groupfile inputs/t\\\${trial}_\\\${stage}.groupfile

                	for lam in \\\${lams[@]};do
                        	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin t\\\${trial}/\\\${lam}_\\\${stage}.rst7
autoimage
trajout t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                        	source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                        	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                       	 	cpptraj < center.in
                        	sleep 1
                        	mv t\\\${trial}/\\\${lam}_\\\${stage}_centered.rst7 t\\\${trial}/\\\${lam}_\\\${stage}.rst7
                 	done
        	fi


        	if [ "\\\${stage}" == "eqP0" ]; then
                	for lam in \\\${endstates[@]};do
                        	box=(\\\$(tail -1 t\\\${trial}/\\\${lam}_\\\${stage}.rst7))
                        	a=\\\$(awk "BEGIN {print ( \\\$a + \\\${box[0]} ) }")
                        	b=\\\$(awk "BEGIN {print ( \\\$b + \\\${box[1]} ) }")
                        	c=\\\$(awk "BEGIN {print ( \\\$c + \\\${box[2]} ) }")
                	done
                	a=\\\$(awk "BEGIN {print ( \\\$a / \\\${#endstates[@]} ) }")
                	b=\\\$(awk "BEGIN {print ( \\\$b / \\\${#endstates[@]} ) }")
                	c=\\\$(awk "BEGIN {print ( \\\$c / \\\${#endstates[@]} ) }")

                	a=\\\$(printf "%8.7f" \\\$a); b=\\\$(printf "%8.7f" \\\$b); c=\\\$(printf "%8.7f" \\\$c)
                	for lam in \\\${endstates[@]};do
                        	sed -e "s/ABOX/\\\${a}/g" -e "s/BBOX/\\\${b}/g" -e "s/CBOX/\\\${c}/g" inputs/\\\${lam}_eqNTP4.mdin.template > inputs/\\\${lam}_eqNTP4.mdin
                	done
                	sleep 1
        	fi
	done
done
EOF

fi


##########
##########



if [ "\${twostate}" != true ]; then
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
endstates=(\${endstates[@]})
lams=(\${lams[@]})
twostate=\${twostate}
eqstage=(${eqstagelist[*]})


# check if AMBERHOME is set
#if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

if [ ! -d current ];then mkdir current; fi

count=-1; alllams=0
for stage in \\\${eqstage[@]}; do
        count=\\\$((\\\${count}+1))
        lastcount=\\\$((\\\${count}-1))
        if [ "\\\${stage}" == "init" ]; then continue; fi
        laststage=\\\${eqstage[\\\${lastcount}]}

        if [ "\\\${stage}" == "minTI" ];then
                for lam in \\\${lams[@]}; do
                	cp current/\\\${endstates[0]}_\\\${laststage}.rst7 current/\\\${lam}_init2.rst7
                done
                laststage="init2"
		alllams=1
        fi

	if [ \\\${alllams} -eq 0 ];then

        	source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
        	# check if pmemd.cuda is present
        	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi

        	export LAUNCH="srun"
        	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

		lam=\\\${endstates[0]}
        	echo "Running \\\$stage for lambda \\\${lam}..."
        	\\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                # check if cpptraj is present
                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                cpptraj < center.in
                sleep 1
		mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7

	elif [ \\\${alllams} -eq 1 ] && [ "\\\${stage}" == "minTI" ];then
                laststage="init2"
		source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                # check if pmemd.cuda is present
                if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
                export LAUNCH="srun"
                export EXE=\\\${AMBERHOME}/bin/pmemd.cuda
                for lam in \\\${lams[@]}; do
			echo "Running \\\$stage for lambda \\\${lam}..."
                        \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                         source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                         if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                         cpptraj < center.in
                         sleep 1
                         mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
                done
        else
                source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
               	# check if pmemd.cuda.MPI is present
                if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                export LAUNCH="mpirun -np \\\${#lams[@]}"
                export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                export MV2_ENABLE_AFFINITY=0
                \\\${LAUNCH} \\\${EXE} -ng \\\${#lams[@]} -groupfile inputs/current_\\\${stage}.groupfile

                for lam in \\\${lams[@]};do
                	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                        source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
                done
        fi

done

EOF

else
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
endstates=(\${endstates[@]})
lams=(\${lams[@]})
twostate=\${twostate}
eqstage=(${eqstagelist[*]})


# check if AMBERHOME is set
#if [ -z "\\\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi

if [ ! -d current ];then mkdir current; fi

count=-1; alllams=0
for stage in \\\${eqstage[@]}; do
        count=\\\$((\\\${count}+1))
        lastcount=\\\$((\\\${count}-1))
        if [ "\\\${stage}" == "init" ]; then continue; fi
        laststage=\\\${eqstage[\\\${lastcount}]}

	if [ "\\\$stage" == "minTI" ]; then
                for lam in \\\${lams[@]};do
                        if [ \\\$(bc -l <<< "\\\${lam} <= 0.5") -eq 1 ]; then
                                cp current/\\\${endstates[0]}_\\\${laststage}.rst7 current/\\\${lam}_init2.rst7
                        else
                                cp current/\\\${endstates[1]}_\\\${laststage}.rst7 current/\\\${lam}_init2.rst7
                        fi
		done
                laststage="init2"
		alllams=1
	fi

	if [ \\\${alllams} -eq 0 ];then

		if [ "\\\$stage" == "min1" ] || [ "\\\$stage" == "min2" ] || [ "\\\$stage" == "minProt2" ] || [ "\\\$stage" == "minProt1" ] || [ "\\\$stage" == "minProt05" ] || [ "\\\$stage" == "minProt025" ] || [ "\\\$stage" == "minProt01" ]; then
			source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi

			export LAUNCH="srun"
			export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

	        	for lam in \\\${endstates[@]};do
        	        	echo "Running \\\$stage for lambda \\\${lam}..."
                		\\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                		cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
				source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                		# check if cpptraj is present
                		if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                		cpptraj < center.in
               	 		sleep 1
                		mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
        		done

		else
			source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
                	# check if pmemd.cuda.MPI is present
                	if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                	export LAUNCH="mpirun -np \\\${#endstates[@]}"
                	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                	export MV2_ENABLE_AFFINITY=0
                	\\\${LAUNCH} \\\${EXE} -ng \\\${#endstates[@]} -groupfile inputs/current_\\\${stage}.groupfile

                	for lam in \\\${endstates[@]};do
                        	cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
				source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
               			if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                		cpptraj < center.in
                		sleep 1
                		mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
			done
		fi

	elif [ "\\\${alllams}" == 1 ] && [ "\\\$stage" == "minTI" ]; then

		source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
		# check if pmemd.cuda is present
		if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
		export LAUNCH="srun"
        	export EXE=\\\${AMBERHOME}/bin/pmemd.cuda

		for lam in \\\${lams[@]};do
        	        echo "Running \\\$stage for lambda \\\${lam}..."
        	        \\\${EXE} -O -p \\\${top}/unisc.parm7 -c current/\\\${lam}_\\\${laststage}.rst7 -i inputs/\\\${lam}_\\\${stage}.mdin -o current/\\\${lam}_\\\${stage}.mdout -r current/\\\${lam}_\\\${stage}.rst7 -ref current/\\\${lam}_\\\${laststage}.rst7
                                cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
			source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
			# check if cpptraj is present
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
			cpptraj < center.in
			sleep 1
			mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
		done
	else
		source /scratch/01969/abirg/Yorkgroup/Gitlab/bashrc_amberDD_gitlab_02-25-2021
		# check if pmemd.cuda.MPI is present
                if ! command -v \\\${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                export LAUNCH="mpirun -np \\\${#lams[@]}"
                export EXE=\\\${AMBERHOME}/bin/pmemd.cuda.MPI
                export MV2_ENABLE_AFFINITY=0
                \\\${LAUNCH} \\\${EXE} -ng \\\${#lams[@]} -groupfile inputs/current_\\\${stage}.groupfile

                for lam in \\\${lams[@]};do
                        cat <<EOF2 > center.in
parm \\\${top}/unisc.parm7
trajin current/\\\${lam}_\\\${stage}.rst7
autoimage
trajout current/\\\${lam}_\\\${stage}_centered.rst7
go
quit
EOF2
                        source /scratch/01969/abirg/Yorkgroup/bashrc_lbsr_dev
                        if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                        cpptraj < center.in
                        sleep 1
                        mv current/\\\${lam}_\\\${stage}_centered.rst7 current/\\\${lam}_\\\${stage}.rst7
                 done
	fi
	

	if [ "\\\${stage}" == "eqP0" ]; then
                for lam in \\\${endstates[@]};do
        	        box=(\\\$(tail -1 current/\\\${lam}_\\\${stage}.rst7))
                	a=\\\$(awk "BEGIN {print ( \\\$a + \\\${box[0]} ) }")
                	b=\\\$(awk "BEGIN {print ( \\\$b + \\\${box[1]} ) }")
                	c=\\\$(awk "BEGIN {print ( \\\$c + \\\${box[2]} ) }")
        	done
        	a=\\\$(awk "BEGIN {print ( \\\$a / \\\${#endstates[@]} ) }")
        	b=\\\$(awk "BEGIN {print ( \\\$b / \\\${#endstates[@]} ) }")
        	c=\\\$(awk "BEGIN {print ( \\\$c / \\\${#endstates[@]} ) }")

        	a=\\\$(printf "%8.7f" \\\$a); b=\\\$(printf "%8.7f" \\\$b); c=\\\$(printf "%8.7f" \\\$c)
                for lam in \\\${endstates[@]};do
                        sed -e "s/ABOX/\\\${a}/g" -e "s/BBOX/\\\${b}/g" -e "s/CBOX/\\\${c}/g" inputs/\\\${lam}_eqNTP4.mdin.template > inputs/\\\${lam}_eqNTP4.mdin
                done
                sleep 1
        fi
	
done
EOF

fi

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

