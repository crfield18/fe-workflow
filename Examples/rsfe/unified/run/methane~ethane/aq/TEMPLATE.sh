#!/bin/bash
top=`pwd`

parmbase=unisc
rstbase=(stateA stateA)
endstates=(0.00000000)
lams=(0.00000000 0.05000000 0.10000000 0.15000000 0.20000000 0.25000000 0.30000000 0.35000000 0.40000000 0.45000000 0.50000000 0.55000000 0.60000000 0.65000000 0.70000000 0.75000000 0.80000000 0.85000000 0.90000000 0.95000000 1.00000000)

#generate initial configurations for lambda windows
mkdir -p current inputs
count=0
for lam in ${endstates[@]};do
        cp ${rstbase[$count]}.rst7  current/${lam}_init.rst7
        #cp ${parmbase}.parm7 unisc.parm7

        init=${lam}_init
        min1=${lam}_min1
        min2=${lam}_min2
	eqpre1P0=${lam}_eqpre1P0
	eqpre2P0=${lam}_eqpre2P0
        eqP0=${lam}_eqP0
        eqNTP4=${lam}_eqNTP4 
        eqV=${lam}_eqV
        eqP=${lam}_eqP
        eqA=${lam}_eqA
        ti=${lam}_ti

        cat<<EOF>inputs/${min1}.mdin
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

noshakemask	= ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
scmask1         = ''
scmask2         = ''

gti_syn_mass	= 0
/
EOF
	cat<<EOF>inputs/${min2}.mdin
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
scmask1         = ''
scmask2         = ''

gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/${eqpre1P0}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''
clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/${eqpre2P0}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/${eqP0}.mdin
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2 
ntf             = 1 
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass	= 0
/
EOF
        cat<<EOF>inputs/${eqNTP4}.mdin.template
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0

&end
 &ewald
   target_a=ABOX,
   target_b=BBOX,
   target_c=CBOX,
   target_n=500,

/
EOF
        cat<<EOF>inputs/${eqV}.mdin
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF

        cat<<EOF>inputs/${eqP}.mdin
&cntrl
imin            = 0
nstlim          = 500000 
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 2 
ntf             = 1 
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF

        cat<<EOF>inputs/${eqA}.mdin
&cntrl
imin            = 0
nstlim          = 2000000
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 2 
ntf             = 1 
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
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

count=$(($count+1))
done


for lam in ${lams[@]}; do

        init2=${lam}_init2
        minTI=${lam}_minTI
        eqpre1P0TI=${lam}_eqpre1P0TI
        eqpre2P0TI=${lam}_eqpre2P0TI
        eqP0TI=${lam}_eqP0TI
        eqATI=${lam}_eqATI
        preTI=${lam}_preTI
        ti=${lam}_ti
        anal=${lam}_analyze

	cat << EOF > inputs/${minTI}.mdin
&cntrl
imin            = 1
maxcyc          = 5000
ntmin           = 2
ntx             = 1
ntxo            = 1
ntpr            = 5


ifsc            = 1
icfe            = 1

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0

/
EOF
        cat<<EOF>inputs/${eqpre1P0TI}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/${eqpre2P0TI}.mdin
&cntrl
imin            = 0
nstlim          = 5000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF
        cat<<EOF>inputs/${eqP0TI}.mdin
&cntrl
imin            = 0
nstlim          = 100000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF
        cat << EOF > inputs/${eqATI}.mdin
&cntrl
imin            = 0
nstlim          = 500000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
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
	cat << EOF > inputs/${preTI}.mdin
&cntrl
imin            = 0
nstlim          = 2000000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
gti_syn_mass    = 0
/
EOF
        cat<< EOF >inputs/${ti}.mdin
&cntrl
imin            = 0
nstlim          = 20
numexchg        = 135000
dt              = 0.001
irest           = 0
ntx             = 1
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 10000
ntwr            = 5000
ntpr            = 1000
cut             = 10
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
mbar_states     = ${#lams[@]}
EOF
ilam=0
for l in ${lams[@]}; do
     ilam=$(( ${ilam} + 1 ))
     echo "mbar_lambda(${ilam})     = ${l}" >> inputs/${ti}.mdin
done
        cat<< EOF >> inputs/${ti}.mdin

noshakemask     = ''

clambda         = ${lam}
timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6

gremd_acyc      = 1
/
 &ewald
 /
EOF

        cat<< EOF >inputs/${lam}_analyze.mdin
&cntrl
imin            = 6
nstlim          = 20
numexchg        = 135000
dt              = 0.001
irest           = 1
ntx             = 5
ntxo            = 1
ntc             = 2
ntf             = 1
ntwx            = 0
ntwr            = 0
ntpr            = 1
cut             = 10
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
mbar_states     = ${#lams[@]}
EOF
ilam=0
for l in ${lams[@]}; do
     ilam=$(( ${ilam} + 1 ))
     echo "mbar_lambda(${ilam})     = ${l}" >> inputs/${ti}.mdin
done
        cat<< EOF >> inputs/${ti}.mdin

noshakemask     = ''

clambda         = ${lam}

timask1         = ''
timask2         = ''
crgmask         = ""
scmask1         = ''
scmask2         = ''
scalpha         = 0.2
scbeta          = 50

gti_cut         = 1
gti_output      = 1
gti_add_sc      = 2
gti_scale_beta  = 0
gti_cut_sc_on   = 8
gti_cut_sc_off  = 10
gti_lam_sch     = 1
gti_ele_sc      = 1
gti_vdw_sc      = 1
gti_cut_sc      = 2
gti_ele_exp     = 2
gti_vdw_exp     = 6
/
 &ewald
 /

EOF
done

truncate -s0 inputs/eqATI.groupfile inputs/preTI.groupfile inputs/ti.groupfile
for lam in ${lams[@]};do
	minTI=${lam}_minTI
        eqpre1P0TI=${lam}_eqpre1P0TI
        eqpre2P0TI=${lam}_eqpre2P0TI
        eqP0TI=${lam}_eqP0TI
	eqATI=${lam}_eqATI
	preTI=${lam}_preTI
	ti=${lam}_ti

	cat<<EOF>> inputs/eqATI.groupfile
-O -p unisc.parm7 -c current/${eqP0TI}.rst7 -i inputs/${eqATI}.mdin -o current/${eqATI}.mdout -r current/${eqATI}.rst7 -x current/${eqATI}.nc -ref current/${minTI}.rst7
EOF
	cat<<EOF>> inputs/preTI.groupfile
-O -p unisc.parm7 -c current/${eqATI}.rst7 -i inputs/${preTI}.mdin -o current/${preTI}.mdout -r current/${preTI}.rst7 -x current/${preTI}.nc -ref current/${eqATI}.rst7
EOF
	cat<<EOF>> inputs/ti.groupfile
-O -p unisc.parm7 -c current/${preTI}.rst7 -i inputs/${ti}.mdin -o current/${ti}.mdout -r current/${ti}.rst7 -x current/${ti}.nc -ref current/${preTI}.rst7
EOF

done

# slurm to submit all trials
###############
###############

cat<<EOF > run_alltrials.slurm
#!/bin/bash
#SBATCH --job-name="eq_methane~ethane.slurm"
#SBATCH --output="eq_methane~ethane.slurm.slurmout"
#SBATCH --error="eq_methane~ethane.slurm.slurmerr"
#SBATCH --partition=general-long-gpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=${#lams[@]}
#SBATCH --gres=gpu:8
#SBATCH --time=1-00:00:00

top=\${PWD}
endstates=(${endstates[@]})
lams=(${lams[@]})
eqstage=(init min1 min2 eqpre1P0 eqpre2P0 eqP0 eqV eqP eqA minTI eqpre1P0TI eqpre2P0TI eqP0TI eqATI preTI)


# check if AMBERHOME is set
if [ -z "\${AMBERHOME}" ]; then echo "AMBERHOME is not set" && exit 0; fi
# check if cpptraj is present
if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi

EXE=\${AMBERHOME}/bin/pmemd.cuda

for trial in \$(seq 1 1 4); do
	
	if [ ! -d t\${trial} ];then mkdir t\${trial}; fi

	count=-1; alllams=0
	for stage in \${eqstage[@]}; do
        	count=\$((\$count+1))
        	lastcount=\$((\$count-1))
		if [ "\${stage}" == "init" ] || [ "\${stage}" == "eqpre1P0TI" ] || [ "\${stage}" == "eqpre2P0TI" ] || [ "\${stage}" == "eqP0TI" ]; then continue; fi
        	laststage=\${eqstage[\${lastcount}]}

		if [ "\${stage}" == "minTI" ];then alllams=1; fi

		
        	if [ \${alllams} -eq 0 ];then


			# check if pmemd.cuda is present
                        if ! command -v \${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
			export LAUNCH="srun"
			export EXE=\${AMBERHOME}/bin/pmemd.cuda

                	lam=${endstates[0]}
                	echo "Running \$stage for lambda \${lam}..."
                	\${EXE} -O -p \${top}/unisc.parm7 -c t\${trial}/\${lam}_\${laststage}.rst7 -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${lam}_\${laststage}.rst7
                	cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin t\${trial}/\${lam}_\${stage}.rst7
autoimage
trajout t\${trial}/\${lam}_\${stage}_centered.rst7
go
quit
EOF2
                	# check if cpptraj is present
                	if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                	cpptraj < center.in
                	sleep 1
                	mv t\${trial}/\${lam}_\${stage}_centered.rst7 t\${trial}/\${lam}_\${stage}.rst7

		elif [ \${alllams} -eq 1 ] && [ "\${stage}" == "minTI" ];then
			# check if pmemd.cuda is present
			if ! command -v \${AMBERHOME}/bin/pmemd.cuda &> /dev/null; then echo "pmemd.cuda is missing." && exit 0; fi
			export LAUNCH="srun"
			export EXE=\${AMBERHOME}/bin/pmemd.cuda
			for i in \${!lams[@]}; do
				lam=\${lams[\$i]}
				if [ "\${i}" -eq 0 ]; then
					init=\${endstates[0]}_eqA.rst7
 				else
					init=\${lams[\$((\$i-1))]}_eqP0TI.rst7
				fi			

				echo "Running \$stage for lambda \${lam}..."

				stage=minTI
                                \${EXE} -O -p \${top}/unisc.parm7 -c t\${trial}/\${init} -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${init}
                                sleep 1

                                laststage=minTI; stage=eqpre1P0TI
                                \${EXE} -O -p \${top}/unisc.parm7 -c t\${trial}/\${lam}_\${laststage}.rst7 -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${lam}_\${laststage}.rst7
                                sleep 1

                                laststage=eqpre1P0TI; stage=eqpre2P0TI
                                \${EXE} -O -p \${top}/unisc.parm7 -c t\${trial}/\${lam}_\${laststage}.rst7 -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${lam}_\${laststage}.rst7
                                sleep 1

                                laststage=eqpre2P0TI; stage=eqP0TI
                                \${EXE} -O -p \${top}/unisc.parm7 -c t\${trial}/\${lam}_\${laststage}.rst7 -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${lam}_\${laststage}.rst7
                                sleep 1

                                cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin t\${trial}/\${lam}_\${stage}.rst7
autoimage
trajout t\${trial}/\${lam}_\${stage}_centered.rst7
go
quit
EOF2
                                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                cpptraj < center.in
                                sleep 1
                                mv t\${trial}/\${lam}_\${stage}_centered.rst7 t\${trial}/\${lam}_\${stage}.rst7
                        done
                        laststage=eqP0TI
		else
                        # check if pmemd.cuda.MPI is present
                        if ! command -v \${AMBERHOME}/bin/pmemd.cuda.MPI &> /dev/null; then echo "pmemd.cuda.MPI is missing." && exit 0; fi

                        export LAUNCH="mpirun -np \${#lams[@]}"
                        export EXE=\${AMBERHOME}/bin/pmemd.cuda.MPI
                        export MV2_ENABLE_AFFINITY=0
                        \${LAUNCH} \${EXE} -ng \${#lams[@]} -groupfile inputs/t\${trial}_\${stage}.groupfile

			for lam in \${lams[@]};do
				cat <<EOF2 > center.in
parm \${top}/unisc.parm7
trajin t\${trial}/\${lam}_\${stage}.rst7
autoimage
trajout t\${trial}/\${lam}_\${stage}_centered.rst7
go
quit
EOF2
                                if ! command -v cpptraj &> /dev/null; then echo "cpptraj is missing." && exit 0; fi
                                cpptraj < center.in
                                sleep 1
                                mv t\${trial}/\${lam}_\${stage}_centered.rst7 t\${trial}/\${lam}_\${stage}.rst7
                        done
		fi
	done

        # run production
        EXE=\${AMBERHOME}/bin/pmemd.cuda.MPI
        echo "running replica ti"
        mpirun -np \${#lams[@]} \${EXE} -rem 3 -ng \${#lams[@]} -groupfile inputs/t\${trial}_ti.groupfile


	cat << EOFP > extract.py
#!/usr/bin/env python3

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
        if abs(param.box[3]-109.471219)<1.e-4 and            abs(param.box[4]-109.471219)<1.e-4 and            abs(param.box[5]-109.471219)<1.e-4:
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
       if abs(param.box[3]-109.471219)<1.e-4 and           abs(param.box[4]-109.471219)<1.e-4 and           abs(param.box[5]-109.471219)<1.e-4:
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

    parser = argparse.ArgumentParser \
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

	mkdir -p ../vac/t\${trial} ../vac/inputs
	for lam in \${lams[@]};do
        	./extract.py -p unisc.parm7 -c t\${trial}/\${lam}_preTI.rst7 -m '!:WAT,Na+,K+,Cl-' -o ../vac/t\${trial}/\${lam}_init
        	sed -e 's/nstlim.*/nstlim          = 500000/g' -e 's/restraint_wt.*/restraint_wt    = 0/g' -e 's/nmropt.*/nmropt          = 0/g' -e '59,65d' -e "s/clambda.*/clambda         = /g" -e 's/irest.*/irest           = 0/g' -e 's/ntx.*/ntx             = 1/g' inputs/0.00000000_eqA.mdin > ../vac/inputs/\${lam}_preTI.mdin
        	sed -e 's/ntb.*/ntb             = 1/g' -e '/barostat.*/d' -e '/ntp.*/d' -e '/pres0.*/d' -e '/taup.*/d' -e '/numexchg/d' -e '/gremd_acyc/d' -e 's/nstlim.*/nstlim          = 2700000/g' inputs/\${lam}_ti.mdin > ../vac/inputs/\${lam}_ti.mdin
                sed -i -e 's/ntwr.*.=*.*/ntwr            = 5000\nntpr            = 1000/' ../vac/inputs/\${lam}_ti.mdin
	done

	cd ../vac
        	cd t\${trial}
                	mv 0.00000000_init.parm7 ../unisc.parm7
                	rm *.parm7
        	cd ../

        	stage=preTI; laststage=init
        	EXE=\${AMBERHOME}/bin/pmemd.cuda
        	LAUNCH="srun --exclusive -N 1 -n 1 -c 1 --gres=gpu:1"
        	for lam in \${lams[@]};do
                	\${LAUNCH} \${EXE} -O -p unisc.parm7 -c t\${trial}/\${lam}_\${laststage}.rst7 -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${lam}_\${laststage}.rst7 &
        	done
        	wait
	cd ../

	for dir in vac; do
		cd \${dir}
			stage=ti; laststage=preTI
			for lam in \${lams[@]};do
				\${LAUNCH} \${EXE} -O -p unisc.parm7 -c t\${trial}/\${lam}_\${laststage}.rst7 -i inputs/\${lam}_\${stage}.mdin -o t\${trial}/\${lam}_\${stage}.mdout -r t\${trial}/\${lam}_\${stage}.rst7 -ref t\${trial}/\${lam}_\${laststage}.rst7 &
			done
		cd ../
	done
	wait

done

EOF

