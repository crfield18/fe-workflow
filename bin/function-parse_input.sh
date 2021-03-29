# parse input from user
function parse_input {
        if [ ! -z "$1" ]; then
                if [ "$1" == "-h" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ] ; then
                        cat << EOFN > input.template
system=Tyk2                     # system
translist=(ejm42~ejm54 ejm42~ejm55 ejm55~ejm54) 	# list of transformations
chargelist=(0~0 0~0 0~0)				# formal charge on ligands in the respective transformations. chargelist should correspond to translist

path_to_input=../initial        # path to folder containing input configuration files 
inputformat=parm		# parm/pdb
nlambda=4                       # number of lambda windows
protocol=unified                # unified protocol for TI

mapmethod=0                     # 0 --> MCSS/ 1 --> MCSS-E
mapinspect=0                    # 0 --> auto/ 1 --> mapinspect=true/ 2 --> mapinspect=checked
mapnetwork=true                 # generate network-wide consistent sc maps

boxbuild=1			# 0 --> dont build boxes for complexes/ 1 --> build boxes/ 2 --> build boxes with same number of water and ions
boxbufcom=15
boxbufaq=20
ionconc=0.15

ntrials=3

cutoff=10                       # non-bonded cutoff
mincyc=2000                     # max minimization cycles in lambda window equilibration
nstlimnvt=298000                # length of 0->298K NVT heating of lambda windows
nstlimnpt=500000                # length of NPT equilibration of lambda windows. Ingored if ticalc is set to 'abfe'
repex=true                      # true/false corresponding to use of replica exchange in TI simulations
nstlimti=5000                   # length of TI simulations
numexchgti=1000                 # number of exchanges in replica exchange TI simulations. if repex=false, numexchgti is ignored
hmr=false                       # "true/false" If hmr=true, dt (timestep) is set to 4fs
scalpha=0.5                     # scalpha
scbeta=1.0                      # scbeta
gti_add_sc=5                    # gti_add_sc
gti_scale_beta=1                # gti_scale_beta
gti_cut=1                       # gti_cut
gti_cut_sc_on=8                 # gti_cut_sc_on
gti_cut_sc_off=10               # gti_cut_sc_off
gti_lam_sch=1                   # gti_lam_sch
gti_ele_sc=1                    # gti_ele_sc
gti_vdw_sc=1                    # gti_vdw_sc
gti_cut_sc=2                    # gti_cut_sc
gti_ele_exp=2                   # gti_ele_exp
gti_vdw_exp=2                   # gti_vdw_exp

twostate=true                   # enables the two state model for generation of initial lambda window configurations.

ticalc=rbfe                     # "rbfe -> relative binding free energy/rsfe -> relative solvation free energy"
stage=setup                     # "setup/run-equil/check-equil/run-TI/check-TI"
setupmode=0                     # 0 --> regular TI/ 1 --> end-point ACES/ 2 --> TI from end-point ACES

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
        if [ -z "${AMBERHOME}" ]; then printf "\n\nAMBERHOME is not set\n\n" && exit 0; fi
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then printf "\n\ncpptraj is missing.\n\n" && exit 0; fi
        # check if parmed is present
        if ! command -v parmed &> /dev/null;  then printf "\n\nparmed is missing.\n\n" && exit 0; fi


        # check input file parameters
        if [ "${protocol}" != "unified" ]; then printf "\n\nScript currently supports only \"unified\" protocol\n\n" && exit 0; fi

	if [ "${inputformat}" != "parm" ] && [ "${inputformat}" != "pdb" ]; then printf "\n\n\"inputformat\" should be set to \"parm\" or \"pdb\"\n\n" && exit 0; fi

	if [ "${inputformat}" == "pdb" ]; then
		if [ "${#chargelist[@]}" -ne ${#translist[@]} ]; then printf "\n\n\"if inputformat is set to pdb, chargelist, a list containing molecular charges of the ligands involved in the transformations listed in translist, should also be provided. For example, an entry A~B in translist should have a corresponding entry CA~CB in chargelist, where CA, CB are net charges of the ligands in A and B, respectively. The lengths of the lists translist and chargelist should be equal." && exit 0; fi   
	fi

	if [ "${mapmethod}" -lt 0 ] && [ "${mapmethod}" -gt 2 ]; then printf "\n\n\"mapmethod\" should be set to 0 for \"MCSS\", 1 for \"MCSS-E\" algorithm for atom-mapping\n\n, 2 for \"MCSS-E\" variant applicable to linear transformations that involve change in mass of atoms" && exit 0; fi

	if [ "${mapinspect}" -lt 0 ] || [ "${mapinspect}" -gt 3 ]; then printf "\n\n\"mapinspect\" should be set to 0 (for creating file infrastructure WITHOUT map manual inspection), 1 (for interrupting setup to allow for map manual inspection), or 2 (for continuing setup assuming map manual inspection has been completed)" && exit 0; fi

        if [ "${mapnetwork}" != "true" ] && [ "${mapnetwork}" != "false" ]; then printf "\n\n\"mapnetwork\" should be set to \"true\" or \"false\"\n\n" && exit 0; fi

        if [ "${boxbuild}" != "skip" ] && [ "${boxbuild}" == 0 ] && [ "${boxbuild}" == 1 ] && [ "${boxbuild}" == 2 ]; then printf "\n\n\"boxbuild\" should be set to \"skip\", \"0\", \"1\" or \"2\". When \"twostate\" = true, \"boxbuild\" is ignored\n\n" && exit 0; fi

        if [ "${repex}" != "true" ] && [ "${repex}" != "false" ]; then printf "\n\n\"repex\" should be set to \"true\" or \"false\"\n\n" && exit 0; fi

        if [ "${hmr}" != "true" ] && [ "${hmr}" != "false" ]; then printf "\n\n\"hmr\" should be set to \"true\" or \"false\"\n\n"   && exit 0; fi

        if [ "${gti_add_sc}" -lt 1 ] || [ "${gti_add_sc}" -gt 5 ]; then printf "\n\nAcceptable values for \"gti_add_sc\" are 1, 2(Recommended), 3, 4, and 5\n\n" && exit 0; fi

        if [ "${gti_lam_sch}" != 0 ] && [ "${gti_lam_sch}" != 1 ]; then printf "\n\nAcceptable values of \"gti_lam_sch\" are 0(default) and 1\n\n" && exit 0; fi

        if [ "${gti_lam_sch}" -eq 1 ]; then
                if [ "${gti_ele_sc}" != 1 ] || [ "${gti_vdw_sc}" != 1 ]; then printf "\n\nIf gti_lam_sch is set to 1, both gti_ele_sc and gti_vdw_sc should be set to 1\n\n" && exit 0; fi
        fi

        if [ "${gti_scale_beta}" != 0 ] && [ "${gti_scale_beta}" != 1 ]; then printf "\n\nAcceptable values of \"gti_scale_beta\" are 0(default) and 1\n\n" && exit 0; fi

        if [ "${gti_scale_beta}" -eq 1 ]; then
                if [ -z "${gti_ele_exp}" ] || [ -z "${gti_vdw_exp}"  ] || [ -z "$scalpha" ] || [ -z "$scbeta" ]; then printf "\n\nIf gti_scale_beta is set to 1, gti_ele_exp gti_vdw_exp scalpha scbeta all needs to be defined. Recommended values are gti_lam_sch = 1,  gti_ele_sc = 1, gti_vdw_sc = 1, gti_ele_exp = 2, gti_vdw_exp = 2,  gti_scale_beta = 1, scalpha = 0.5, scbeta = 1.0\n\n" && exit 0; fi
        fi

        if [ "${gti_cut}" != 0 ] && [ "${gti_cut}" != 1 ]; then printf "\n\nAcceptable values of \"gti_cut\" are 0(default) and 1\n\n" && exit 0; fi

        if [ "${gti_cut_sc}" != 0 ] && [ "${gti_cut_sc}" != 1 ] && [ "${gti_cut_sc}" != 2 ]; then printf "\n\nAcceptable values of \"gti_cut_sc\" are 0(default), 1, 2\n\n" && exit 0; fi

        if [ "${gti_cut_sc}" -eq 1 ] || [ "${gti_cut_sc}" -eq 2 ]; then
                if [ -z "${gti_cut_sc_on}" ] || [ -z "${gti_cut_sc_off}" ]; then printf "\n\nif \"gti_cut_sc_on\" is set to 1 or 2, gti_cut_sc_on and gti_cut_sc_off must be defined.\n\n" && exit 0; fi
        fi

        if [ "${cutoff}" -lt "${gti_cut_sc_on}" ] || [ "${cutoff}" -lt "${gti_cut_sc_off}" ] || [ "${gti_cut_sc_off}" -lt "${gti_cut_sc_on}" ]; then print "\n\nShould be \"cutoff\" >= \"gti_cut_sc_off\" > \"gti_cut_sc_on\" \n\n" && exit 0; fi

	if [ "${ticalc}" == "rbfe" ]; then slist=(com aq); else slist=(aq); fi

        # ensure path_to_input is absolute and check if input directories are present
        if [[ "${path_to_input}" != /* ]]; then path_to_input=${path}"/"${path_to_input}; fi

	# check input files
	for i in "${!translist[@]}";do
                stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
                listA+=("${stA}"); listB+=("${stB}")
        done
        listligs+=(${listA[@]} ${listB[@]})
        uniqueligs=($(echo "${listligs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
	if [ "${inputformat}" == "parm" ]; then
		parmmissing=0
		for i in "${!uniqueligs[@]}";do
			molname=${uniqueligs[$i]}
			if [ ! -f ${path_to_input}/${system}/${molname}.parm7 ] ||  [ ! -f ${path_to_input}/${system}/${molname}.rst7 ]; then
				printf "\n\ninputformat is set to parm. ${molname}.parm7 or ${molname}.rst7 is missing\n\n"
				parmmissing=1
			fi
		done
		if [ "${parmmissing}" -eq 1 ]; then exit 0; fi
	else
		pdbmissing=0
		for i in "${!uniqueligs[@]}";do
                        molname=${uniqueligs[$i]}
                        if [ ! -f ${path_to_input}/${system}/${molname}.pdb ]; then
                                printf "\n\ninputformat is set to pdb. ${molname}.pdb is missing\n\n"
                                pdbmissing=1
                        fi
                done
                if [ "${pdbmissing}" -eq 1 ]; then exit 0; fi


	fi	
}


