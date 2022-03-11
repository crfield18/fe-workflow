# parse input from user
function parse_input {
        if [ ! -z "$1" ]; then
                if [ "$1" == "-h" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ] ; then
                        cat << EOFN > ${path}/input.template
#######################################################################
# 1. NAME OF SYSTEM, INPUT STRUCTURES, and TYPE OF CALCULATION
#######################################################################
# Directory that contains within it a subdirectory names "system"
# containing initial structure/parameter files.
# For example, a single folder "initial" may contain multiple
# subdirectories containing initial files for different systems
path_to_input=initial


# Subdirectory containing initial structure/parameter files
# of "system"
# For example,
#system=smallMols
system=CDK2


# setupmode determines the calculation.
# setupmode=0 --> regular TI
# setupmode=1 --> end-point ACES
# setupmode=2 --> TI from end-point ACES.
# currenty, only setupmode=0 is implemented. 
setupmode=0                     

# ticalc determines the nature of TI calculation.
# ticalc=rbfe --> relative binding free energy
# ticalc=rsfe --> relative solvation free energy
# ticalc=asfe --> absolute solvation free energy
ticalc=rbfe                     

# stage controls the action of the script
# stage=setup 		--> setup of TI  simulations
# stage=analysis 	--> analysis of TI simulations using FE-ToolKit
stage=setup                     


# List of desired transformations or edges
# For example, RBFE or RSFE calculations should have a list
# in which each entry consists of two molnames separated by the
# character "~". Initial structure/parameter files of these
# molnames should be provided in ${path_to_input}/${system}
# For RBFE calculations, the PDB file of protein-ligand complex
# and mol2,lib,frcmod files of ligand are expected.
# For RSFE calculations, mol2,lib,frcmod files of ligand are expected.
# example,
# translist=(1h1q~1h1r 1h1q~1h1s)
#
# For ASFE calculations, "translist" should contain a list of
# molnames.
# mol2,lib,frcmod files of these molnames are expected in
# ${path_to_input}/${system}
# example,
# translist=(mobley_1527293 mobley_3034976)
#translist=(1h1q~1h1r 1h1r~1h1s 1h1s~1oiu 1oiu~1h1q 1h1r~1oiu 1h1s~1h1q)
translist=(1h1q~1h1r 1h1r~1h1s)
#######################################################################


#######################################################################
# 2. ATOM MAPPING 
#######################################################################
#
# mapmethod determines the algorithm using which cc and sc regions
# will be determined.
# mapmethod=0 --> MCSS
# mapmethod=1 --> MCSS-E
# mapmethod=2 --> MCSS-E2
mapmethod=0

# mapinspect determines if there is need of manual inspection of the
# atom maps
# mapinspect=0 --> no-inspection. generate the atom maps using
# algorithm specified by "mapmethod", and then proceed to generate
# file infrastructure
# mapinspect=1 --> manual inspection. stop after generating the
# atom maps.
# mapinspect=2 --> resume generation of file infrastructure assuming
# map inspection has been completed.
# mapinspect=2 expects necessary atom map files to be present in the "setup" folder
mapinspect=0


# mapnetwork determines whether network-wide consistent cc and sc regions
# will be generated.
# mapnetwork=true ensures that in a given network of transformations, cc and sc
# regions of each ligand is identical in every transformation in which is participates
mapnetwork=false
#######################################################################

#######################################################################
# 3. MD BOX BUILDING
#######################################################################
#
# boxbuild determines if and how MD boxes will be built
# "skip" 	--> skip box building
# 0      	--> for RBFE calculations, do not build boxes for "complex" state, 
#		    only for "aqueous" state.
# 1 		--> build boxes for both "complex" and "aqueous" states
# 		    for RSFE and ASFE calculations, boxbuild=0 and boxbuild=1 are 
#		    identical.
# 2 		--> build boxes for both "complex" and "aqueous" states with same 
#		    number of water and ions
boxbuild=2
boxbufcom=16                    # MD box buffer for "complex" states
boxbufaq=20                     # MD box buffer for "aqueous" states
ionconc=0.15                    # Ion concentration in MD box
pff=ff14SB                      # Protein force field
lff=gaff2                       # Ligand forcefield
wm=tip4pew                      # Water model
mdboxshape=cubic                # Shape of MD box
#######################################################################


#######################################################################
# 4. GENERAL SETTINGS OF TI SIMULATIONS
#######################################################################
#
nlambda=25                      # number of lambda windows
lamschedule=yes
lams=(0 0.176834 0.229764 0.269379 0.302697 0.33229 0.359436 0.384886 0.40913 0.432518 0.455318 0.477748 0.5 0.522252 0.544682 0.567482 0.59087 0.615114 0.640564 0.66771 0.697303 0.730621 0.770236 0.823166 1)
protocol=unified                # unified protocol for TI

ntrials=3                       # Number of independent trials

cutoff=10                       # non-bonded cutoff
repex=true
nstlimti=20                     # length of TI simulations
numexchgti=250000               # number of exchanges in replica exchange TI simulations. if repex=true
hmr=false
notrajectory=true               # when true, no output trajectories are generated
scalpha=0.5                     # scalpha
scbeta=1.0                      # scbeta
gti_add_sc=5
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

# twostate determines the protocol to be used for equilibration of protein-ligand complex systems.
# twostate=false directs script to setup the equilibration file infrastructure
# in an "1-state" way in which for a given transformation P:A --> P:B, only the P:A structure is
# considered and the ligand B is superimposed on ligand A.
# twostate=true directs script to setup the equilibration file infrastructure
# in a "2-state" way in which for a given transformation P:A --> P:B, both P:A and P:B structures
# considered and represents the two end states.
twostate=true
bidirection_com=false
bidirection_aq=false
#######################################################################



#######################################################################
# 5. SETTINGS RELATED TO JOB SUBMISSION
#######################################################################
#
# job submission related
partition=general-long-gpu      # name of specific partition on HPC. Use "null" is not relevant
nnodes=1                        # number of nodes to be used for each transformation
ngpus=8                         # number of gpus/node to be used for each transformation
wallclock=3-00:00:00            # wallclock for individual jobs
#######################################################################

#######################################################################
# 6. ANALYSIS
#######################################################################
#
# analysis related
# path to production runs. default path_to_input="system"/"protocol"/run
# exptdatafile is an optional text file containing experimental free energies.
# exptdatafile can be set to "skip" or if provided, should be a file containing 2 columns.
# col 1 should be ligand name (identical to ligand name in translist) and col2 should be
# relative experimental free energy
#
path_to_data=data
exptdatafile=cdk2_expt.dat
bar=true
ccc=true
start=0.0
stop=100.0
check_convergence=true
#######################################################################

EOFN
                        echo "Script expects a file named \"input\" in working directory. Check input.template for details."
                        exit 0
                fi
        else
                if [ ! -f ${path}/input ]; then 
			printf "Script expects a file named \"input\" in \$pathTOWRKDIR directory\n" 
			printf "currently set to ${path}\n" 
			printf "Run script with -h/-help to generate template input file \n" 
			exit 0 
		fi
        fi

        # read input file
        read_input ${path}/input
	
        # check if AMBERHOME is set
        if [ -z "${AMBERHOME}" ]; then printf "\n\nAMBERHOME is not set\n\n" && exit 0; fi
        # check if cpptraj is present
        if ! command -v cpptraj &> /dev/null; then printf "\n\ncpptraj is missing.\n\n" && exit 0; fi
        # check if parmed is present
        if ! command -v parmed &> /dev/null;  then printf "\n\nparmed is missing.\n\n" && exit 0; fi


        # check input file parameters

	# lambda schedule related
	if [ "${lamschedule}" == "yes" ]; then
		printf "\n\n User defined lambda schedule... \n"
		if [ "${nlambda}" -ne "${#lams[@]}" ]; then printf "\n\nThe list of lambda values should be provided as a list using \"lams\". The value of \"nlambda\" should be equal to the total number of entries in the list \"lams\" \n\n" && exit 0; fi
		# format lambda values
		for l in ${!lams[@]}; do
			lams[$l]=$(printf "%0.8f" ${lams[$l]})
		done
	else
		printf "\n\n Auto generated lambda schedule... \n"
		lams=($(gen_lambdas $nlambda))
	fi
	printf "\n The following lambda schedule will be used \n"
	printf "%s\n\n" "${lams[*]}"

        if [ "${protocol}" != "unified" ]; then printf "\n\nScript currently supports only \"unified\" protocol\n\n" && exit 0; fi

	if [ "${mapmethod}" -lt 0 ] && [ "${mapmethod}" -gt 2 ]; then printf "\n\n\"mapmethod\" should be set to 0 for \"MCSS\", 1 for \"MCSS-E\" algorithm for atom-mapping\n\n, 2 for \"MCSS-E\" variant applicable to linear transformations that involve change in mass of atoms" && exit 0; fi

	if [ "${mapinspect}" -lt 0 ] || [ "${mapinspect}" -gt 3 ]; then printf "\n\n\"mapinspect\" should be set to 0 (for creating file infrastructure WITHOUT map manual inspection), 1 (for interrupting setup to allow for map manual inspection), or 2 (for continuing setup assuming map manual inspection has been completed)" && exit 0; fi

        if [ "${mapnetwork}" != "true" ] && [ "${mapnetwork}" != "false" ]; then printf "\n\n\"mapnetwork\" should be set to \"true\" or \"false\"\n\n" && exit 0; fi

        if [ "${boxbuild}" != "skip" ] && [ "${boxbuild}" == 0 ] && [ "${boxbuild}" == 1 ] && [ "${boxbuild}" == 2 ]; then printf "\n\n\"boxbuild\" should be set to \"skip\", \"0\", \"1\" or \"2\". When \"twostate\" = true, \"boxbuild\" is ignored\n\n" && exit 0; fi

        if [ "${repex}" != "true" ] && [ "${repex}" != "false" ]; then printf "\n\n\"repex\" should be set to \"true\" or \"false\"\n\n" && exit 0; fi

        if [ "${hmr}" != "true" ] && [ "${hmr}" != "false" ]; then printf "\n\n\"hmr\" should be set to \"true\" or \"false\"\n\n"   && exit 0; fi

        if [ "${notrajectory}" != "true" ] && [ "${notrajectory}" != "false" ]; then printf "\n\n\"notrajectory\" should be set to \"true\" or \"false\"\n\n"   && exit 0; fi

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

	if [ "${ticalc}" != "rbfe" ] && [ "${ticalc}" != "rsfe" ]  && [ "${ticalc}" != "asfe" ]; then printf "\n\n\"ticalc\" should be set to either \"rbfe\" or \"rsfe\" or \"asfe\"\n\n" && exit 0; fi
	if [ "${ticalc}" == "asfe" ] &&  [ "${twostate}" == "true" ]; then printf "\n\n\"ticalc\"=asfe is not compatible with \"twostate\"=true \n\n" && exit 0; fi 
	if [ "${twostate}" != "true" ] && [ "${twostate}" != "false" ]; then printf "\n\n\"twostate\" can either be \"true\" or \"false\" \n\n" && exit 0; fi 
	if [ "${bidirection_aq}" != "true" ] && [ "${bidirection_aq}" != "false" ]; then printf "\n\n\"bidirection_aq\" can either be \"true\" or \"false\" \n\n" && exit 0; fi 
	if [ "${bidirection_com}" != "true" ] && [ "${bidirection_com}" != "false" ]; then printf "\n\n\"bidirection_com\" can either be \"true\" or \"false\" \n\n" && exit 0; fi 
	
	if [ "${stage}" == "setup" ]; then

		# ensure path_to_input is absolute and check if input directories are present
		if [[ "${path_to_input}" != /* ]]; then path_to_input=${path}"/"${path_to_input}; fi
		
		########################
		########################
		# check input files
		if [ "${ticalc}" == "rbfe" ] || [ "${ticalc}" == "rsfe" ]; then
			for i in "${!translist[@]}";do
		        	stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		        	listA+=("${stA}"); listB+=("${stB}")
				translistrev+=("${stB}~${stA}")
			done
			listligs+=(${listA[@]} ${listB[@]})
			uniqueligs=($(echo "${listligs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
		else
			for i in "${!translist[@]}";do
				if grep -q '~' <<< "${translist[$i]}"; then printf "\n\n The character '~' is present is ${translist[$i]}. For \"ticalc\"=asfe, translist should contain a list of ligand molecules.\n\n" && exit 0; fi
				listligs+=("${translist[$i]}")
				uniqueligs=($(echo "${listligs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
			done
		fi
		########################
		
		if [ "${ticalc}" == "rbfe" ]; then
			pdbmissing=0; ligmissing=0; slist=(com aq)
			for i in "${!uniqueligs[@]}";do
		        	molname=${uniqueligs[$i]}
		        	if [ ! -f ${path_to_input}/${system}/${molname}.pdb ]; then
					printf "\n\n!!!! ERROR !!!!\n\n"
		                	printf "\n\n ${molname}.pdb is missing\n\n"
		                	pdbmissing=1
		        	elif ! grep -Fq 'LIG' ${path_to_input}/${system}/${molname}.pdb; then
					printf "\n\n!!!! ERROR !!!!\n\n"
					printf "\n\n ${molname}.pdb does not contain residue named "LIG". The ligand molecule in ${molname}.pdb must be named "LIG" \n\n"
		                	pdbmissing=1
				fi
		
		        	if [ ! -f ${path_to_input}/${system}/${molname}_0.mol2 ] || [ ! -f ${path_to_input}/${system}/${molname}_0.frcmod ] || [ ! -f ${path_to_input}/${system}/${molname}_0.lib ]; then
					printf "\n\n!!!! ERROR !!!!\n\n"
		                	printf "\n\n One or more of ligand parameter files present in ${molname}.pdb is missing. \n"
					printf "Parameter files of the ligand present in ${molname}.pdb must be provided as ${molname}_0.mol2, ${molname}_0.frcmod, ${molname}_0.lib. \n"
					printf "If ${molname}.pdb contain additional nonstandard residues, paramater files associated with these nonstandard residues must also be provided as \n"
					printf "${molname}_1.frcmod/${molname}_1.lib, ${molname}_2.frcmod/${molname}_2.lib etc \n\n"
		                	ligmissing=1
				elif ! grep -Fq 'LIG' ${path_to_input}/${system}/${molname}_0.mol2; then
					printf "\n\n!!!! ERROR !!!!\n\n"
					printf "\n\n resname of the ligand in ${molname}_0.mol2 must be \"LIG\" \n\n"
					ligmissing=1
				elif ! grep -Fq 'LIG' ${path_to_input}/${system}/${molname}_0.lib; then
					printf "\n\n!!!! ERROR !!!!\n\n"
			        	printf "\n\n resname 0of the ligand in ${molname}_0.lib must be \"LIG\" \n\n"
					ligmissing=1
		        	fi
		
				if [ "${ligmissing}" -eq 1 ]; then exit 0; fi
		
				frcmods=$(ls -l ${path_to_input}/${system}/${molname}_?.frcmod | wc -l)
				libs=$(ls -l ${path_to_input}/${system}/${molname}_?.lib | wc -l)
			
				if [ "${frcmods}" -ne "${libs}" ]; then 
					printf "\n\n!!!! ERROR !!!!\n\n"
					printf "\n\n Each ${molname}_?.lib file should have a corresponding ${molname}_?.frcmod file \n\n"
					exit 0
				else
					frcmodlist+=($frcmods)
					liblist+=($libs)
				fi
			done
			if [ "${pdbmissing}" -eq 1 ] || [ "${ligmissing}" -eq 1 ]; then exit 0; fi
		fi
		
		if [ "${ticalc}" == "rsfe" ] || [ "${ticalc}" == "asfe" ]; then
		        pdbmissing=0; ligmissing=0; slist=(aq)
		        for i in "${!uniqueligs[@]}";do
		                molname=${uniqueligs[$i]}
		                if [ ! -f ${path_to_input}/${system}/${molname}_0.mol2 ] || [ ! -f ${path_to_input}/${system}/${molname}_0.frcmod ] || [ ! -f ${path_to_input}/${system}/${molname}_0.lib ]; then
		                        printf "\n\n!!!! ERROR !!!!\n\n"
		                        printf "\n\n One or more of ligand parameter files present in ${molname}.pdb is missing. \n"
		                        printf "Parameter files of the ligand present in ${molname}.pdb must be provided as ${molname}_0.mol2, ${molname}_0.frcmod, ${molname}_0.lib. \n"
		                        printf "If ${molname}.pdb contain additional nonstandard residues, paramater files associated with these nonstandard residues must also be provided as \n"
		                        printf "${molname}_1.frcmod/${molname}_1.lib, ${molname}_2.frcmod/${molname}_2.lib etc \n\n"
		                        ligmissing=1
		                elif ! grep -Fq 'LIG' ${path_to_input}/${system}/${molname}_0.mol2; then
		                        printf "\n\n!!!! ERROR !!!!\n\n"
		                        printf "\n\n resname of the ligand in ${molname}_0.mol2 must be \"LIG\" \n\n"
		                        ligmissing=1
		                elif ! grep -Fq 'LIG' ${path_to_input}/${system}/${molname}_0.lib; then
		                        printf "\n\n!!!! ERROR !!!!\n\n"
		                        printf "\n\n resname of the ligand in ${molname}_0.lib must be \"LIG\" \n\n"
		                        ligmissing=1
		                fi
		
				if [ "${ligmissing}" -eq 1 ]; then exit 0; fi
		
		                frcmods=$(ls -l ${path_to_input}/${system}/${molname}_?.frcmod | wc -l)
		                libs=$(ls -l ${path_to_input}/${system}/${molname}_?.lib | wc -l)
		
		                if [ "${frcmods}" -ne "${libs}" ]; then
		                        printf "\n\n!!!! ERROR !!!!\n\n"
		                        printf "\n\n Each ${molname}_?.lib file should have a corresponding ${molname}_?.frcmod file \n\n"
		                        exit 0
		                else
		                        frcmodlist+=($frcmods)
		                        liblist+=($libs)
		                fi
		        done
		fi
	fi

	# analysis keywords
	if [ "${stage}" == "analysis" ]; then
                if [ ! -d ${pathTOFEToolKit}/local/bin ]; then
                        printf "%s\n\n" "!!! ERROR !!!"
                        printf "%s\n" "${pathTOFEToolKit}/local/bin is missing. The location of \"\${pathTOFEToolKit}\" may not be correct. Exiting analysis."
                        exit 0
                fi
		if [ ! -d "${path_to_data}" ]; then 
			printf "\n\n!!!! ERROR !!!!\n\n"
			printf "\n\n!!!! ${path_to_data} does not exist \n\n"
			exit 0
		fi
		if [ ! -f "${exptdatafile}" ] && [ "${exptdatafile}" != "skip" ]; then
                        printf "\n\n!!!! ERROR !!!!\n\n"
                        printf "\n\n!!!! ${exptdatafile} does not exist \n\n"
                        exit 0
                fi
		if [ "${bar}" != "true" ] && [ "${bar}" != "false" ]; then 
			printf "\n\n!!!! ERROR !!!!\n\n"
			printf "\n\n${bar} must be set to \"true\" or \"false\"\n\n"
		fi
                if [ "${ccc}" != "true" ] && [ "${ccc}" != "false" ]; then
                        printf "\n\n!!!! ERROR !!!!\n\n"
                        printf "\n\n${ccc} must be set to \"true\" or \"false\"\n\n"
                fi

		if (( $(echo "$start < 0" | bc -l) )) || (( $(echo "$stop < 0" | bc -l) )) || (( $(echo "$start > 100" | bc -l) )) || (( $(echo "$stop > 100" | bc -l) )) || (( $(echo "$start > $stop" | bc -l) )); then
			printf "\n\n!!!! ERROR !!!!\n\n"
			printf "\n\n \"start\" and \"stop\" should have values between 0 to 100 \n\n"
			printf "\n\n \"start\" should have a value less than \"stop\" \n\n"
		fi
                if [ "${check_convergence}" != "true" ] && [ "${check_convergence}" != "false" ]; then
                        printf "\n\n!!!! ERROR !!!!\n\n"
                        printf "\n\n${check_convergence} must be set to \"true\" or \"false\"\n\n"
                fi
		if [ "${ticalc}" == "rbfe" ]; then 
			slist=(com aq)
		elif [ "${ticalc}" == "rsfe" ]; then
			slist=(aq vac)
		else
			slist=(aq)
		fi
	fi

		

	############################
	############################

}


