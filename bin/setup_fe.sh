#!/bin/bash

##########################################
path=`pwd`
pathTObin=${path}/bin			  		  	# path to bin folder
export PATH="$PATH:${pathTObin}"
###########################################

##########################################
# load modules
source ${pathTObin}/function-read_input.sh 					# read input file
source ${pathTObin}/function-parse_input.sh					# parse input file
source ${pathTObin}/function-gen_lambda.sh					# generate lambda values
source ${pathTObin}/function-write_template_rbfe.sh				# template for rbfe calculations
source ${pathTObin}/function-write_template_rsfe.sh				# template for rsfe calculations
source ${pathTObin}/function-createbox.sh
source ${pathTObin}/function-preparePDBs.sh


###########################################
# Main program 

#read and parse input data
parse_input $1
####

#################################
##BEGIN STAGE=setup
if [ "$stage" == "setup" ]; then
#################################

	mkdir -p ${system}/setup

	for i in "${!translist[@]}";do
		stA=$(basename ${translist[$i]}); stB="${stA##*~}"; stA="${stA%~*}"
		chgA=$(basename ${chargelist[$i]}); chgB="${chgA##*~}"; chgA="${chgA%~*}"
		listA+=("${stA}"); listB+=("${stB}")
		listchgA+=("${chgA}"); listchgB+=("${chgB}")
	done
	listligs+=(${listA[@]} ${listB[@]})
	listchgs+=(${listchgA[@]} ${listchgB[@]})
	uniqueligs=($(echo "${listligs[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
	for i in "${!uniqueligs[@]}";do
        	for j in "${!listligs[@]}";do
                	if [ "${uniqueligs[$i]}" == "${listligs[$j]}" ]; then uniquechgs[$i]="${listchgs[$j]}" && break; fi
        	done
        done

	# File Organization
	source ${pathTObin}/section-FileOrganization.sh

	# Atom Mapping
	source ${pathTObin}/section-AtomMapping.sh

	# Construction of MD boxes
	source ${pathTObin}/section-MDboxes.sh

	# setup mode 0 correspond to regular TI setup
	if [ "${setupmode}" == 0 ]; then
		source ${pathTObin}/section-setupmodezero.sh	
	fi
	# END of setupmode=0

	# setup mode 1 correspond to end-point ACES setup
	if [ "${setupmode}" == 1 ]; then
		source ${pathTObin}/section-setupmodeone.sh	
	fi
	# END of setupmode=0
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
		for s in ${slist[@]}; do
			cd ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
				echo "Submitting TI equil jobs in ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}"
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
		for s in ${slist[@]}; do
			cd ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
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
		for s in ${slist[@]}; do
			cd ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
                                echo "Submitting TI production jobs in ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}"
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
		for s in ${slist[@]}; do
			if [ ! -d ${path}/${system}/${protocol}/run/${stA}~${stB}/${s} ]; then echo "Folder ${path}/${system}/${protocol}/run/${stA}~${stB}/${s} missing" && continue; fi
                        cd ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}
                                for(( t=1;t<=${ntrials};t++));do
                                      	 for i in "${!lams[@]}";do
                                                	if ! grep -Eq 'EPtot      = **************|NaN' t${t}/${lams[$i]}_ti.mdout; then
                                               		if grep -Eq 'Final Performance Info' t${t}/${lams[$i]}_ti.mdout; then continue; fi
                                                	fi
                                                	echo "check ${path}/${system}/${protocol}/run/${stA}~${stB}/${s}/t${t}/${lams[$i]}_ti.mdout"
                                        	done
                                done
                        cd ${path}
                done
        done
#################################
##END STAGE=check-TI
fi
#################################


