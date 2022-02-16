#!/bin/bash
#  This script helps prepare the executable setup_fe
#  It is recommended that all three repositories of the AMBER DD Boost package
# FE-MDEngine, FE-ToolKit, and FE-Workflow 
# are placed in a single folder, for example GitLab
#
#  For additional help reach out 
#  to abir.ganguly@rutgers.edu 


path=`pwd`

#####################################################################
# function to read compile directives
# the list "varlist" contains the variables that are defined from the input file
function read_input {

while read line; do
        varlist=(MDEngine ToolKit Workflow)
        IFS=$'\t| |=' read -ra args <<< $line
        if [[ "${args[0]}" =~ ^#.* ]]; then continue; fi
        keyword=${args[0]}; value=${args[1]}
        for var in ${varlist[@]}; do
                declare -n arr="$var"
                if [ "$var" == "$keyword" ]; then
                        arr=$value
                fi
        done


done < $1
}
######################################################################


# MAIN

path=`pwd`
# check if the file "setup_directives" is absent in the current directory
# if yes, write the default options
if [ ! -f "${path}/setup_directives" ]; then

	# default locations
	dir=$(basename $(dirname $(dirname `pwd`)))
	if [ -d $(dirname `pwd`)/FE-MDEngine ]; then
		MDEngine=$(dirname `pwd`)/FE-MDEngine/install_serial
	else
		read -p "Where is FE-MDEngine installed (serial)? " MDEngine
	fi
	if [ -d $(dirname `pwd`)/FE-ToolKit ]; then
		ToolKit=$(dirname `pwd`)/FE-ToolKit
	else
		read -p "What is FE-ToolKit installed? " ToolKit	
	fi
	if [ -d $(dirname `pwd`)/FE-Workflow ]; then
		Workflow=$(dirname `pwd`)/FE-Workflow
        else
                read -p "What is FE-Workflow installed? " Workflow
        fi
	MDEngine=$(echo "$(cd "$(dirname "${MDEngine}")"; pwd)/$(basename "${MDEngine}")")
	ToolKit=$(echo "$(cd "$(dirname "${ToolKit}")"; pwd)/$(basename "${ToolKit}")")
	Workflow=$(echo "$(cd "$(dirname "${Workflow}")"; pwd)/$(basename "${Workflow}")")
        cat << EOF2 > ${path}/setup_directives
MDEngine ${MDEngine}
ToolKit ${ToolKit}
Workflow ${Workflow}
EOF2

fi


# read setup directives
read_input setup_directives

# double check settings with user
printf "%s \n" "*************************************************************************"
printf "%s \n" "The location of MDEngine (serial installation)			is set to    ${MDEngine}"
printf "%s \n" "The location of ToolKit  					is set to    ${ToolKit}"
printf "%s \n" "The location of Workflow			 		is set to    ${Workflow}"
printf "%s \n\n"
printf "%s \n" "If these settings do not look correct, please modify \"${path}/setup_directives\" accordingly and re-run this script."
printf "%s \n" "*************************************************************************"

read -p "Would you like to continue (Y/N)? " cont
if [ "${cont}" != "Y" ] && [ "${cont}" != "y" ]; then exit 0; fi

if [ -f ${MDEngine}/amber.sh ]; then
	amberhome=${MDEngine}
elif [ ! -z ${AMBERHOME} ]; then
	amberhome=${AMBERHOME}
	printf "%s \n" "The amber.sh file at ${MDEngine} is missing." 
	printf "%s \n" "AMBERHOME is currently set to ${AMBERHOME}."
	printf "%s \n" "If these settings do not look correct, please modify \"${path}/setup_directives\" accordingly and re-run this script."
else
	read -p "The amber.sh file at ${MDEngine} is missing and AMBERHOME is currently not defined. Please provide an alternate location where AMBER is installed. " amberhome
	if [ ! -f ${amberhome}/amber.sh ]; then
		printf "%s \n" "The amber.sh file at ${amberhome} seems to be missing."
		printf "%s \n" "AMBER seems to be not installed at ${amberhome}"
		printf "%s \n" "Please install AMBER (${MDEngine} is recommended) and re-try setup"
		printf "%s \n" "Exiting"
		exit 0
	fi
fi

# write setup_fe
cat << EOF2 > ${Workflow}/bin/setup_fe
#!/bin/bash

##########################################
pathhere=\`pwd\`
path=\`pwd\`
pathTOWFToolKit=${Workflow}
pathTOFEToolKit=${ToolKit}
export PATH="\$PATH:\${Workflow}/bin"

# set AMBERHOME
source ${amberhome}/amber.sh
###########################################

##########################################
# load modules
source \${pathTOWFToolKit}/bin/function-read_input.sh 					# read input file
source \${pathTOWFToolKit}/bin/function-parse_input.sh					# parse input file
source \${pathTOWFToolKit}/bin/function-gen_lambda.sh						# generate lambda values
source \${pathTOWFToolKit}/bin/function-write_template_rbfe.sh				# template for rbfe calculations
source \${pathTOWFToolKit}/bin/function-write_template_rsfe.sh				# template for rsfe calculations
source \${pathTOWFToolKit}/bin/function-createbox.sh
source \${pathTOWFToolKit}/bin/function-preparePDBs.sh
source \${pathTOWFToolKit}/bin/function-analyze.sh


###########################################
# Main program

# change to working directory
cd \$path

	#read and parse input data
	parse_input \$1

	# initial array initialization
	for i in "\${!translist[@]}";do
		stA=\$(basename \${translist[\$i]}); stB="\${stA##*~}"; stA="\${stA%~*}"
		listA+=("\${stA}"); listB+=("\${stB}")
	done
	listligs+=(\${listA[@]} \${listB[@]})
	uniqueligs=(\$(echo "\${listligs[@]}" | tr ' ' '\\n' | sort -u | tr '\\n' ' '))

	#################################
	##BEGIN STAGE=setup
	if [ "\$stage" == "setup" ]; then
	#################################
	
		mkdir -p \${system}/setup

		# File Organization
		source \${pathTOWFToolKit}/bin/section-FileOrganization.sh

		# Atom Mapping
		source \${pathTOWFToolKit}/bin/section-AtomMapping.sh

		# Construction of MD boxes
		source \${pathTOWFToolKit}/bin/section-MDboxes.sh

		# setup mode 0 correspond to regular TI setup
		if [ "\${setupmode}" == 0 ]; then
			source \${pathTOWFToolKit}/bin/section-setupmodezero.sh
		fi
		# END of setupmode=0

		# setup mode 1 correspond to end-point ACES setup
		if [ "\${setupmode}" == 1 ]; then
			source \${pathTOWFToolKit}/bin/section-setupmodeone.sh
		fi
		# END of setupmode=0
	#################################
	##END STAGE=setup
	fi
	#################################

	#################################
	##BEGIN STAGE=analysis
	if [ "\$stage" == "analysis" ]; then
	#################################
		source \${pathTOWFToolKit}/bin/section-analysis.sh
	#################################
	##END STAGE=analysis
	fi
	#################################
cd \$pathhere

EOF2

chmod a+x ${Workflow}/bin/setup_fe
#printf "%s \n" "source ${MDEngine}/amber.sh" 			>  ${path}/FE-Workflow.bashrc
#printf "%s \n" "export PATH=\$PATH:${Workflow}/bin"		>> ${path}/FE-Workflow.bashrc
#printf "%s \n" "export PATH=\$PATH:${ToolKit}/local/bin" 	>> ${path}/FE-Workflow.bashrc

cat << EOF3 > ${path}/FE-Workflow.bashrc
#!/bin/bash

printf "%s \n\n" "Adding ${Workflow}/bin to \\\$PATH..."
export PATH=\$PATH:${Workflow}/bin

printf "%s \n\n" "Sourcing amber.sh available in ${MDEngine}..."
source ${MDEngine}/amber.sh

if [ -f "\$MODULEPATH/fetoolkit.module" ]; then 
	printf "%s \n\n" "Loading fetoolkit.module available in \$MODULEPATH..."
	module load fetoolkit
elif [ -f ${ToolKit}/fetoolkit.bashrc ]; then
	printf "%s \n\n" "Sourcing fetoolkit.bashrc available in ${ToolKit}..."
	source ${ToolKit}/fetoolkit.bashrc
else
	printf "%s \n" "Unable to load ${ToolKit} environment. Currently ${ToolKit} is unusable."
	printf "%s \n" "After successful ${ToolKit} installation, a ${ToolKit}/fetoolkit.bashrc file is created."
	printf "%s \n" "The ${ToolKit}/fetoolkit.bashrc file will also have information on how use ${ToolKit} as a module."
	printf "%s \n" "In that case, the file fetoolkit.module should be available in \$MODULEPATH."
fi
EOF3

printf "%s \n" " "
printf "%s \n" "***************************************************************************************************************************"
printf "%s \n" "Issue the following command before using the FE-Workflow (consider adding this line to your login startup script, e.g. ~/.bashrc)"
printf "%s \n" "source ${path}/FE-Workflow.bashrc"
printf "%s \n" "***************************************************************************************************************************"
printf "%s \n" " "


