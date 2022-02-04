#!/bin/bash

# check if ${pathTOWRKDIR}, ${pathTOWFToolKit}, ${pathTOFEToolKit} paths are set

if [ -z ${pathTOWRKDIR} ]; then
	cat << EOF

***** The variable \${pathTOWRKDIR} has not been set. It should be set to the intended working directory
By default, \${pathTOWRKDIR} is set to the path of current directory.

To change the path of working directory, either re-run this make script after setting \${pathTOWRKDIR}, OR,
Manually change the variable 'path' at the top of the 'setup_fe' script

EOF
	pathTOWRKDIR=`pwd`
else
	cat << EOF

***** The current working directory (\${pathTOWRKDIR}) is set to ${pathTOWRKDIR}

To change the path of working directory, either re-run this make script after setting \${pathTOWRKDIR}, OR,
Manually change the variable 'path' at the top of the 'setup_fe' script.

EOF

fi


if [ -z ${pathTOWFToolKit} ]; then
	cat << EOF

***** The variable \${pathTOWFToolKit} has not been set. 
\${pathTOWFToolKit} should be set to the path of the 'alchemical_fe' folder available within the AMBER_DD_BOOST package.
By default, \${pathTOWFToolKit} is set to the path of current directory.

To change the path of working directory, either re-run this make script after setting \${pathTOWFToolKit}, OR,
Manually change the variable 'pathTOWFToolKit' at the top of the 'setup_fe' script

EOF
	pathTOWFToolKit=`pwd`
else
        cat << EOF

***** The current path to workflow toolkit (\${pathTOWFToolKit}) is set to ${pathTOWFToolKit}
\${pathTOWFToolKit} should be set to the path of the 'alchemical_fe' folder available within the AMBER_DD_BOOST package.

To change the path to workflow toolkit (alchemical_fe) either re-run this make script after setting \${pathTOWFToolKit}, OR,
Manually change the variable 'pathTOWFToolKit' at the top of the 'setup_fe' script

EOF

fi


if [ -z ${pathTOFEToolKit} ]; then
        cat << EOF

***** The variable \${pathTOFEToolKit} has not been set.
\${pathTOFEToolKit} should be set to the path of the 'FE-ToolKit' folder available within the AMBER_DD_BOOST package.
By default, \${pathTOFEToolKit} is set to the path of current directory.

To change the path of working directory, either re-run this make script after setting \${pathTOFEToolKit}, OR,
Manually change the variable 'pathTOFEToolKit' at the top of the 'setup_fe' script

EOF
        pathTOWFToolKit=`pwd`
else
        cat << EOF

***** The current path to FE-ToolKit (\${pathTOFEToolKit}) is set to ${pathTOFEToolKit}
\${pathTOFEToolKit} should be set to the path of the 'FE-ToolKit' folder available within the AMBER_DD_BOOST package.

To change the path to FE-ToolKit either re-run this make script after setting \${pathTOFEToolKit}, OR,
Manually change the variable 'pathTOFEToolKit' at the top of the 'setup_fe' script

EOF

fi


if [ ! -d $pathTOWFToolKit/bin ] ||  [ ! -d $pathTOFEToolKit/local/bin ]; then
	printf "!!! ERROR bin directories located in \$pathTOWFToolKit OR \$pathTOFEToolKit/local are missing! \n"
	printf "!!! ERROR check paths and re-run make script again. \n"

else
	# write setup_fe
	cat << EOF > setup_fe
#!/bin/bash

##########################################
pathhere=\`pwd\`
path=${pathTOWRKDIR}
pathTOWFToolKit=${pathTOWFToolKit}
pathTOFEToolKit=${pathTOFEToolKit}
export PATH="\$PATH:\${pathTOWFToolKit}/bin"
###########################################

##########################################
# load modules
source \${pathTOWFToolKit}/bin/function-read_input.sh 					# read input file
source \${pathTOWFToolKit}/bin/function-parse_input.sh					# parse input file
source \${pathTOWFToolKit}/bin/function-gen_lambda.sh					# generate lambda values
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

EOF

	chmod a+x setup_fe
fi

