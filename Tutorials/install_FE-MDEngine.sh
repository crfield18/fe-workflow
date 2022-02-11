#!/bin/bash
#
#  This scripts attempts to install AMBER using cmake with
#  some default options that should work for most cases.
#  However, in specific cases, this script needs to be 
#  modified. 
#
#  For information on how to get cmake, visit this page:
#  https://ambermd.org/pmwiki/pmwiki.php/Main/CMake-Quick-Start
#
#  For information on common options for cmake, visit this page:
#  http://ambermd.org/pmwiki/pmwiki.php/Main/CMake-Common-Options
#
#  For additional help related to installation, you can reach out 
#  to abir.ganguly@rutgers.edu 


path=`pwd`

#####################################################################
# function to read compile directives
# the list "varlist" contains the variables that are defined from the input file
function read_input {

while read line; do
	varlist=(amber_src amber_build amber_install cuda install)
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

# check if the file "compile_directives" is absent in the current directory
# if yes, write the default options
if [ ! -f "${path}/compile_directives" ]; then

	if [ -z "${amber_src}" ]; then
        	read -p "What is the location of the FE-MDEngine directory? " amber_src
	fi
	# default options
	amber_src=$(echo "$(cd "$(dirname "${amber_src}")"; pwd)/$(basename "${amber_src}")")
	amber_build=${path}/build
	amber_install=${path}/install
	cuda=no
	install=serial
	cat << EOF > ${path}/compile_directives
amber_src ${amber_src}
amber_build ${amber_build}
amber_install ${amber_install}
cuda no
install serial
EOF

fi

# read input directives
read_input compile_directives

# double check settings with user
printf "%s \n" "*************************************************************************"
printf "%s \n" "The AMBER source directory is set to 	${amber_src}"
printf "%s \n" "The AMBER build directory is set to 	${amber_build}"
printf "%s \n" "The AMBER install directory is set to 	${amber_install}"

if [ "${cuda,,}" == "no" ]; then 
	printf "%s \n" "This installation will be without CUDA"
	cudavalue="FALSE"
elif [ "${cuda,,}" == "yes" ]; then
	printf "%s \n" "This installation will be with CUDA"
	cudavalue="TRUE"
else
	printf "%s \n" "The flag \"cuda\" should be set to either no or yes. Exiting installation." 
	exit 0
fi	

if [ "${install,,}" == "serial" ]; then 
        printf "%s \n" "This will be a SERIAL install" 
	mpivalue="FALSE"
elif [ "${install,,}" == "parallel" ]; then
        printf "%s \n" "This will be a PARALLEL install"    
	mpivalue="TRUE"
else
        printf "%s \n" "The flag \"install\" should be set to either serial or parallel."
	printf "%s \n" "Please correct \"${path}/compile_directives\" and re-run script. Exiting installation."     
        exit 0
fi

printf "%s \n\n"
printf "%s \n" "If these settings do not look correct, please modify \"${path}/compile_directives\" accordingly and re-run this script."
printf "%s \n" "*************************************************************************"

read -p "Would you like to continue (Y/N)" cont
if [ "${cont}" != "Y" ] && [ "${cont}" != "y" ]; then exit 0; fi

# ensure paths are absolute
amber_src=$(echo "$(cd "$(dirname "${amber_src}")"; pwd)/$(basename "${amber_src}")")
amber_build=$(echo "$(cd "$(dirname "${amber_build}")"; pwd)/$(basename "${amber_build}")")
amber_install=$(echo "$(cd "$(dirname "${amber_install}")"; pwd)/$(basename "${amber_install}")")

# create the intended build and install directories
mkdir -p ${amber_build} ${amber_install}

# go to the build directory and run cmake
cd ${amber_build}
	if [ `uname -s|awk '{print $1}'` = "Darwin" ]; then

		#macOS
  		if [ -x /Applications/CMake.app/Contents/bin/cmake ]; then
     			cmake=/Applications/CMake.app/Contents/bin/cmake
  		else
     			cmake=cmake
  		fi

  		${cmake} 	-DCMAKE_INSTALL_PREFIX=${amber_install} \
    				-DCOMPILER=CLANG  \
				-DBLA_VENDOR=Apple \
    				-DMPI=${mpivalue} \
				-DCUDA=${cudavalue} \
				-DINSTALL_TESTS=TRUE \
    				-DDOWNLOAD_MINICONDA=TRUE \
				-DMINICONDA_USE_PY3=TRUE \
				${amber_src} 2>&1 | tee cmake.log

	else
		#Linux OS
		cmake 		-DCMAKE_INSTALL_PREFIX=${amber_install} \
				-DCOMPILER=GNU  \
				-DMPI=${mpi_value} \
				-DCUDA=${cuda_value} \
				-DINSTALL_TESTS=TRUE \
				-DDOWNLOAD_MINICONDA=TRUE \
				-DMINICONDA_USE_PY3=TRUE \
				${amber_src} 2>&1 | tee  cmake.log
	fi
	
	if [ ! -s cmake.log ]; then
  		printf "%s \n" "Error:  No cmake.log file was created: you may need to edit run_cmake"
  		exit 0
	fi

	printf "\n\n\n\n"
	printf "%s \n" "*****************************************************************************************************************"
	printf "%s \n" "Done with compilation."
	printf "%s \n" "If the cmake build report looks OK, issue the following commands:"
	printf "\n"
	printf "%s \n" "1) cd ${amber_build}"
	printf "%s \n" "2) make install"
	printf "%s \n" "3) source ${amber_install}/amber.sh"
	printf "\n"
	printf "%s \n" "Consider adding the last line "source ${amber_install}/amber.sh" to your login startup script, e.g. ~/.bashrc"
	printf "\n"
	printf "%s \n" "*****************************************************************************************************************"

cd ${path}


