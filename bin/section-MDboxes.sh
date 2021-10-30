# generation of MD boxes
##########################
##########################
##########################
if [ "${ticalc}" == "rbfe" ]; then
	cd $system/setup
		# one-state setup
		if [ "${twostate}" == "false" ]; then
			printf "\n\n setting up files in \"1-state\" mode... \n\n"
			if [ "${boxbuild}" == "skip" ]; then
				flag=0; flag2=0
				for i in "${!translist[@]}"; do
					for s in aq com; do
						if [ "${s}" == "aq" ]; then bidirection=${bidirection_aq}; else bidirection=${bidirection_com}; fi
						if [ ! -f "${translist[$i]}_${s}.parm7" ] || [ ! -f "${translist[$i]}_${s}.rst7" ] || [ ! -f "${translist[$i]}.scmask1" ] || [ ! -f "${translist[$i]}.scmask2" ] || [ ! -f "${translist[$i]}.timask1" ] || [ ! -f "${translist[$i]}.timask2" ]; then
							printf "\n\n One of the following files are missing from $system/setup...\n"
							printf "${translist[$i]}_${s}.parm7, ${translist[$i]}_${s}.rst7, ${translist[$i]}.scmask1 ${translist[$i]}.scmask2 ${translist[$i]}.timask1 ${translist[$i]}.timask2 \n"
							flag=1
						fi

						if [ "${bidirection}" == "true" ]; then
							if [ ! -f "${translistrev[$i]}_${s}.parm7" ] || [ ! -f "${translistrev[$i]}_${s}.rst7" ] || [ ! -f "${translistrev[$i]}.scmask1" ] || [ ! -f "${translistrev[$i]}.scmask2" ] || [ ! -f "${translistrev[$i]}.timask1" ] || [ ! -f "${translistrev[$i]}.timask2" ]; then
								printf "\n\n One of the following files are missing from $system/setup...\n"
								printf "${translistrev[$i]}_${s}.parm7, ${translistrev[$i]}_${s}.rst7, ${translistrev[$i]}.scmask1 ${translistrev[$i]}.scmask2 ${translistrev[$i]}.timask1 ${translistrev[$i]}.timask2 \n"
								flag2=1
							fi

						fi
					done
				done
				if [ "${flag}" -eq 1 ]; then 
					printf "\n\n ${boxbuild} was set to 0. Above files must be present in $system/setup folder. Exiting...\n"
					exit 0 
				fi
				if [ "${flag2}" -eq 1 ]; then 
					printf "\n\n ${boxbuild} was set to 0 and bidirection_aq/_com was set to \"true\". Above files must be present in $system/setup folder. Exiting...\n"
					exit 0 
				fi
			else
			
				if [ "${boxbuild}" == "0" ]; then
					printf "\n\nWater and Ions will be added only to \"aq\" systems. Existing solvent configurations of \"com\" systems will be used.\n\n"
				elif [ "${boxbuild}" == "1" ]; then
					printf "\n\nWater and Ions will be added to both \"aq\" and \"com\" systems. In case of \"com\" systems, existing solvent configurations will be stripped before resolvating.\n\n"
				elif [ "${boxbuild}" == "2" ]; then
					printf "\n\nFor all transformations identical number of water and ions will be added to both \"aq\" and \"com\" systems. In case of \"com\" systems, existing solvent configurations will be stripped before resolvating. \n\n"
				fi
				for s in aq com; do
					if [ "${s}" == "aq" ]; then 
						rbuf=${boxbufaq}
						bidirection=${bidirection_aq}
				       	else 
						rbuf=${boxbufcom}
						bidirection=${bidirection_com}
				       	fi
					printf "\n\nBuilding topology+parameter files for \"${s}\" systems. Water and Ions being added...\n\n"
					
					if [ "${bidirection}" == "true" ]; then
						printf "\n\nbidirection_${s} set to true. One state simulations for \"${s}\" systems will be setup for both directions... \n\n"
					fi

					create_box_rbfe "${pff}" "${lff}" "${wm}" "${boxbuild}" "${mdboxshape}" "${rbuf}" "${ionconc}" "molname-ligname.mapping" "${s}" "${ticalc}" "${bidirection}" "${translist[@]}"
					printf "\n\nDone...\n\n"
				done
			fi
		fi

		# two-state setup
		if [ "${twostate}" == "true" ]; then
			printf "\n\n setting up files in \"2-state\" mode... \n\n"
                        if [ "${boxbuild}" == "skip" ]; then
                                flag=0; flag2=0
                                for i in "${!translist[@]}"; do
                                        if [ ! -f "${translist[$i]}-1_aq.parm7" ] || [ ! -f "${translist[$i]}-1_aq.rst7" ] || [ ! -f "${translist[$i]}-2_aq.parm7" ] || [ ! -f "${translist[$i]}-2_aq.rst7" ][ ! -f "${translist[$i]}-1_com.parm7" ] || [ ! -f "${translist[$i]}-1_com.rst7" ] || [ ! -f "${translist[$i]}-2_com.parm7" ] || [ ! -f "${translist[$i]}-2_com.rst7" ] || [ ! -f "${translist[$i]}.scmask1" ] || [ ! -f "${translist[$i]}.scmask2" ] || [ ! -f "${translist[$i]}.timask1" ] || [ ! -f "${translist[$i]}.timask2" ]; then
                                                printf "\n\n One of the following files are missing from $system/setup...\n"
                                                printf "${translist[$i]}-1_aq.parm7, ${translist[$i]}-1_aq.rst7, ${translist[$i]}-2_aq.parm7, ${translist[$i]}-2_aq.rst7, ${translist[$i]}-1_com.parm7, ${translist[$i]}-1_com.rst7, ${translist[$i]}-2_com.parm7, ${translist[$i]}-2_com.rst7, ${translist[$i]}.scmask1 ${translist[$i]}.scmask2 ${translist[$i]}.timask1 ${translist[$i]}.timask2 \n"                
                                                flag=1
                                        fi
						
					if [ "${bidirection_aq}" == "true" ]; then
						printf "\n\n bidirection_aq=true will be ignored in \"2-state\" setup of \"com\" systems... \n\n"
					fi

					if [ "${bidirection_com}" == "true" ]; then
						printf "\n\n bidirection_com=true will be ignored in \"2-state\" setup of \"com\" systems... \n\n"
					fi	
                                done 
				if [ "${flag}" -eq 1 ]; then
                                        printf "\n\n ${boxbuild} was set to 0. Above files must be present in $system/setup folder. Exiting...\n"
                                        exit 0
                                fi
                        else
                        
                                if [ "${boxbuild}" == "0" ]; then
                                        printf "\n\nWater and Ions will be added only to \"aq\" systems. Existing solvent configurations of \"com\" systems will be used.\n\n"                  
                                elif [ "${boxbuild}" == "1" ]; then
                                        printf "\n\nWater and Ions will be added to both \"aq\" and \"com\" systems. In case of \"com\" systems, existing solvent configurations will be stripped before resolvating.\n\n"
                                elif [ "${boxbuild}" == "2" ]; then
                                        printf "\n\nFor all transformations identical number of water and ions will be added to both \"aq\" and \"com\" systems. In case of \"com\" systems, existing solvent configurations will be stripped before resolvating. \n\n"
                                fi
                                
				printf "\n\nBuilding topology+parameter files for \"aq\" systems. Water and Ions being added...\n\n"


				for s in aq com; do
					if [ "${s}" == "aq" ]; then boxbuf=${boxbufaq}; else boxbuf=${boxbufcom}; fi
					create_box_twostate "${pff}" "${lff}" "${wm}" "${boxbuild}" "${mdboxshape}" "${boxbuf}" "${ionconc}" "molname-ligname.mapping" "${s}" "${translist[@]}"
                                	printf "\n\nDone...\n\n"
				done


                        fi
		fi			

	cd ${path}
fi
##########################
##########################
##########################



# generation of MD boxes for RSFE calculations
##########################
##########################
##########################
if [ "${ticalc}" == "rsfe" ]; then
        cd $system/setup
                # one-state setup
                if [ "${twostate}" == "false" ]; then
                        printf "\n\n setting up files in \"1-state\" mode... \n\n"
                        if [ "${boxbuild}" == "skip" ]; then
                                flag=0; flag2=0
                                for i in "${!translist[@]}"; do
                                	if [ ! -f "${translist[$i]}_aq.parm7" ] || [ ! -f "${translist[$i]}_aq.rst7" ] || [ ! -f "${translist[$i]}.scmask1" ] || [ ! -f "${translist[$i]}.scmask2" ] || [ ! -f "${translist[$i]}.timask1" ] || [ ! -f "${translist[$i]}.timask2" ]; then
                                	        printf "\n\n One of the following files are missing from $system/setup...\n"
                                	        printf "${translist[$i]}_aq.parm7, ${translist[$i]}_aq.rst7, ${translist[$i]}.scmask1 ${translist[$i]}.scmask2 ${translist[$i]}.timask1 ${translist[$i]}.timask2 \n"
                                	        flag=1
                                	fi
					if [ "${bidirection_aq}" == "true" ]; then
						if [ ! -f "${translistrev[$i]}_aq.parm7" ] || [ ! -f "${translistrev[$i]}_aq.rst7" ] || [ ! -f "${translistrev[$i]}.scmask1" ] || [ ! -f "${translistrev[$i]}.scmask2" ] || [ ! -f "${translistrev[$i]}.timask1" ] || [ ! -f "${translistrev[$i]}.timask2" ]; then
							printf "\n\n One of the following files are missing from $system/setup...\n"
							printf "${translistrev[$i]}_aq.parm7, ${translistrev[$i]}_aq.rst7, ${translistrev[$i]}.scmask1 ${translistrev[$i]}.scmask2 ${translistrev[$i]}.timask1 ${translistrev[$i]}.timask2 \n"
							flag2=1
						fi
					fi

                                done
				if [ "${flag}" -eq 1 ]; then
					printf "\n\n ${boxbuild} was set to 0. Above files must be present in $system/setup folder. Exiting...\n"
					exit 0
				fi
				if [ "${flag2}" -eq 1 ]; then
					printf "\n\n ${boxbuild} was set to 0 and bidirection_aq was set to \"true\". Above files must be present in $system/setup folder. Exiting...\n"
					exit 0
				fi

                        else

                                if [ "${boxbuild}" == "0" ]; then
                                        printf "\n\nFor rsfe calclulations, Boxbuild=0 and Boxbuild=1 are equivalent. water and Ions will be added to \"aq\" systems.\n\n"
                                elif [ "${boxbuild}" == "1" ]; then
                                        printf "\n\nFor rsfe calclulations, Boxbuild=0 and Boxbuild=1 are equivalent. water and Ions will be added to \"aq\" systems.\n\n"
                                elif [ "${boxbuild}" == "2" ]; then
                                        printf "\n\nFor all transformations identical number of water and ions will be added to \"aq\" systems.\n\n"
                                fi

                                printf "\n\nBuilding topology+parameter files for \"aq\" systems. Water and Ions being added...\n\n"
				if [ "${bidirection_aq}" == "true" ]; then
					printf "\n\nbidirection_aq set to true. One state simulations for \"aq\" systems will be setup for both directions... \n\n"
				fi

                                create_box_rsfe "${pff}" "${lff}" "${wm}" "${boxbuild}" "${mdboxshape}" "${boxbufaq}" "${ionconc}" "molname-ligname.mapping" "aq" "${ticalc}" "${bidirection_aq}" "${translist[@]}" 
                                printf "\n\nDone...\n\n"
                        fi
                fi

                # two-state setup
                if [ "${twostate}" == "true" ]; then
                        printf "\n\n setting up files in \"2-state\" mode... \n\n"
                        if [ "${boxbuild}" == "skip" ]; then
                                flag=0; flag2=0
                                for i in "${!translist[@]}"; do
                                        if [ ! -f "${translist[$i]}-1_aq.parm7" ] || [ ! -f "${translist[$i]}-1_aq.rst7" ] || [ ! -f "${translist[$i]}-2_aq.parm7" ] || [ ! -f "${translist[$i]}-2_aq.rst7" ] || [ ! -f "${translist[$i]}.scmask1" ] || [ ! -f "${translist[$i]}.scmask2" ] || [ ! -f "${translist[$i]}.timask1" ] || [ ! -f "${translist[$i]}.timask2" ]; then
                                                printf "\n\n One of the following files are missing from $system/setup...\n"
                                                printf "${translist[$i]}-1_aq.parm7, ${translist[$i]}-1_aq.rst7, ${translist[$i]}-2_aq.parm7, ${translist[$i]}-2_aq.rst7, ${translist[$i]}.scmask1 ${translist[$i]}.scmask2 ${translist[$i]}.timask1 ${translist[$i]}.timask2 \n"
                                                flag=1
                                        fi

                                        if [ "${bidirection_aq}" == "true" ]; then
                                                printf "\n\n bidirection_aq=true will be ignored in \"2-state\" setup of \"com\" systems... \n\n"
                                        fi
                                done
                                if [ "${flag}" -eq 1 ]; then
                                        printf "\n\n ${boxbuild} was set to 0. Above files must be present in $system/setup folder. Exiting...\n"
                                        exit 0
                                fi
                        else

                                if [ "${boxbuild}" == "0" ]; then
                                        printf "\n\nFor rsfe calclulations, Boxbuild=0 and Boxbuild=1 are equivalent. water and Ions will be added to \"aq\" systems.\n\n"
                                elif [ "${boxbuild}" == "1" ]; then
                                        printf "\n\nFor rsfe calclulations, Boxbuild=0 and Boxbuild=1 are equivalent. water and Ions will be added to \"aq\" systems.\n\n"
                                elif [ "${boxbuild}" == "2" ]; then
                                        printf "\n\nFor all transformations identical number of water and ions will be added to \"aq\" systems.\n\n"
                                fi

                                printf "\n\nBuilding topology+parameter files for \"aq\" systems. Water and Ions being added...\n\n"

                                create_box_twostate "${pff}" "${lff}" "${wm}" "${boxbuild}" "${mdboxshape}" "${boxbufaq}" "${ionconc}" "molname-ligname.mapping" "aq" "${translist[@]}"
                                printf "\n\nDone...\n\n"

                        fi
		fi
		
        cd ${path}
fi
##########################
##########################
##########################


# generation of MD boxes for ASFE calculations
##########################
##########################
##########################
if [ "${ticalc}" == "asfe" ]; then
        cd $system/setup
        	printf "\n\n setting up files in \"1-state\" mode... \n\n"
        	if [ "${boxbuild}" == "skip" ]; then
        	        flag=0
        	        for i in "${!translist[@]}"; do
        	                if [ ! -f "${translist[$i]}_aq.parm7" ] || [ ! -f "${translist[$i]}_aq.rst7" ] || [ ! -f "${translist[$i]}.scmask1" ] || [ ! -f "${translist[$i]}.scmask2" ] || [ ! -f "${translist[$i]}.timask1" ] || [ ! -f "${translist[$i]}.timask2" ]; then
        	                        printf "\n\n One of the following files are missing from $system/setup...\n"
        	                        printf "${translist[$i]}_aq.parm7, ${translist[$i]}_aq.rst7, ${translist[$i]}.scmask1 ${translist[$i]}.scmask2 ${translist[$i]}.timask1 ${translist[$i]}.timask2 \n"
        	                        flag=1
        	                fi
        	        done
        	        if [ "${flag}" -eq 1 ]; then exit 0; fi
        	else

        	        if [ "${boxbuild}" == "0" ]; then
        	                printf "\n\nFor asfe calclulations, Boxbuild=0 and Boxbuild=1 are equivalent. water and Ions will be added to \"aq\" systems.\n\n"
        	        elif [ "${boxbuild}" == "1" ]; then
        	                printf "\n\nFor asfe calclulations, Boxbuild=0 and Boxbuild=1 are equivalent. water and Ions will be added to \"aq\" systems.\n\n"
        	        elif [ "${boxbuild}" == "2" ]; then
        	                printf "\n\nFor all transformations identical number of water and ions will be added to \"aq\" systems.\n\n"
        	        fi

        	        printf "\n\nBuilding topology+parameter files for \"aq\" systems. Water and Ions being added...\n\n"
        	        create_box_asfe "${pff}" "${lff}" "${wm}" "${boxbuild}" "${mdboxshape}" "${boxbufaq}" "${ionconc}" "molname-ligname.mapping" "aq" "${ticalc}" "${translist[@]}"
        	        printf "\n\nDone...\n\n"
        	fi

        cd ${path}
fi
##########################
##########################
##########################

