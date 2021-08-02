# Organization of initial files organization
# RBFE calculations
if [ "${ticalc}" == "rbfe" ]; then
        cd ${system}/setup
                if [[ -f molname-ligname.mapping && $(cat molname-ligname.mapping | wc -l) -eq "${#uniqueligs[@]}"  ]]; then
                        printf "\n\n\"molname-ligname.mapping\" present in working directory. Skipping initial organization of input files."
                        printf "To enforce this step, delete \"molname-ligname.mapping\" from current directory and re-run script\n\n"
		else
                        preparePDBs "${#uniqueligs[@]}" "${path_to_input}/${system}" "${boxbuild}" "${uniqueligs[@]}"


                       #printf "\n\nBuilding boxes for both aq and com systems...\n\n"
                       #create_box "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                       #printf "\n\nbox building for aq systems complete\n\n"
                       #create_box "com" "${boxbufcom}" "${ionconc}" "${comlist[@]}"
                       #printf "\n\nbox building for com systems complete\n\n"
                       #if [ "${boxbuild}" -ne 2 ]; then boxbuild="skip"; fi
                fi
        cd ${path}
fi

# Organization of initial files organization
# RSFE or ASFE calculations
if [ "${ticalc}" == "rsfe" ] || [ "${ticalc}" == "asfe" ]; then
        cd ${system}/setup
                if [[ -f molname-ligname.mapping && $(cat molname-ligname.mapping | wc -l) -eq "${#uniqueligs[@]}"  ]]; then
                        printf "\n\n\"molname-ligname.mapping\" present in working directory. Skipping initial organization of input files."
                        printf "To enforce this step, delete \"molname-ligname.mapping\" from current directory and re-run script\n\n"
		else
                        truncate -s0 molname-ligname.mapping
			for i in "${!uniqueligs[@]}";do
                		molname=${uniqueligs[$i]}
                		ligname=$(printf "L%02d\n" $i)

				cp ${path_to_input}/${system}/${molname}_0.mol2 .
                		cp ${path_to_input}/${system}/${molname}_0.frcmod .
                		cp ${path_to_input}/${system}/${molname}_0.lib .
                		sed -i "s/LIG/${ligname}/g" ${molname}_0.mol2 ${molname}_0.frcmod ${molname}_0.lib

                		printf "${molname} ${ligname}\n" >> molname-ligname.mapping
			done
		fi
        cd ${path}
fi

