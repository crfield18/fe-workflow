# Organization of initial files organization
# RBFE calculations
if [ "${ticalc}" == "rbfe" ]; then
        cd ${system}/setup
                if [[ -f molname-ligname.mapping && $(cat molname-ligname.mapping | wc -l) -eq "${#uniqueligs[@]}"  ]]; then
                        printf "\n\n\"molname-ligname.mapping\" present in working directory. Skipping initial organization of input files."
                        printf "To enforce this step, delete \"molname-ligname.mapping\" from current directory and re-run script\n\n"
		else
                        preparePDBs "${#uniqueligs[@]}" "${path_to_input}/${system}" "${boxbuild}" "${uniqueligs[@]}"
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
                                sed -i "s/${molname}/${ligname}/g" ${molname}_0.lib # RS Edit

                		printf "${molname} ${ligname}\n" >> molname-ligname.mapping
			done
		fi
        cd ${path}
fi

