# Organization of initial files organization
# RBFE calculations
if [ "${ticalc}" == "rbfe" ]; then
        comlist=("${uniqueligs[@]/%/_com_dry}")
        liglist=("${uniqueligs[@]/%/_lig_dry}")
        cd ${system}/setup
                if [[ -f molname-ligname.mapping && $(cat molname-ligname.mapping | wc -l) -eq "${#uniqueligs[@]}"  ]]; then
                        printf "\n\n\"molname-ligname.mapping\" present in working directory. Skipping initial organization of parm/rst/mol2/seq files."
                        printf "To enforce this step, delete \"molname-ligname.mapping\" from current directory and re-run script\n\n"
                elif [ "${inputformat}" == "parm" ]; then
                        truncate -s0 molname-ligname.mapping
                        for i in "${!uniqueligs[@]}";do
                                molname=${uniqueligs[$i]}
                                ligname=$(printf "L%02d\n" $i)
                                parmutils-Organize.py -p ${path_to_input}/${system}/${molname}.parm7 -c ${path_to_input}/${system}/${molname}.rst7 -m ":${molname},L1" -n "${ligname}" -o "${molname}_com" -e1 ":${ligname}" -o1 "${molname}_lig_dry" -e2 '!:WAT,K+,Na+,Cl-' -o2 "${molname}_com_dry"
                                echo "${molname} ${ligname}" >> molname-ligname.mapping
                        done
                elif [ "${inputformat}" == "pdb" ]; then
                        preparePDBs "${#uniqueligs[@]}" "${path_to_input}/${system}" "${uniqueligs[@]}" "${uniquechgs[@]}"
                        printf "\n\nBuilding boxes for both aq and com systems...\n\n"
                        create_box "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                        printf "\n\nbox building for aq systems complete\n\n"
                        create_box "com" "${boxbufcom}" "${ionconc}" "${comlist[@]}"
                        printf "\n\nbox building for com systems complete\n\n"
                        if [ "${boxbuild}" -ne 2 ]; then boxbuild="skip"; fi
                fi
        cd ${path}
fi

# Organization of initial files organization
# RBFE calculations
if [ "${ticalc}" == "rsfe" ]; then
        liglist=("${uniqueligs[@]/%/_lig_dry}")
        cd ${system}/setup
                if [[ -f molname-ligname.mapping && $(cat molname-ligname.mapping | wc -l) -eq "${#uniqueligs[@]}"  ]]; then
                        printf "\n\n\"molname-ligname.mapping\" present in working directory. Skipping initial organization of parm/rst/mol2/seq files."
                        printf "To enforce this step, delete \"molname-ligname.mapping\" from current directory and re-run script\n\n"
                elif [ "${inputformat}" == "parm" ]; then
                        truncate -s0 molname-ligname.mapping
                        for i in "${!uniqueligs[@]}";do
                                molname=${uniqueligs[$i]}
                                ligname=$(printf "L%02d\n" $i)
                                parmutils-Organize.py -p ${path_to_input}/${system}/${molname}.parm7 -c ${path_to_input}/${system}/${molname}.rst7 -m ":${molname},L1" -n "${ligname}" -o "${molname}_aq" -e1 ":${ligname}" -o1 "${molname}_lig_dry"
                                echo "${molname} ${ligname}" >> molname-ligname.mapping
                        done
                elif [ "${inputformat}" == "pdb" ]; then
                        preparePDBs "${#uniqueligs[@]}" "${path_to_input}/${system}" "${uniqueligs[@]}" "${uniquechgs[@]}"
                        printf "\n\nBuilding boxes for both aq and com systems...\n\n"
                        create_box_noseq "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                        printf "\n\nbox building for aq systems complete\n\n"
                        if [ "${boxbuild}" -ne 2 ]; then boxbuild="skip"; fi
                fi
        cd ${path}
fi

