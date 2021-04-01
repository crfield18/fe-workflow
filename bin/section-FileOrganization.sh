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
				# check if the flag -- "ATOMS_PER_MOLECULE" is present in user-provided parmfile
				# if absent, re-build parm with setbox command to generate the flag
				# the flag is essential for parmutils-timutate to parse parm files
				if ! grep -Fq "ATOMS_PER_MOLECULE" ${molname}_com.parm7;then
					cat << EOFT > tleap.in
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams frcmod.ff14SB

#source leaprc.water.tip4pew
#loadamberparams frcmod.tip4pew
#loadAmberParams frcmod.ionsjc_tip4pew
#loadoff tip4pewbox.off

loadamberparams ${molname}_com.frcmod
loadoff ${molname}_com.lib

m = loadPdbUsingSeq ${molname}_com.pdb { $(cat ${molname}_com.seq) }
setbox m vdw 
saveamberparm m tmp.parm7 tmp.rst7
quit

EOFT
					tleap -s -f tleap.in
					mv tmp.parm7 ${molname}_com.parm7
					mv tmp.rst7  ${molname}_com.rst7
				fi
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
                                if ! grep -Fq "ATOMS_PER_MOLECULE" ${molname}_aq.parm7;then
					printf "\n\nThe flag -- "ATOMS_PER_MOLECULE" is missing in the file ${system}/${molname}.parm7"
					printf "Rebuilding parm using tleap...\n\n"
                                        cat << EOFT > tleap.in
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams frcmod.ff14SB

#source leaprc.water.tip4pew
#loadamberparams frcmod.tip4pew
#loadAmberParams frcmod.ionsjc_tip4pew
#loadoff tip4pewbox.off

loadamberparams ${molname}_aq.frcmod
loadoff ${molname}_aq.lib

m = loadPdbUsingSeq ${molname}_aq.pdb { $(cat ${molname}_aq.seq) }
setbox m vdw
saveamberparm m tmp.parm7 tmp.rst7
quit

EOFT
                                        tleap -s -f tleap.in > output
                                        mv tmp.parm7 ${molname}_aq.parm7
                                        mv tmp.rst7  ${molname}_aq.rst7
                                fi

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

