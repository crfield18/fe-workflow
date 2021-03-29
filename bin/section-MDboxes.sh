# generation of MD boxes
##########################
##########################
##########################
if [ "${ticalc}" == "rbfe" ]; then
        if [ "${boxbuild}" != "skip" ]; then
                cd $system/setup
                        if [ "${twostate}" == true ]; then
                                printf "\n\ntwostate=true. The flag \"boxbuild\" is ignored. Merged boxes of endstates will be built containing identican number of water and ions\n\n"
                                printf "\n\nBuilding boxes for only aq systems...\n\n"
                                create_box "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                                printf "\n\nbox building for aq systems complete\n\n"

                        else
                                if [ "${boxbuild}" -eq 0 ]; then
                                        printf "\n\nBuilding boxes for only aq systems...\n\n"
                                        create_box "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                                        printf "\n\nbox building for aq systems complete\n\n"
                                elif [ "${boxbuild}" -eq 1 ]; then
                                        printf "\n\nBuilding boxes for both aq and com systems...\n\n"
                                        create_box "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                                        printf "\n\nbox building for aq systems complete\n\n"
                                        create_box "com" "${boxbufcom}" "${ionconc}" "${comlist[@]}"
                                        printf "\n\nbox building for com systems complete\n\n"
                                elif [ "${boxbuild}" -eq 2 ]; then
                                        printf "\n\nBuilding boxes with same numbers of solvent molecules for aq and com systems...\n\n"
                                        create_box_equal "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                                        printf "\n\nbox building for aq systems complete\n\n"
                                        create_box_equal "com" "${boxbufcom}" "${ionconc}" "${comlist[@]}"
                                        printf "\n\nbox building for com systems complete\n\n"
                                fi
                        fi
                cd $path
        else
                printf "\n\nSkipping box building completely\n\n"
        fi
else
        if [ "${boxbuild}" != "skip" ]; then
                cd $system/setup
                        printf "\n\nBuilding boxes with same numbers of solvent molecules for aq systems...\n\n"
                        create_box_equal "aq" "${boxbufaq}" "${ionconc}" "${liglist[@]}"
                        printf "\n\nbox building for aq systems complete\n\n"
                cd $path
        fi
fi
##########################
##########################
##########################

